// To do: More formats?
// To do: Make stdin work on MacOS
// To do: Implement multithreading. One thread reads input and moves view etc,
// another draws. Drawings can be "skipped" if there are still queued inputs, perhaps?
// To do: Possibly show deviation from consensus?

mod constants;
mod data;

use constants::HEADER_LINES;
use data::View;

use std::cmp::min;
use std::io::{stdout, BufReader, Write};
use std::path::Path;

use anyhow::Result;

use regex::RegexBuilder;

use unicode_segmentation::UnicodeSegmentation;

use crossterm::{
    cursor,
    event::{self, Event, KeyCode, KeyEvent},
    execute, queue,
    style::{Color, Print, ResetColor, SetBackgroundColor, SetForegroundColor},
    terminal::{self, ClearType},
};

fn move_view_and_redraw<T: Write>(view: &mut View, io: &mut T, dy: isize, dx: isize) -> Result<()> {
    let old_rowstart = view.rowstart;
    let old_colstart = view.colstart;

    view.move_view(dy, dx);

    // Only update the view if the view was actually moved.
    if view.rowstart != old_rowstart {
        draw_names(io, view)?;
        io.flush()?;
    }
    if view.colstart != old_colstart {
        draw_ruler(io, view)?;
        io.flush()?;
    }
    if view.rowstart != old_rowstart || view.colstart != old_colstart {
        draw_sequences(io, view)?;
        io.flush()?;
    }
    Ok(())
}

fn get_color_background_dna(byte: u8) -> Option<Color> {
    match byte {
        b'a' | b'A' => Some(Color::AnsiValue(228)), // yellow
        b'c' | b'C' => Some(Color::AnsiValue(77)),  // green
        b'g' | b'G' => Some(Color::AnsiValue(39)),  // blue
        b't' | b'T' | b'u' | b'U' => Some(Color::AnsiValue(168)), // pink
        _ => None,
    }
}

fn get_color_background_aa(byte: u8) -> Option<Color> {
    match byte {
        // Negative (reds)
        b'e' | b'E' => Some(Color::AnsiValue(198)),
        b'd' | b'D' => Some(Color::AnsiValue(161)),

        // Positive (blues)
        b'r' | b'R' => Some(Color::AnsiValue(39)),
        b'k' | b'K' => Some(Color::AnsiValue(26)),
        b'h' | b'H' => Some(Color::AnsiValue(32)),

        // Aromatic (yellows)
        b'f' | b'F' => Some(Color::AnsiValue(184)),
        b'w' | b'W' => Some(Color::AnsiValue(228)),
        b'y' | b'Y' => Some(Color::AnsiValue(186)),

        // Aliphatic (greys)
        b'a' | b'A' => Some(Color::AnsiValue(244)),
        b'v' | b'V' => Some(Color::AnsiValue(246)),
        b'l' | b'L' => Some(Color::AnsiValue(248)),
        b'i' | b'I' => Some(Color::AnsiValue(250)),
        b'm' | b'M' => Some(Color::AnsiValue(252)),

        // Neutral (greens)
        b's' | b'S' => Some(Color::AnsiValue(40)),
        b't' | b'T' => Some(Color::AnsiValue(42)),
        b'n' | b'N' => Some(Color::AnsiValue(76)),
        b'q' | b'Q' => Some(Color::AnsiValue(78)),

        b'c' | b'C' => Some(Color::AnsiValue(112)),
        b'p' | b'P' => Some(Color::AnsiValue(114)),
        b'g' | b'G' => Some(Color::AnsiValue(149)),
        _ => None,
    }
}

fn display<T: Write>(io: &mut T, view: &mut View) -> Result<()> {
    terminal::enable_raw_mode()?;
    execute!(io, terminal::EnterAlternateScreen, cursor::Hide,)?;
    draw_all(io, view)?;
    io.flush()?;

    loop {
        let event = event::read()?;

        // Break on Q or Esc or Control-C
        if event == Event::Key(KeyCode::Esc.into())
            || event == Event::Key(KeyCode::Char('q').into())
            || event
                == Event::Key(KeyEvent {
                    code: KeyCode::Char('c'),
                    modifiers: event::KeyModifiers::CONTROL,
                })
        {
            break;
        }

        match event {
            Event::Key(kevent) => {
                let delta = match kevent {
                    KeyEvent {
                        code: KeyCode::Left,
                        modifiers: event::KeyModifiers::NONE,
                    } => Some((0, -1)),
                    KeyEvent {
                        code: KeyCode::Right,
                        modifiers: event::KeyModifiers::NONE,
                    } => Some((0, 1)),
                    KeyEvent {
                        code: KeyCode::Down,
                        modifiers: event::KeyModifiers::NONE,
                    } => Some((1, 0)),
                    KeyEvent {
                        code: KeyCode::Up,
                        modifiers: event::KeyModifiers::NONE,
                    } => Some((-1, 0)),

                    // SHIFT: Move by 10
                    KeyEvent {
                        code: KeyCode::Left,
                        modifiers: event::KeyModifiers::SHIFT,
                    } => Some((0, -10)),
                    KeyEvent {
                        code: KeyCode::Right,
                        modifiers: event::KeyModifiers::SHIFT,
                    } => Some((0, 10)),
                    KeyEvent {
                        code: KeyCode::Down,
                        modifiers: event::KeyModifiers::SHIFT,
                    } => Some((10, 0)),
                    KeyEvent {
                        code: KeyCode::Up,
                        modifiers: event::KeyModifiers::SHIFT,
                    } => Some((-10, 0)),

                    // CONTROL: Move to end
                    KeyEvent {
                        code: KeyCode::Left,
                        modifiers: event::KeyModifiers::CONTROL,
                    } => Some((0, isize::MIN)),
                    KeyEvent {
                        code: KeyCode::Right,
                        modifiers: event::KeyModifiers::CONTROL,
                    } => Some((0, isize::MAX)),
                    KeyEvent {
                        code: KeyCode::Down,
                        modifiers: event::KeyModifiers::CONTROL,
                    } => Some((isize::MAX, 0)),
                    KeyEvent {
                        code: KeyCode::Up,
                        modifiers: event::KeyModifiers::CONTROL,
                    } => Some((isize::MIN, 0)),
                    _ => None,
                };
                if let Some((dy, dx)) = delta {
                    move_view_and_redraw(view, io, dy, dx)?;
                };
                let name_move = match kevent {
                    KeyEvent {
                        code: KeyCode::Char(','),
                        modifiers: event::KeyModifiers::NONE,
                    } => Some(-1),
                    KeyEvent {
                        code: KeyCode::Char('.'),
                        modifiers: event::KeyModifiers::NONE,
                    } => Some(1),
                    _ => None,
                };

                // If the names were actually moved, re-draw the screen.
                if let Some(delta) = name_move {
                    let old_namewidth = view.namewidth;
                    view.resize_names(delta);
                    if old_namewidth != view.namewidth {
                        draw_all(io, view)?;
                    }
                }

                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('f'),
                        modifiers: event::KeyModifiers::CONTROL,
                    })
                {
                    enter_search_mode(io, view)?;
                    draw_default_footer(io, view)?;
                    io.flush()?;
                };

                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('j'),
                        modifiers: event::KeyModifiers::CONTROL,
                    })
                {
                    enter_jumpcol_mode(io, view)?;
                    draw_default_footer(io, view)?;
                    io.flush()?;
                };

                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('r'),
                        modifiers: event::KeyModifiers::NONE,
                    })
                {
                    draw_all(io, view)?;
                    io.flush()?;
                };
            }
            Event::Resize(ncols, nrows) => {
                view.resize(ncols, nrows);
                draw_all(io, view)?;
            }
            _ => (),
        };
    }
    Ok(())
}

// TODO: Refactor - this and search mode?
fn enter_jumpcol_mode<T: Write>(io: &mut T, view: &mut View) -> Result<()> {
    let mut query = String::new();
    let mut invalid_column = false;
    loop {
        draw_jump_footer(io, &query, view, invalid_column)?;
        io.flush()?;
        let event = event::read()?;

        if event == Event::Key(KeyCode::Esc.into())
            || event
                == Event::Key(KeyEvent {
                    code: KeyCode::Char('c'),
                    modifiers: event::KeyModifiers::CONTROL,
                })
        {
            break;
        };

        match event {
            Event::Key(KeyEvent {
                code: KeyCode::Char(c),
                modifiers: _,
            }) => {
                if let '0'..='9' = c {
                    query.push(c)
                };
            }
            Event::Key(KeyEvent {
                code: KeyCode::Backspace,
                modifiers: _,
            }) => {
                query.pop();
            }
            Event::Key(KeyEvent {
                code: KeyCode::Enter,
                modifiers: _,
            }) => {
                match query.parse::<usize>() {
                    Ok(n) => {
                        if n < 1 || n > view.ncols() {
                            invalid_column = true;
                        } else {
                            // minus one because our alignment should be 1-indexed
                            move_view_and_redraw(
                                view,
                                io,
                                0,
                                n as isize - view.colstart as isize - 1,
                            )?;
                            break;
                        }
                    }
                    Err(_) => {
                        unreachable!(); // should never happen, since we only accept numerical inputs
                    }
                }
            }
            _ => (),
        }
    }
    Ok(())
}

fn draw_jump_footer<T: Write>(
    io: &mut T,
    query: &str,
    view: &View,
    invalid_column: bool,
) -> Result<()> {
    let mut text = "[Esc: Quit] ".to_owned();
    text.push_str(if invalid_column {
        "Invalid number: "
    } else {
        "Jump to column: "
    });
    let background_color = if invalid_column {
        Color::Red
    } else {
        Color::Grey
    };
    text.push_str(query);
    draw_footer(io, view, &text, background_color)?;
    Ok(())
}

fn enter_search_mode<T: Write>(io: &mut T, view: &mut View) -> Result<()> {
    let mut query = String::new();
    let mut last_result = None; // we change the prompt based on the last result
    loop {
        draw_search_footer(io, view, &query, last_result)?;
        io.flush()?;
        let event = event::read()?;

        // Quit on escape or Ctrl-C
        if event == Event::Key(KeyCode::Esc.into())
            || event
                == Event::Key(KeyEvent {
                    code: KeyCode::Char('c'),
                    modifiers: event::KeyModifiers::CONTROL,
                })
        {
            break;
        };

        match event {
            Event::Key(KeyEvent {
                code: KeyCode::Char(c),
                modifiers: _,
            }) => {
                query.push(c);
            }
            Event::Key(KeyEvent {
                code: KeyCode::Backspace,
                modifiers: _,
            }) => {
                query.pop();
            }
            Event::Key(KeyEvent {
                code: KeyCode::Enter,
                modifiers: _,
            }) => {
                let search_result = search_query(view, &query);
                match search_result {
                    // Invalid result: Just continue
                    SearchResult::RegexErr | SearchResult::NoMatch => {
                        last_result = Some(search_result);
                        continue;
                    }
                    // Else, go to the correct row.
                    SearchResult::MatchHeader(nrow) => {
                        let dy = nrow as isize
                            - (view.rowstart + (view.seq_nrows_display() / 2)) as isize;
                        move_view_and_redraw(view, io, dy, 0)?;
                        let header_row = (nrow.saturating_sub(view.rowstart) + HEADER_LINES) as u16;
                        draw_highlight(
                            io,
                            header_row,
                            0,
                            view.namewidth,
                            &view.graphemes()[nrow].string,
                        )?;
                        return Ok(());
                    }
                    // Else go to correct seq.
                    SearchResult::MatchSeq { row, start, stop } => {
                        let dy = row as isize
                            - (view.rowstart + (view.seq_nrows_display() / 2)) as isize;
                        let dx = start as isize
                            - (view.colstart + (view.seq_ncols_display() / 2)) as isize;
                        move_view_and_redraw(view, io, dy, dx)?;
                        let seq_row = (row.saturating_sub(view.rowstart) + HEADER_LINES) as u16;
                        let seq_col =
                            start.saturating_sub(view.colstart) as u16 + (view.namewidth + 1);
                        let highlight_str = unsafe {
                            std::str::from_utf8_unchecked(&view.seqs()[row][start..stop])
                        };
                        draw_highlight(io, seq_row, seq_col, view.term_ncols, highlight_str)?;
                        return Ok(());
                    }
                }
            }
            _ => (),
        }
    }
    Ok(())
}

#[derive(Clone, Copy)]
enum SearchResult {
    MatchSeq {
        row: usize,
        start: usize,
        stop: usize,
    },
    MatchHeader(usize),
    RegexErr,
    NoMatch,
}

fn search_query(view: &View, query: &str) -> SearchResult {
    let re = match RegexBuilder::new(query).case_insensitive(true).build() {
        Ok(regex) => regex,
        Err(_) => {
            return SearchResult::RegexErr; // TODO: Somehow print that it didn't work
        }
    };

    // First search the headers
    for (n_header, header) in view.graphemes().iter().enumerate() {
        if re.is_match(&header.string) {
            return SearchResult::MatchHeader(n_header);
        }
    }

    // Then search the sequences
    for (rowno, seq) in view.seqs().iter().enumerate() {
        // We have checked on instantiation that this is OK, and do not provide
        // any functionality to mutate the sequences
        let string = unsafe { std::str::from_utf8_unchecked(seq) };
        if let Some(regex_match) = re.find(string) {
            return SearchResult::MatchSeq {
                row: rowno,
                start: regex_match.start(),
                stop: regex_match.end(),
            };
        }
    }
    SearchResult::NoMatch
}

fn draw_search_footer<T: Write>(
    io: &mut T,
    view: &View,
    query: &str,
    last_result: Option<SearchResult>,
) -> Result<()> {
    let mut footer = String::from("[Esc: Quit]");
    let (background_color, message) = match last_result {
        None => (Color::Grey, " Enter query | "),
        Some(SearchResult::NoMatch) => (Color::Red, " Not found | "),
        Some(SearchResult::RegexErr) => (Color::Red, " Bad regex | "),
        _ => unreachable!(), // The last two are hits - they should never happen!
    };
    footer.push_str(message);
    footer.push_str(query);

    draw_footer(io, view, &footer, background_color)
}

fn draw_all<T: Write>(io: &mut T, view: &View) -> Result<()> {
    execute!(io, ResetColor, terminal::Clear(ClearType::All),)?;

    if view.term_ncols < 2 || view.term_nrows < 4 {
        // This seems silly, but I have it because it allows me to assume a minimal
        // terminal size when drawing the regular alignment
        execute!(io, cursor::MoveTo(0, 0), Print(":("),)?;
    } else {
        draw_ruler(io, view)?;
        draw_names(io, view)?;
        draw_default_footer(io, view)?;
        draw_sequences(io, view)?;
    }
    Ok(io.flush()?)
}

fn draw_names<T: Write>(io: &mut T, view: &View) -> Result<()> {
    if let Some(range) = view.seq_row_range() {
        for (i, nameindex) in range.into_iter().enumerate() {
            let termrow = (i + HEADER_LINES) as u16;
            let name = if view.namewidth == 0 {
                "".to_owned()
            } else {
                let graphemes = &view.graphemes()[nameindex];
                let elide = graphemes.len() > view.namewidth as usize;
                let namelen = min(graphemes.len(), view.namewidth as usize - elide as usize);
                let mut name = graphemes.get_n_graphemes(namelen).unwrap().to_owned();
                if elide {
                    name.push('…')
                } else {
                    let missing_graphemes = view.namewidth as usize - graphemes.len();
                    name.push_str(&" ".repeat(missing_graphemes));
                }
                name
            };
            queue!(
                io,
                cursor::MoveTo(0, termrow),
                Print(name),
                cursor::MoveTo(view.namewidth, termrow),
                Print('│'),
            )?;
        }
    }
    Ok(())
}

fn draw_ruler<T: Write>(io: &mut T, view: &View) -> Result<()> {
    // Make the top line with the numbers.
    let num_string = {
        let mut string = String::new();
        let start = 10 * (view.colstart / 10);
        let stop = view.colstart + view.seq_ncols_display() + 1;
        for i in (start..stop).step_by(10) {
            let add = i.to_string();
            string.push_str(&add);
            string.push_str(&" ".repeat(10 - add.len()));
        }
        let start_index = view.colstart - start;
        string[start_index..=(start_index + view.seq_ncols_display())].to_owned()
    };

    // Make the bottom line with the ticks
    let tick_string = {
        let aln_range = view.colstart..=(view.colstart + view.seq_ncols_display() - 1);
        let mut tick_string = "┌".to_owned();
        for alncol in aln_range {
            tick_string.push(if (alncol + 1) % 10 == 0 { '┴' } else { '─' })
        }
        tick_string
    };

    Ok(queue!(
        io,
        cursor::MoveTo(view.namewidth, 0),
        Print(num_string),
        cursor::MoveTo(view.namewidth, 1),
        Print(tick_string),
    )?)
}

fn draw_default_footer<T: Write>(io: &mut T, view: &View) -> Result<()> {
    draw_footer(
        io,
        view,
        "q/Esc: Quit   ←/→/↑/↓ + None/Shift/Ctrl: Move alignment   ./,: Adjust names   Ctrl+f: Find    Ctrl+j Jump",
        Color::Grey
    )
}

fn draw_footer<T: Write>(io: &mut T, view: &View, text: &str, background: Color) -> Result<()> {
    // First we create the full footer, then we truncate, if needed
    let mut footer = text.to_owned();

    // Pad or truncate footer to match num columns
    let nchars = footer.chars().count();
    let ncols = view.term_ncols as usize;

    if nchars > ncols {
        footer = UnicodeSegmentation::graphemes(footer.as_str(), true)
            .take(ncols)
            .collect::<String>();
    } else {
        footer.push_str(" ".repeat(ncols - nchars).as_str())
    }

    Ok(queue!(
        io,
        SetBackgroundColor(background),
        SetForegroundColor(Color::Black),
        cursor::MoveTo(0, (view.term_nrows - 1) as u16),
        Print(footer),
        ResetColor,
    )?)
}

// This function is the performance bottleneck, so I try to queue as little as
// possible to the terminal by checking every single operation to see if it's
// necessary
fn draw_sequences<T: Write>(io: &mut T, view: &View) -> Result<()> {
    let row_range = match view.seq_row_range() {
        Some(n) => n,
        None => return Ok(()),
    };
    let col_range = match view.seq_col_range() {
        Some(n) => n,
        None => return Ok(()),
    };

    let mut oldcolor: Option<Color> = None;
    let mut is_foreground_black = false;
    for (i, alnrow) in row_range.enumerate() {
        let termrow = (i + HEADER_LINES) as u16;
        queue!(io, cursor::MoveTo(view.namewidth + 1, termrow),)?;
        for byte in view.seqs()[alnrow][col_range.clone()].iter() {
            let color = if view.is_aa() {
                get_color_background_aa(*byte)
            } else {
                get_color_background_dna(*byte)
            };
            // if color is same as before, don't queue up color changing operations
            if color != oldcolor {
                // If the current cell needs a background color, set it
                if let Some(clr) = color {
                    // only set foreground to black if it isn't already
                    if !is_foreground_black {
                        queue!(io, SetForegroundColor(Color::Black))?;
                        is_foreground_black = true;
                    };
                    queue!(io, SetBackgroundColor(clr))?;
                // Otherwise, reset to default color scheme
                } else {
                    queue!(io, ResetColor)?;
                    is_foreground_black = false;
                }
            };
            queue!(io, Print(*byte as char))?;
            oldcolor = color;
        }
    }

    Ok(queue!(io, ResetColor,)?)
}

fn draw_highlight<T: Write>(
    io: &mut T,
    screenrow: u16,
    screencol: u16,
    maxcols: u16,
    text: &str,
) -> Result<()> {
    let max_len = (maxcols - screencol) as usize;
    let truncated = if text.len() > max_len {
        if text.is_ascii() {
            &text[0..max_len]
        } else {
            let lastindex = if let Some((index, _)) =
                UnicodeSegmentation::grapheme_indices(text, true).nth(max_len)
            {
                index - 1
            } else {
                text.len()
            };
            &text[0..=lastindex]
        }
    } else {
        text
    };
    Ok(queue!(
        io,
        SetBackgroundColor(Color::White),
        SetForegroundColor(Color::Black),
        cursor::MoveTo(screencol, screenrow),
        Print(truncated),
        ResetColor,
    )?)
}

fn clean_terminal<T: Write>(io: &mut T) -> Result<()> {
    execute!(io, ResetColor, cursor::Show, terminal::LeaveAlternateScreen)?;
    Ok(terminal::disable_raw_mode()?)
}

fn main() {
    let args = clap::App::new("Alen")
        .version("0.1")
        .author("Jakob Nybo Nissen <jakobnybonissen@gmail.com>")
        .about("Simple terminal alignment viewer")
        .arg(
            clap::Arg::with_name("alignment")
                .help("Input alignment in FASTA format")
                .takes_value(true)
                .required(true),
        )
        .arg(
            clap::Arg::with_name("uppercase")
                .short("u")
                .takes_value(false)
                .help("Displays sequences in uppercase"),
        )
        .arg(
            clap::Arg::with_name("aminoacids")
                .short("a")
                .takes_value(false)
                .help("Force parsing as amino acids"),
        )
        .get_matches();

    let filename = args.value_of("alignment").unwrap();

    // Check if file exists
    if filename != "-" && !Path::new(filename).is_file() {
        println!("ERROR: Filename not found: \"{}\"", filename);
        std::process::exit(1);
    }
    let file = std::fs::File::open(filename).unwrap_or_else(|err| {
        println!("ERROR reading file: {}", err);
        std::process::exit(1)
    });

    let buffered_io = BufReader::new(file);
    let uppercase = args.is_present("uppercase");
    let must_aa = args.is_present("aminoacids");
    let view = View::from_reader(BufReader::new(buffered_io), uppercase, must_aa);
    match view {
        Err(e) => {
            println!("ERROR when loading FASTA: {}", e);
            std::process::exit(1);
        }
        Ok(mut view) => {
            let mut io = stdout();
            if let Err(e) = display(&mut io, &mut view) {
                println!("Error: {}", e);
            }
            clean_terminal(&mut io).unwrap();
            std::process::exit(0);
        }
    }
}
