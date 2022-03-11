// To do: More formats?
// To do: Make stdin work on MacOS
// To do: Implement multithreading. One thread reads input and moves view etc,
// another draws. Drawings can be "skipped" if there are still queued inputs, perhaps?

mod constants;
mod data;

use constants::HEADER_LINES;
use data::{Graphemes, View};

use std::cmp::min;
use std::ffi::OsString;
use std::io::{stdout, BufReader, Write};
use std::path::Path;

use anyhow::Result;

use clap::Parser;

use regex::RegexBuilder;

use unicode_segmentation::UnicodeSegmentation;

use crossterm::{
    cursor,
    event::{self, Event, KeyCode, KeyEvent},
    execute, queue,
    style::{Color, Print, ResetColor, SetBackgroundColor, SetForegroundColor},
    terminal::{self, ClearType},
};

// This struct represents the mutable terminal state
// It's important to cache the current colors, because the major bottleneck
// in this app is switching between colors, so we keep that to a minimum.
struct TerminalIO<T> {
    io: T,
    // If has_color is true, color is always None
    color: Option<Color>,
    has_color: bool,
}

fn set_terminal_color<T: Write>(io: &mut TerminalIO<T>, color: Option<Color>) -> Result<()> {
    if color != io.color && io.has_color {
        if let Some(clr) = color {
            queue!(io.io, SetBackgroundColor(clr))?;
            if color.is_some() {
                queue!(io.io, SetForegroundColor(Color::Black))?;
            }
        } else {
            queue!(io.io, ResetColor)?;
        }
        io.color = color;
    }
    Ok(())
}

// ================================= Footers =================================
fn draw_footer_text<T: Write>(
    io: &mut TerminalIO<T>,
    view: &View,
    text: &str,
    background: Color,
) -> Result<()> {
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
    set_terminal_color(io, Some(background))?;

    Ok(queue!(
        io.io,
        cursor::MoveTo(0, (view.term_nrows - 1) as u16),
        Print(footer),
    )?)
}

fn draw_default_footer<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    draw_footer_text(
        io,
        view,
        "q/Esc: Quit | [^⇧] + ←/→/↑/↓: Move | ./,: Adjust names | ^f: Find | ^j: Jump | ^s: Select | c: Consensus | r: Redraw",
        Color::Grey
    )
}

fn draw_select_footer<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    draw_footer_text(
        io,
        view,
        "q/Esc: Quit | [^⇧] + ↑/↓: Select | [⇧] + k/j: Move up/down | t/b: To top/bottom",
        Color::Grey,
    )
}

fn draw_search_footer<T: Write>(
    io: &mut TerminalIO<T>,
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
    draw_footer_text(io, view, &footer, background_color)
}

fn draw_jump_footer<T: Write>(
    io: &mut TerminalIO<T>,
    query: &str,
    view: &View,
    invalid_column: bool,
) -> Result<()> {
    let mut text = "[Esc: Quit] ".to_owned();
    text.push_str(if invalid_column {
        "Invalid column: "
    } else {
        "Jump to column: "
    });
    let background_color = if invalid_column {
        Color::Red
    } else {
        Color::Grey
    };
    text.push_str(query);
    draw_footer_text(io, view, &text, background_color)?;
    Ok(())
}

// ================================= Ruler =================================

fn draw_ruler<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    set_terminal_color(io, None)?;

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
        let mut tick_string = if view.colstart % 10 == 0 {
            "├"
        } else {
            "┌"
        }
        .to_owned();
        for alncol in aln_range {
            tick_string.push(if (alncol + 1) % 10 == 0 { '┴' } else { '─' })
        }
        tick_string
    };

    Ok(queue!(
        io.io,
        cursor::MoveTo(view.namewidth, 0),
        Print(num_string),
        cursor::MoveTo(view.namewidth, 1),
        Print(tick_string),
    )?)
}

// ================================= Names =================================

fn draw_names<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    set_terminal_color(io, None)?;
    if view.consensus {
        draw_consensus_names(io, view)?;
    } else {
        draw_nonconsensus_names(io, view)?;
    }
    Ok(())
}

fn draw_nonconsensus_names<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    if let Some(range) = view.seq_row_range() {
        for (i, nameindex) in range.into_iter().enumerate() {
            let termrow = (i + HEADER_LINES) as u16;
            draw_name(io, view.namewidth, &view.graphemes()[nameindex], termrow)?;
        }
    }
    Ok(())
}

fn draw_consensus_names<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    if view.term_nrows > HEADER_LINES as u16 {
        draw_name(
            io,
            view.namewidth,
            &Graphemes::new("consensus"),
            HEADER_LINES as u16,
        )?;
    }
    if let Some(range) = view.seq_row_range() {
        for (i, nameindex) in range.enumerate() {
            let termrow = (i + 1 + HEADER_LINES) as u16;
            draw_name(io, view.namewidth, &view.graphemes()[nameindex], termrow)?;
        }
    }
    Ok(())
}

fn draw_name<T: Write>(
    io: &mut TerminalIO<T>,
    namewidth: u16,
    graphemes: &Graphemes,
    termrow: u16,
) -> Result<()> {
    let name = if namewidth == 0 {
        "│".to_owned()
    } else {
        let elide = graphemes.len() > namewidth as usize;
        let namelen = min(graphemes.len(), namewidth as usize - elide as usize);
        let mut name = graphemes.get_n_graphemes(namelen).unwrap().to_owned();
        if elide {
            name.push('…')
        } else {
            let missing_graphemes = namewidth as usize - graphemes.len();
            name.push_str(&" ".repeat(missing_graphemes));
        }
        name.push('│');
        name
    };
    queue!(io.io, cursor::MoveTo(0, termrow), Print(name),)?;
    Ok(())
}

fn highlight_name<T: Write>(io: &mut TerminalIO<T>, view: &mut View, nrow: usize) -> Result<()> {
    let header_row = view.display_row_of_index(nrow).unwrap();
    draw_highlight(
        io,
        header_row,
        0,
        view.namewidth,
        &view.graphemes()[nrow].string,
    )
}

// ================================= Names =================================

fn draw_sequences<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    if view.consensus {
        draw_consensus_sequences(io, view)?;
    } else {
        draw_nonconsensus_sequences(io, view)?;
    }
    Ok(())
}

fn draw_nonconsensus_sequences<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    let row_range = match view.seq_row_range() {
        Some(n) => n,
        None => return Ok(()),
    };
    let col_range = match view.seq_col_range() {
        Some(n) => n,
        None => return Ok(()),
    };

    for (i, alnrow) in row_range.enumerate() {
        let termrow = (i + HEADER_LINES) as u16;
        let seq = &view.seqs()[alnrow][col_range.clone()];
        draw_sequence(io, view.namewidth + 1, view.is_aa(), seq, termrow)?;
    }
    Ok(())
}

fn draw_sequence<T: Write>(
    io: &mut TerminalIO<T>,
    colstart: u16,
    is_aa: bool,
    seq: &[u8],
    termrow: u16,
) -> Result<()> {
    queue!(io.io, cursor::MoveTo(colstart, termrow))?;
    for byte in seq {
        let color = if is_aa {
            get_color_background_aa(*byte)
        } else {
            get_color_background_dna(*byte)
        };
        set_terminal_color(io, color)?;
        queue!(io.io, Print(*byte as char))?;
    }
    Ok(())
}

fn draw_consensus_sequences<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    let col_range = match view.seq_col_range() {
        Some(n) => n,
        None => return Ok(()),
    };

    // First draw top row
    let x = view.consensus().unwrap();
    let cons_seq = &x[col_range.clone()];
    //let cons_seq = &view.consensus().unwrap().as_ref()[col_range.clone()];
    draw_top_consensus(io, view.namewidth + 1, view.is_aa(), cons_seq)?;

    // Then draw rest, if applicable
    if let Some(alnrows) = view.seq_row_range() {
        for (i, alnrow) in alnrows.enumerate() {
            let termrow = (i + HEADER_LINES + 1) as u16;
            let seq = &view.seqs()[alnrow][col_range.clone()];
            draw_consensus_other_seq(io, view.namewidth + 1, termrow, view.is_aa(), seq, cons_seq)?
        }
    }
    Ok(())
}

fn draw_top_consensus<T: Write>(
    io: &mut TerminalIO<T>,
    colstart: u16,
    is_aa: bool,
    seq: &[Option<u8>],
) -> Result<()> {
    queue!(io.io, cursor::MoveTo(colstart, HEADER_LINES as u16),)?;
    for maybe_base in seq {
        let (background_color, symbol) = if let Some(byte) = maybe_base {
            let bc = if is_aa {
                get_color_background_aa(*byte)
            } else {
                get_color_background_dna(*byte)
            };
            (bc, *byte as char)
        } else {
            (None, ' ')
        };
        set_terminal_color(io, background_color)?;
        queue!(io.io, Print(symbol))?;
    }
    Ok(())
}

fn draw_consensus_other_seq<T: Write>(
    io: &mut TerminalIO<T>,
    colstart: u16,
    termrow: u16,
    is_aa: bool,
    seq: &[u8],
    cons: &[Option<u8>],
) -> Result<()> {
    queue!(io.io, cursor::MoveTo(colstart, termrow))?;
    for (byte, maybe_cons) in seq.iter().zip(cons.iter()) {
        let (color, symbol) =
            if maybe_cons.is_some() && maybe_cons.unwrap() & 0b11011111 == byte & 0b11011111 {
                (None, ' ')
            } else {
                let color = if is_aa {
                    get_color_background_aa(*byte)
                } else {
                    get_color_background_dna(*byte)
                };
                (color, *byte as char)
            };
        set_terminal_color(io, color)?;
        queue!(io.io, Print(symbol))?;
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
        b'k' | b'K' => Some(Color::AnsiValue(27)),
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

// =========================================== Misc ===========================

fn draw_highlight<T: Write>(
    io: &mut TerminalIO<T>,
    screenrow: u16,
    screencol: u16,
    maxcols: u16,
    text: &str,
) -> Result<()> {
    let max_len = (maxcols - screencol) as usize;
    if max_len == 0 {
        return Ok(());
    }

    // Truncate to max_len
    let truncated = if text.len() > max_len {
        let mut before_dots = if text.is_ascii() {
            text[0..max_len - 1].to_owned()
        } else {
            let lastindex = if let Some((index, _)) =
                UnicodeSegmentation::grapheme_indices(text, true).nth(max_len - 1)
            {
                index - 1
            } else {
                text.len()
            };
            text[0..=lastindex].to_owned()
        };
        before_dots.push('…');
        before_dots
    } else {
        text.to_owned()
    };
    set_terminal_color(io, Some(Color::White))?;
    Ok(queue!(
        io.io,
        cursor::MoveTo(screencol, screenrow),
        Print(truncated),
    )?)
}

fn move_view_and_redraw<T: Write>(
    io: &mut TerminalIO<T>,
    view: &mut View,
    dy: isize,
    dx: isize,
) -> Result<()> {
    let old_rowstart = view.rowstart;
    let old_colstart = view.colstart;

    view.move_view(dy, dx);

    // Only update the view if the view was actually moved.
    if view.rowstart != old_rowstart {
        draw_names(io, view)?;
    }
    if view.colstart != old_colstart {
        draw_ruler(io, view)?;
    }
    if view.rowstart != old_rowstart || view.colstart != old_colstart {
        draw_sequences(io, view)?;
        io.io.flush()?;
    }
    Ok(())
}

// ================================================== Loops =======================

fn default_loop<T: Write>(io: &mut TerminalIO<T>, view: &mut View) -> Result<()> {
    draw_default_mode_screen(io, view)?;
    loop {
        let event = event::read()?;

        // Break on Q or Esc or Control-C
        if event == Event::Key(KeyCode::Esc.into())
            || event == Event::Key(KeyCode::Char('q').into())
            || event == Event::Key(KeyCode::Char('Q').into())
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
                // Arrow keys: Move
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
                    move_view_and_redraw(io, view, dy, dx)?;
                    continue;
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
                        draw_default_mode_screen(io, view)?;
                        continue;
                    }
                }

                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('f'),
                        modifiers: event::KeyModifiers::CONTROL,
                    })
                    || kevent
                        == (KeyEvent {
                            code: KeyCode::Char('F'),
                            modifiers: event::KeyModifiers::CONTROL,
                        })
                {
                    search_loop(io, view)?;
                    // Return from search mode
                    draw_default_mode_screen(io, view)?;
                    continue;
                };

                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('s'),
                        modifiers: event::KeyModifiers::CONTROL,
                    })
                    || kevent
                        == (KeyEvent {
                            code: KeyCode::Char('S'),
                            modifiers: event::KeyModifiers::CONTROL,
                        })
                {
                    // We can't enter select mode if there are not enough screen space to
                    // show any rows.
                    if let Some(range) = view.seq_row_range() {
                        select_loop(io, view, *range.start())?;
                        // Return from select mode
                        draw_default_mode_screen(io, view)?;
                    }
                    continue;
                };

                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('j'),
                        modifiers: event::KeyModifiers::CONTROL,
                    })
                    || kevent
                        == (KeyEvent {
                            code: KeyCode::Char('J'),
                            modifiers: event::KeyModifiers::CONTROL,
                        })
                {
                    jumpcol_loop(io, view)?;

                    // Return from jumpcol mode
                    draw_default_footer(io, view)?;
                    io.io.flush()?;
                    continue;
                };

                // Redraw
                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('r'),
                        modifiers: event::KeyModifiers::NONE,
                    })
                    || kevent
                        == (KeyEvent {
                            code: KeyCode::Char('R'),
                            modifiers: event::KeyModifiers::NONE,
                        })
                {
                    draw_default_mode_screen(io, view)?;
                    continue;
                };

                // Shift to/from consensus view
                if kevent
                    == (KeyEvent {
                        code: KeyCode::Char('c'),
                        modifiers: event::KeyModifiers::NONE,
                    })
                    || kevent
                        == (KeyEvent {
                            code: KeyCode::Char('C'),
                            modifiers: event::KeyModifiers::NONE,
                        })
                {
                    if view.consensus().is_none() {
                        execute!(
                            io.io,
                            ResetColor,
                            terminal::Clear(ClearType::All),
                            cursor::MoveTo(0, 0),
                            Print("Calculating consensus..."),
                        )?;
                        view.calculate_consensus()
                    }
                    view.consensus = !view.consensus;

                    // Setting consensus moves everything a tick down, so we
                    // compensate by moving the view. Also, this prevents us
                    // from going out of bounds when toggling consensus
                    let delta = if view.consensus { 1 } else { -1 };
                    view.move_view(delta, 0);
                    draw_default_mode_screen(io, view)?;
                    continue;
                };
            }
            Event::Resize(ncols, nrows) => {
                view.resize(ncols, nrows);
                draw_default_mode_screen(io, view)?;
                continue;
            }
            _ => (),
        };
    }
    Ok(())
}

fn jumpcol_loop<T: Write>(io: &mut TerminalIO<T>, view: &mut View) -> Result<()> {
    let mut query = String::new();
    let mut invalid_column = false;
    loop {
        draw_jump_footer(io, &query, view, invalid_column)?;
        io.io.flush()?;
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
                                io,
                                view,
                                0,
                                n as isize - view.colstart as isize - 1,
                            )?;
                            break;
                        }
                    }
                    Err(_) => {
                        // Since we only accepts input keys 0-9, this should only
                        // happen if the user tries to input a number larger than
                        // the max usize.
                        invalid_column = true;
                    }
                }
            }
            _ => (),
        }
    }
    Ok(())
}

fn search_loop<T: Write>(io: &mut TerminalIO<T>, view: &mut View) -> Result<()> {
    let mut query = String::new();
    let mut last_result = None; // we change the prompt based on the last result
    loop {
        draw_search_footer(io, view, &query, last_result)?;
        io.io.flush()?;
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
                        move_view_and_redraw(io, view, dy, 0)?;
                        select_loop(io, view, nrow)?;
                        return Ok(());
                    }
                    // Else go to correct seq.
                    SearchResult::MatchSeq { row, start, stop } => {
                        let dy = row as isize
                            - (view.rowstart + (view.seq_nrows_display() / 2)) as isize;
                        let dx = start as isize
                            - (view.colstart + (view.seq_ncols_display() / 2)) as isize;
                        move_view_and_redraw(io, view, dy, dx)?;
                        let seq_row = view.display_row_of_index(row).unwrap();
                        let seq_col =
                            start.saturating_sub(view.colstart) as u16 + (view.namewidth + 1);
                        let highlight_str = unsafe {
                            std::str::from_utf8_unchecked(&view.seqs()[row][start..stop])
                        };
                        draw_highlight(io, seq_row, seq_col, view.term_ncols, highlight_str)?;
                        select_loop(io, view, row)?;
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
            return SearchResult::RegexErr;
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

// swaps rows around
fn select_loop<T: Write>(
    io: &mut TerminalIO<T>,
    view: &mut View,
    row_selected: usize,
) -> Result<()> {
    draw_select_footer(io, view)?;
    let mut row_selected = row_selected;
    loop {
        draw_names(io, view)?;
        highlight_name(io, view, row_selected)?;
        io.io.flush()?;
        if let Event::Key(kevent) = event::read()? {
            let delta = match kevent {
                // Quit to default view on escape or Ctrl-C or q
                KeyEvent {
                    code: KeyCode::Esc,
                    modifiers: event::KeyModifiers::NONE,
                }
                | KeyEvent {
                    code: KeyCode::Char('q'),
                    modifiers: event::KeyModifiers::NONE,
                }
                | KeyEvent {
                    code: KeyCode::Char('c'),
                    modifiers: event::KeyModifiers::CONTROL,
                } => break,

                // Move selection up/down
                KeyEvent {
                    code: KeyCode::Up,
                    modifiers: event::KeyModifiers::NONE,
                } => Some(-1),
                KeyEvent {
                    code: KeyCode::Down,
                    modifiers: event::KeyModifiers::NONE,
                } => Some(1),
                KeyEvent {
                    code: KeyCode::Up,
                    modifiers: event::KeyModifiers::SHIFT,
                } => Some(-10),
                KeyEvent {
                    code: KeyCode::Down,
                    modifiers: event::KeyModifiers::SHIFT,
                } => Some(10),
                KeyEvent {
                    code: KeyCode::Up,
                    modifiers: event::KeyModifiers::CONTROL,
                } => Some(isize::MIN),
                KeyEvent {
                    code: KeyCode::Down,
                    modifiers: event::KeyModifiers::CONTROL,
                } => Some(isize::MAX),

                _ => None,
            };
            if let Some(delta) = delta {
                row_selected = ((row_selected as isize).saturating_add(delta).max(0) as usize)
                    .min(view.nrows() - 1);
                shift_highlight_move(io, view, row_selected)?;
                continue;
            }
            let shift_delta = match kevent {
                KeyEvent {
                    code: KeyCode::Char('j'),
                    modifiers: event::KeyModifiers::NONE,
                } => 1,
                KeyEvent {
                    code: KeyCode::Char('k'),
                    modifiers: event::KeyModifiers::NONE,
                } => -1,
                KeyEvent {
                    code: KeyCode::Char('J'),
                    modifiers: event::KeyModifiers::SHIFT,
                } => 10,
                KeyEvent {
                    code: KeyCode::Char('K'),
                    modifiers: event::KeyModifiers::SHIFT,
                } => -10,
                KeyEvent {
                    code: KeyCode::Char('t'),
                    modifiers: event::KeyModifiers::NONE,
                } => isize::MIN,
                KeyEvent {
                    code: KeyCode::Char('b'),
                    modifiers: event::KeyModifiers::NONE,
                } => isize::MAX,
                _ => 0,
            };
            if shift_delta != 0 {
                let new_index = ((row_selected as isize).saturating_add(shift_delta).max(0)
                    as usize)
                    .min(view.nrows() - 1);
                view.move_row(row_selected, new_index);
                if shift_delta > -100 && shift_delta < 100 {
                    row_selected = new_index;
                    shift_highlight_move(io, view, row_selected)?;
                }
                draw_sequences(io, view)?;
            }
        }
    }
    Ok(())
}

fn shift_highlight_move<T: Write>(
    io: &mut TerminalIO<T>,
    view: &mut View,
    row: usize,
) -> Result<()> {
    if let Some(range) = view.seq_row_range() {
        if &row < range.start() {
            let dy = row as isize - *range.start() as isize;
            move_view_and_redraw(io, view, dy, 0)?;
        } else if &row > range.end() {
            let dy = row as isize - *range.end() as isize;
            move_view_and_redraw(io, view, dy, 0)?;
        };
    }
    Ok(())
}

fn draw_default_mode_screen<T: Write>(io: &mut TerminalIO<T>, view: &View) -> Result<()> {
    // We do need to reset the colors here, else I think the clearing of the
    // terminal will "clear" to the current color.
    if io.has_color {
        queue!(io.io, ResetColor)?;
    }
    execute!(io.io, terminal::Clear(ClearType::All),)?;

    if view.term_ncols < 2 || view.term_nrows < 4 {
        // This seems silly, but I have it because it allows me to assume a minimal
        // terminal size when drawing the regular alignment
        execute!(io.io, cursor::MoveTo(0, 0), Print(":("),)?;
    } else {
        draw_ruler(io, view)?;
        draw_default_footer(io, view)?;
        draw_names(io, view)?;
        draw_sequences(io, view)?;
    }
    Ok(io.io.flush()?)
}

#[derive(Parser)]
#[clap(version, author, about)]
struct AlenOptions {
    /// Path to alignment
    alignment: OsString,

    /// Display sequences in uppercase
    #[clap(short)]
    uppercase: bool,

    /// Force parsing as amino acids
    #[clap(short)]
    aminoacids: bool,

    /// Disable colors (monochrome, may improve lag)
    #[clap(short)]
    monochrome: bool,
}

fn main() {
    let args = AlenOptions::parse();

    // Check if file exists
    if args.alignment != "-" && !Path::new(&args.alignment).is_file() {
        println!("ERROR: Filename not found: {:?}", args.alignment);
        std::process::exit(1);
    }
    let file = std::fs::File::open(args.alignment).unwrap_or_else(|err| {
        println!("ERROR reading file: {}", err);
        std::process::exit(1)
    });

    let buffered_io = BufReader::new(file);
    let view = View::from_reader(BufReader::new(buffered_io), args.uppercase, args.aminoacids);
    match view {
        Err(e) => {
            println!("ERROR when loading FASTA: {}", e);
            std::process::exit(1);
        }
        Ok(mut view) => {
            let mut io = TerminalIO {
                io: stdout(),
                color: None,
                has_color: !args.monochrome,
            };
            terminal::enable_raw_mode().unwrap();
            execute!(io.io, terminal::EnterAlternateScreen, cursor::Hide,).unwrap();
            if let Err(e) = default_loop(&mut io, &mut view) {
                println!("Error: {}", e);
            }
            if io.has_color {
                queue!(io.io, ResetColor).unwrap();
            }
            execute!(io.io, cursor::Show, terminal::LeaveAlternateScreen).unwrap();
            terminal::disable_raw_mode().unwrap();
            std::process::exit(0);
        }
    }
}
