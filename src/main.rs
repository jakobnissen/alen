// To do: Auto-detect format and sequence kind
// To do: Make stdin work on MacOS
// To do: Colors for AA
// To do: Better error handling in general
// To do: Factorize to files

use std::convert::TryInto;
use std::io::{stdin, BufRead, BufReader, Write};
use std::ops::RangeInclusive;
use bio::alphabets;
use bio::io::fasta;
use std::path::Path;
use std::cmp::{min, max};

use clap;
use unicode_segmentation::UnicodeSegmentation;

use crossterm::{
    cursor,
    event::{self, read, Event, KeyCode, KeyEvent},
    execute, queue,
    style::{self, Print, Color, SetBackgroundColor, SetForegroundColor, ResetColor},
    terminal::{self, disable_raw_mode, enable_raw_mode, ClearType},
};

const HEADER_LINES: usize = 2;
const FOOTER_LINES: usize = 1;

fn get_color_background(byte: u8) -> Option<Color> {
    match byte {
        b'a' | b'A' => Some(Color::AnsiValue(228)), // yellow
        b'c' | b'C' => Some(Color::AnsiValue(77)),  // green
        b'g' | b'G' => Some(Color::AnsiValue(39)),  // blue
        b't' | b'T' | b'u' | b'U' => Some(Color::AnsiValue(168)),  // pink
        _ => None
    }
}

// TODO: Protein/DNA alignment?
// TODO: Remove this debug
#[derive(Debug)]
struct Alignment {
    // TODO: Make graphemes vec vec str?
    graphemes: Vec<Vec<String>>,
    seqs: Vec<Vec<u8>>,
}

impl Alignment {
    fn nrows(&self) -> usize {
        self.graphemes.len()
    }

    fn ncols(&self) -> usize {
        self.seqs[0].len()
    }

    fn new<T: BufRead>(file: T, uppercase: bool) -> Alignment {
        let mut seqs = Vec::new();
        let mut graphemes: Vec<Vec<String>> = Vec::new();
        let reader = fasta::Reader::new(file);
        // Due to logic later, these MUST be ASCII
        let alphabet = alphabets::Alphabet::new(b"-ACMGRSVTUWYHKDBNacmgrsvtuwyhkdbn");
        let mut seqlength: Option<usize> = None;
        for result in reader.records() {
            // TODO: Better error message - file nume and record number, perhaps
            let record = result.expect("Error during FASTA parsing");
            if !alphabet.is_word(record.seq()) {
                panic!("Error: Sequence cannot be understood as DNA") // TODO: More precise error
            }
            let seq = record.seq().iter()
                .map(|&byte| {
                    // if uppercase, convert to char, uppercase, back to u8.
                    // we already check the alphabet, so we can be very lax about safety
                    if uppercase {
                        (byte as char).to_uppercase().next().unwrap() as u8
                    } else {
                        byte
                    }
                }).collect::<Vec<_>>();

            // Check identical sequence lengths
            if let Some(len) = seqlength {
                if seq.len() != len {
                    panic!("Error: Sequence lengths uneven") // TODO: More precise error
                }
            } else {
                seqlength = Some(seq.len())
            }
            seqs.push(seq);
            graphemes.push(UnicodeSegmentation::graphemes(record.id(), true).map(|s| s.to_owned()).collect());
        }
        // TODO: Warn if two headers are identical
        if seqlength.map_or(true, |i| i < 1) {
            panic!("Error: Empty alignment") // TODO: More precise error
        }
        Alignment{graphemes, seqs}
    }
}

/// A view object contains all information of what to draw to the screen
struct View {
    rowstart: usize,
    colstart: usize,
    term_nrows: u16, // obtained from terminal
    term_ncols: u16, // obtained from terminal
    namewidth: u16,
    aln: Alignment
}

fn calculate_start(
    current: usize,
    delta: isize,
    displaysize: usize,
    n_rows_cols: usize,
) -> usize {
    let last_index = n_rows_cols.saturating_sub(displaysize);
    let moveto = (current as isize).saturating_add(delta);
    if moveto < 0 {
        return 0
    } else if (moveto as usize) > last_index {
        return last_index
    } else {
        return moveto.try_into().unwrap()
    }
}

impl View {
    fn new(aln: Alignment) -> View {
        let mut view = View {
            rowstart: 0,
            colstart: 0,
            term_nrows: 0,
            term_ncols: 0,
            namewidth: 0,
            aln,
        };
        let (ncols, nrows) = terminal::size().unwrap();
        view.resize(ncols, nrows);
        return view
    }

    /// Resize the view to the current terminal window, but do not draw anything
    fn resize(&mut self, ncols: u16, nrows: u16) {
        // Get terminal size, and set it.
        let (oldcol, oldrow) = (self.colstart, self.rowstart);
        self.term_nrows = nrows;
        self.term_ncols = ncols;
        
        // Calculate new starts (if you zoom out)
        self.rowstart = calculate_start(
            oldrow, 0, 
            self.seq_nrows_display(), self.aln.nrows()
        );
        self.colstart = calculate_start(
            oldcol, 0, 
            self.seq_ncols_display(), self.aln.ncols()
        );

        // Set namewidth and padded names
        // TODO: Better calculation of namewidth, so it doesn't reset if you modify it
        self.resize_names(min(30, (ncols >> 2) as isize));
    }

    fn move_view<T: Write>(&mut self, io: &mut T, dy: isize, dx: isize) {
        self.rowstart = calculate_start(
            self.rowstart, dy, 
            self.term_nrows as usize, self.aln.nrows()
        );
        self.colstart = calculate_start(
            self.colstart, dx, 
            self.term_ncols as usize, self.aln.ncols()
        );
        if dy != 0 {
            draw_names(io, &self);
        }
        if dx != 0 {
            draw_ruler(io, &self)
        }
        draw_sequences(io, &self);
        io.flush().unwrap();
    }

    fn resize_names(&mut self, delta: isize) {
        let mut namewidth = (self.namewidth as isize) + delta;
        namewidth = max(0, namewidth); // not negative
        namewidth = min(namewidth, (self.term_ncols as isize).saturating_sub(2)); // do not exceed bounds
        // do not exceed longest name (TODO: CACHE THIS MAX?)
        namewidth = min(namewidth, self.aln.graphemes.iter().map(|g| g.len()).max().unwrap() as isize);
        self.namewidth = namewidth.try_into().unwrap();
    }

    fn seq_nrows_display(&self) -> usize {
        (self.term_nrows as usize).saturating_sub(HEADER_LINES + FOOTER_LINES)
    }

    fn seq_row_range(&self) -> Option<RangeInclusive<usize>> {
        match self.seq_nrows_display() {
            0 => None,
            nrows => Some(
                self.rowstart..=(min(self.aln.nrows() - 1, self.rowstart + nrows - 1))
            )
        }
    }

    fn seq_ncols_display(&self) -> usize {
        self.term_ncols.saturating_sub(self.namewidth + 1).into() // one '|' char
    }

    fn seq_col_range(&self) -> Option<RangeInclusive<usize>> {
        match self.seq_ncols_display() {
            0 => None,
            ncols => Some(
                self.colstart..=(min(self.aln.ncols() - 1, self.colstart + ncols - 1))
            )
        }
    }
}

fn display(view: &mut View) {
    let mut io = std::io::stdout();
    enable_raw_mode().unwrap();
    execute!(io,
        terminal::EnterAlternateScreen,
        cursor::Hide,
    ).unwrap();

    draw_all(&mut io, &view);

    io.flush().unwrap();
    loop {
        let event = read().unwrap();

        // Break on Q or Esc
        if event == Event::Key(KeyCode::Esc.into()) || event == Event::Key(KeyCode::Char('q').into()) {
            break;
        }

        match event {
            Event::Key(kevent) => {
                let delta = match kevent {
                    KeyEvent{code: KeyCode::Left, modifiers: event::KeyModifiers::NONE} => Some((0, -1)),
                    KeyEvent{code: KeyCode::Right, modifiers: event::KeyModifiers::NONE} => Some((0, 1)),
                    KeyEvent{code: KeyCode::Down, modifiers: event::KeyModifiers::NONE} => Some((1, 0)),
                    KeyEvent{code: KeyCode::Up, modifiers: event::KeyModifiers::NONE} => Some((-1, 0)),

                    // SHIFT: Move by 10
                    KeyEvent{code: KeyCode::Left, modifiers: event::KeyModifiers::SHIFT} => Some((0, -10)),
                    KeyEvent{code: KeyCode::Right, modifiers: event::KeyModifiers::SHIFT} => Some((0, 10)),
                    KeyEvent{code: KeyCode::Down, modifiers: event::KeyModifiers::SHIFT} => Some((10, 0)),
                    KeyEvent{code: KeyCode::Up, modifiers: event::KeyModifiers::SHIFT} => Some((-10, 0)),

                    // CONTROL: Move to end
                    KeyEvent{code: KeyCode::Left, modifiers: event::KeyModifiers::CONTROL} => Some((0, isize::MIN)),
                    KeyEvent{code: KeyCode::Right, modifiers: event::KeyModifiers::CONTROL} => Some((0, isize::MAX)),
                    KeyEvent{code: KeyCode::Down, modifiers: event::KeyModifiers::CONTROL} => Some((isize::MAX, 0)),
                    KeyEvent{code: KeyCode::Up, modifiers: event::KeyModifiers::CONTROL} => Some((isize::MIN, 0)),
                    _ => None
                };
                if let Some((dy, dx)) = delta {
                    view.move_view(&mut io, dy, dx)
                };
                let name_move = match kevent {
                    KeyEvent{code: KeyCode::Char(','), modifiers: event::KeyModifiers::NONE} => -1,
                    KeyEvent{code: KeyCode::Char('.'), modifiers: event::KeyModifiers::NONE} => 1,
                    _ => 0
                };
                if name_move != 0 {
                    view.resize_names(name_move);
                    draw_all(&mut io, &view);
                }

            },
            Event::Resize(ncols, nrows) => {
                view.resize(ncols, nrows);
                draw_all(&mut io, &view);
            }
            _ => ()
        };
    }

    execute!(
        io,
        style::ResetColor,
        cursor::Show,
        terminal::LeaveAlternateScreen
    ).unwrap();
    disable_raw_mode().unwrap();
}

fn draw_all<T: Write>(io: &mut T, view: &View) {
    execute!(
        io,
        style::ResetColor,
        terminal::Clear(ClearType::All),
    ).unwrap();
    if view.term_ncols < 3 || view.term_nrows < 4 {
        draw_easter_egg(io);
    } else {
        draw_ruler(io, view);
        draw_names(io, view);
        draw_footer(io, view);
        draw_sequences(io, view);
    }
    io.flush().unwrap();
}

// This seems silly, but I have it because it allows me to assume a minimal
// terminal size when drawing the regular alignment
fn draw_easter_egg<T: Write>(io: &mut T) {
    execute!(
        io,
        cursor::MoveTo(0, 0),
        style::Print(":("),
    ).unwrap();
}

fn draw_names<T: Write>(io: &mut T, view: &View) {
    if let Some(range) = view.seq_row_range() {
        for (i, nameindex) in range.into_iter().enumerate() {
            let termrow = (i + HEADER_LINES) as u16;
            let name = if view.namewidth == 0 {
                "".to_owned()
            } else {
                let graphemes = &view.aln.graphemes[nameindex];
                let elide = graphemes.len() > view.namewidth as usize;
                let namelen = min(graphemes.len(), view.namewidth as usize - elide as usize);
                let mut name = graphemes[0..namelen].join("");
                if elide {
                    name.push('…')
                }
                name
            };
            queue!(
                io,
                cursor::MoveTo(0, termrow),
                Print(name),
                cursor::MoveTo(view.namewidth, termrow),
                Print('│'),
            ).unwrap();
        }
    }
}

// TODO: This function is abhorrently complicated, make sure to rewrite it cleaner.
fn draw_ruler<T: Write>(io: &mut T, view: &View) {
    // Get tick positions
    let term_range = (view.namewidth + 1)..=(view.term_ncols-1);
    let aln_range = view.colstart..=(view.colstart+view.seq_ncols_display() - 1);

    // Check they must be same length (TODO: debug statement here?)
    let (aln_low, aln_high) = aln_range.clone().into_inner();
    assert!((aln_high - aln_low) + 1 == term_range.len());

    // In this loop we build the strings.
    let mut line_string = "┌".to_owned();
    let mut num_string = " ".to_owned();
    let mut beginning = true;
    for alncol in aln_range {
        // draw tick
        if (alncol + 1) % 10 == 0 {
            line_string.push('┴');
            let add = (alncol + 1).to_string();
            num_string.push_str(&add);
            num_string.push_str(&" ".repeat(10 - add.len()));
            beginning = false
        } else {
            line_string.push('─');
            if beginning {num_string.push(' ')};
        }
    }

    // Make sure it's not too long! The final byte is the leading space.
    num_string.truncate(term_range.len() + 1);

    queue!(
        io,
        cursor::MoveTo(view.namewidth, 0),
        Print(num_string),
        cursor::MoveTo(view.namewidth, 1),
        Print(line_string),
    ).unwrap();
}

fn draw_footer<T: Write>(io: &mut T, view: &View) {
    // First we create the full footer, then we truncate, if needed
    let mut footer = String::from(
        "q/Esc: Quit   ←/→/↑/↓ + None/Shift/Ctrl: Move alignment   ./,: Move names"
    );
    // Pad or truncate footer to match num columns
    let nchars = footer.chars().count();
    let ncols = view.term_ncols as usize;

    if nchars > ncols {
        footer = UnicodeSegmentation::graphemes(footer.as_str(), true)
            .take(ncols).collect::<String>();
    } else {
        footer.push_str(" ".repeat(ncols - nchars).as_str())
    }

    queue!(
        io,
        SetBackgroundColor(Color::Grey),
        SetForegroundColor(Color::Black),
        cursor::MoveTo(0, (view.term_nrows - 1) as u16),
        style::Print(footer),
        ResetColor,
    ).unwrap();
}

fn draw_sequences<T: Write>(io: &mut T, view: &View) {
    let row_range = match view.seq_row_range() {
        Some(n) => n,
        None => return
    };
    let col_range = match view.seq_col_range() {
        Some(n) => n,
        None => return
    };

    for (i, alnrow) in row_range.enumerate() {
        let termrow = (i + HEADER_LINES) as u16;
        queue!(
            io,
            cursor::MoveTo(view.namewidth + 1, termrow),
        ).unwrap();
        for col in col_range.clone() {
            let byte = view.aln.seqs[alnrow][col];
            let color = get_color_background(byte);
            match color {
                Some(clr) => queue!(
                    io,
                    SetForegroundColor(Color::Black),
                    SetBackgroundColor(clr),
                ).unwrap(),
                None => queue!(io, ResetColor).unwrap(),
            };
            queue!(io, Print(byte as char)).unwrap();
        }
    }

    queue!(
        io,
        ResetColor,
    ).unwrap();
}

fn main() {
    let args = clap::App::new("alen")
        .version("0.1")
        .author("Jakob Nybo Nissen <jakobnybonissen@gmail.com>")
        .about("Simple alignment viewer")
        .arg(clap::Arg::with_name("alignment")
            .help("Input alignment in FASTA format (- for stdin)")
            .takes_value(true)
            .required(true)
        ).arg(clap::Arg::with_name("uppercase")
            .short("u")
            .takes_value(false)
            .help("Displays sequences in uppercase")
        ).get_matches();

    let filename = args.value_of("alignment").unwrap();

    // Check if file exists
    if filename != "-" && !Path::new(filename).is_file() {
        println!("Error: Filename not found: \"{}\"", filename);
        std::process::exit(1);
    }

    let buffered_io: Box<dyn BufRead> = if filename == "-" {
        Box::new(BufReader::new(stdin()))
    } else {
        // TODO: Better error message?
        Box::new(BufReader::new(std::fs::File::open(filename).unwrap()))
    };
    let uppercase = args.is_present("uppercase");
    let aln = Alignment::new(BufReader::new(buffered_io), uppercase);
    let mut view = View::new(aln);
    display(&mut view);
}
