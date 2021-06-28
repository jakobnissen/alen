// To do: Better error handling in general
// To do: Factorize to files
// To do: More formats?
// To do: Jump to column?
// To do: Improve performance:
// First improve drawing sequence: Queue only the minimal amount of ops, only
// change color when there is actually a color change.
// Then implement multithreading. One thread reads input and moves view etc,
// another draws. Drawings can be "skipped" if there are still queued inputs, perhaps?
// To do: Possibly show deviation from consensus?
// To do: Add search - verify input to match current alphabet when typing, and then just do findfirst? Or Regex with i?

use std::cmp::{max, min};
use std::convert::TryInto;
use std::io::{BufRead, BufReader, Write};
use std::ops::RangeInclusive;
use std::path::Path;

use bio::alphabets::Alphabet;
use bio::io::fasta;

use unicode_segmentation::UnicodeSegmentation;

use crossterm::{
    cursor,
    event::{self, Event, KeyCode, KeyEvent},
    execute, queue,
    style::{Color, Print, ResetColor, SetBackgroundColor, SetForegroundColor},
    terminal::{self, ClearType},
};

const HEADER_LINES: usize = 2;
const FOOTER_LINES: usize = 1;

/// Panics if not valid biosequence, else returns true (aa) or false (dna)
fn verify_alphabet(seqs: &[Vec<u8>], graphemes_vec: &[Vec<String>], must_aa: bool) -> bool {
    let dna_alphabet = Alphabet::new(b"-ACMGRSVTUWYHKDBNacmgrsvtuwyhkdbn");
    let aa_alphabet = Alphabet::new(b"*-ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

    let mut valid_dna = true;
    for (seq, graphemes) in seqs.iter().zip(graphemes_vec) {
        if !must_aa {
            valid_dna &= dna_alphabet.is_word(seq);
        }
        // DNA alphabet is a subset of AA alphabet, so we panic if it can't even be AA
        if !aa_alphabet.is_word(seq) {
            println!(
                "ERROR:Sequence \"{}\" cannot be understood as amino acids.",
                graphemes.join("")
            );
            std::process::exit(1);
        }
    }
    must_aa | !valid_dna
}

fn make_uppercase(seqs: &mut [Vec<u8>]) {
    // We exploit the fact that only [A-Za-z\-\*] is allowed. Uppercasing this
    // means setting the third-to-top bit to 0. For - or *, we don't change the bit.
    for seq in seqs.iter_mut() {
        for byte in seq.iter_mut() {
            *byte &= !(((*byte >= b'A') as u8) << 5)
        }
    } 

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

struct Alignment {
    // TODO: Make graphemes vec vec str?
    graphemes: Vec<Vec<String>>,
    longest_name: usize,
    seqs: Vec<Vec<u8>>,
    is_aa: bool,
}

impl Alignment {
    fn nrows(&self) -> usize {
        self.seqs.len()
    }

    fn ncols(&self) -> usize {
        self.seqs[0].len()
    }

    fn new<T: BufRead>(file: T, uppercase: bool, must_aa: bool) -> Alignment {
        let mut seqs = Vec::new();
        let mut graphemes: Vec<Vec<String>> = Vec::new();
        let reader = fasta::Reader::new(file);
        let mut seqlength: Option<usize> = None;
        for result in reader.records() {
            let record = match result {
                Ok(r) => r,
                Err(e) => {
                    println!("ERROR: During FASTA parsing, found error:");
                    println!("{:?}", e);
                    std::process::exit(1);
                }
            };
            let seq = record.seq().iter().copied().collect::<Vec<_>>();

            // Check identical sequence lengths
            if let Some(len) = seqlength {
                if seq.len() != len {
                    println!(
                        "ERROR: Not all input sequences are the same length. \
                    Expected sequence length {}, found {}.",
                        len,
                        seq.len()
                    );
                    std::process::exit(1);
                }
            } else {
                seqlength = Some(seq.len())
            }
            seqs.push(seq);
            graphemes.push(
                UnicodeSegmentation::graphemes(record.id(), true)
                    .map(|s| s.to_owned())
                    .collect(),
            );
        }

        // Verify alphabet
        let is_aa = verify_alphabet(&seqs, &graphemes, must_aa);

        // Turn uppercase if requested
        if uppercase {
            make_uppercase(&mut seqs);
        }

        if seqlength.map_or(true, |i| i < 1) {
            println!("ERROR: Alignment has zero sequences, or has length 0.");
            std::process::exit(1);
        }

        let longest_name = graphemes.iter().map(|v| v.len()).max().unwrap();
        Alignment {
            graphemes,
            longest_name,
            seqs,
            is_aa,
        }
    }
}

/// A view object contains all information of what to draw to the screen
struct View {
    rowstart: usize, // zero-based index
    colstart: usize,
    term_nrows: u16, // obtained from terminal
    term_ncols: u16, // obtained from terminal
    namewidth: u16,
    aln: Alignment,
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
        view
    }

    /// Resize the view to the current terminal window, but do not draw anything
    fn resize(&mut self, ncols: u16, nrows: u16) {
        // Get terminal size, and set it.
        let (oldcol, oldrow) = (self.colstart, self.rowstart);
        self.term_nrows = nrows;
        self.term_ncols = ncols;

        // Calculate new starts (if you zoom out)
        self.rowstart = calculate_start(oldrow, 0, self.seq_nrows_display(), self.aln.nrows());
        self.colstart = calculate_start(oldcol, 0, self.seq_ncols_display(), self.aln.ncols());

        // Set namewidth and padded names
        // TODO: Better calculation of namewidth, so it doesn't reset if you modify it
        self.resize_names(min(30, (ncols >> 2) as isize));
    }

    fn move_view<T: Write>(&mut self, io: &mut T, dy: isize, dx: isize) {
        let old_rowstart = self.rowstart;
        let old_colstart = self.colstart;

        self.rowstart = calculate_start(
            self.rowstart,
            dy,
            self.seq_nrows_display() as usize,
            self.aln.nrows(),
        );
        self.colstart = calculate_start(
            self.colstart,
            dx,
            self.seq_ncols_display() as usize,
            self.aln.ncols(),
        );

        // Only update the view if the view was actually moved.
        if self.rowstart != old_rowstart {
            draw_names(io, self);
        }
        if self.colstart != old_colstart {
            draw_ruler(io, self)
        }
        if self.rowstart != old_rowstart || self.colstart != old_colstart {
            draw_sequences(io, self);
            io.flush().unwrap();
        }
    }

    fn resize_names(&mut self, delta: isize) {
        let mut namewidth = (self.namewidth as isize) + delta;
        namewidth = max(0, namewidth); // not negative
        namewidth = min(namewidth, (self.term_ncols as isize).saturating_sub(2)); // do not exceed bounds
        // do not exceed longest name shown on screen
        namewidth = min(namewidth, self.aln.longest_name as isize);
        self.namewidth = namewidth.try_into().unwrap();
    }

    fn seq_nrows_display(&self) -> usize {
        (self.term_nrows as usize).saturating_sub(HEADER_LINES + FOOTER_LINES)
    }

    fn seq_row_range(&self) -> Option<RangeInclusive<usize>> {
        match self.seq_nrows_display() {
            0 => None,
            nrows => Some(self.rowstart..=(min(self.aln.nrows() - 1, self.rowstart + nrows - 1))),
        }
    }

    fn seq_ncols_display(&self) -> usize {
        self.term_ncols.saturating_sub(self.namewidth + 1).into() // one '|' char
    }

    fn seq_col_range(&self) -> Option<RangeInclusive<usize>> {
        match self.seq_ncols_display() {
            0 => None,
            ncols => Some(self.colstart..=(min(self.aln.ncols() - 1, self.colstart + ncols - 1))),
        }
    }
}

fn calculate_start(current: usize, delta: isize, displaysize: usize, n_rows_cols: usize) -> usize {
    let last_index = n_rows_cols.saturating_sub(displaysize);
    let moveto = (current as isize).saturating_add(delta);
    if moveto < 0 {
        0
    } else if (moveto as usize) > last_index {
        last_index
    } else {
        moveto.try_into().unwrap()
    }
}

fn display(view: &mut View) {
    let mut io = std::io::stdout();
    terminal::enable_raw_mode().unwrap();
    execute!(io, terminal::EnterAlternateScreen, cursor::Hide,).unwrap();

    draw_all(&mut io, view);
    io.flush().unwrap();

    loop {
        let event = event::read().unwrap(); // TODO: This will error taking stdin

        // Break on Q or Esc
        if event == Event::Key(KeyCode::Esc.into())
            || event == Event::Key(KeyCode::Char('q').into())
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
                    view.move_view(&mut io, dy, dx)
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
                        draw_all(&mut io, view);
                    }
                }
            }
            Event::Resize(ncols, nrows) => {
                view.resize(ncols, nrows);
                draw_all(&mut io, view);
            }
            _ => (),
        };
    }

    execute!(io, ResetColor, cursor::Show, terminal::LeaveAlternateScreen).unwrap();
    terminal::disable_raw_mode().unwrap();
}

fn draw_all<T: Write>(io: &mut T, view: &View) {
    execute!(io, ResetColor, terminal::Clear(ClearType::All),).unwrap();

    if view.term_ncols < 2 || view.term_nrows < 4 {
        // This seems silly, but I have it because it allows me to assume a minimal
        // terminal size when drawing the regular alignment
        execute!(io, cursor::MoveTo(0, 0), Print(":("),).unwrap();
    } else {
        draw_ruler(io, view);
        draw_names(io, view);
        draw_footer(io, view);
        draw_sequences(io, view);
    }
    io.flush().unwrap();
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
            )
            .unwrap();
        }
    }
}

// TODO: This function is abhorrently complicated, make sure to rewrite it cleaner.
fn draw_ruler<T: Write>(io: &mut T, view: &View) {
    // Get tick positions
    let term_range = (view.namewidth + 1)..=(view.term_ncols - 1);
    let aln_range = view.colstart..=(view.colstart + view.seq_ncols_display() - 1);

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
            if beginning {
                num_string.push(' ')
            };
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
    )
    .unwrap();
}

fn draw_footer<T: Write>(io: &mut T, view: &View) {
    // First we create the full footer, then we truncate, if needed
    let mut footer =
        String::from("q/Esc: Quit   ←/→/↑/↓ + None/Shift/Ctrl: Move alignment   ./,: Adjust names");
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

    queue!(
        io,
        SetBackgroundColor(Color::Grey),
        SetForegroundColor(Color::Black),
        cursor::MoveTo(0, (view.term_nrows - 1) as u16),
        Print(footer),
        ResetColor,
    )
    .unwrap();
}

// This function is the performance bottleneck, so I try to queue as little as
// possible to the terminal by checking every single operation to see if it's
// necessary
fn draw_sequences<T: Write>(io: &mut T, view: &View) {
    let row_range = match view.seq_row_range() {
        Some(n) => n,
        None => return,
    };
    let col_range = match view.seq_col_range() {
        Some(n) => n,
        None => return,
    };

    let mut oldcolor: Option<Color> = None;
    let mut is_foreground_black = false;
    for (i, alnrow) in row_range.enumerate() {
        let termrow = (i + HEADER_LINES) as u16;
        queue!(io, cursor::MoveTo(view.namewidth + 1, termrow),).unwrap();
        for byte in view.aln.seqs[alnrow][col_range.clone()].iter() {
            let color = if view.aln.is_aa {
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
                        queue!(io, SetForegroundColor(Color::Black)).unwrap();
                        is_foreground_black = true;
                    };
                    queue!(io, SetBackgroundColor(clr)).unwrap();
                // Otherwise, reset to default color scheme
                } else {
                    queue!(io, ResetColor).unwrap();
                    is_foreground_black = false;
                }
            };
            queue!(io, Print(*byte as char)).unwrap();
            oldcolor = color;
        }
    }

    queue!(io, ResetColor,).unwrap();
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

    let buffered_io = BufReader::new(std::fs::File::open(filename).unwrap());
    let uppercase = args.is_present("uppercase");
    let must_aa = args.is_present("aminoacids");
    let aln = Alignment::new(BufReader::new(buffered_io), uppercase, must_aa);
    let mut view = View::new(aln);
    display(&mut view);
}
