// To do: Better error handling in general
// To do: Factorize to files
// To do: More formats?
// To do: Make stdin work on MacOS
// To do: Jump to column?
// To do: Implement multithreading. One thread reads input and moves view etc,
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
fn verify_alphabet(seqs: &[Vec<u8>], graphemes_vec: &[Graphemes], must_aa: bool) -> bool {
    let dna_alphabet = Alphabet::new(b"-ACMGRSVTUWYHKDBNacmgrsvtuwyhkdbn");
    let aa_alphabet = Alphabet::new(b"*-ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

    let mut valid_dna = true;
    for (seq, graphemes) in seqs.iter().zip(graphemes_vec) {
        if !must_aa && valid_dna {
            valid_dna &= dna_alphabet.is_word(seq);
        }
        // DNA alphabet is a subset of AA alphabet, so we panic if it can't even be AA
        if !aa_alphabet.is_word(seq) {
            println!(
                "ERROR:Sequence \"{}\" cannot be understood as amino acids.",
                graphemes.string
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

/// A string that is separated into its constituent graphemes, used for printing
/// with uniform width.
/// The point of this struct is to avoid doing grapheme computations on ASCII
/// strings, and to only compute grapheme offsets once for each string.

struct Graphemes {
    string: String,
    /// If the string is ASCII (it usually is), we don't bother saving this.
    grapheme_stop_indices: Option<Vec<usize>>,
}

impl Graphemes {
    fn new(st: &str) -> Graphemes {
        let string = st.to_owned();
        // If string is ASCII, we save only the string itself, and don't bother
        // to do grapheme identification, since 1 grapheme == 1 byte == 1 char
        let grapheme_stop_indices = if string.is_ascii() {
            None
        } else {
            // Else we add in the LAST byte of each graphemes in vector V,
            // such that the first N graphemes of the string are encoded by
            // the bytes 0..=V[N-1].
            Some({
                let mut v: Vec<usize> =
                    UnicodeSegmentation::grapheme_indices(string.as_str(), true)
                        .skip(1)
                        // The iterator gives start indices, I assume end indices of the
                        // previous grapheme is the previous byte
                        .map(|(index, _grapheme)| index - 1)
                        .collect();
                // End byte of last grapheme is just the last byte-index of the string
                v.push(string.len());
                v
            })
        };
        Graphemes {
            string,
            grapheme_stop_indices,
        }
    }

    /// Number of graphemes in string.
    fn len(&self) -> usize {
        match &self.grapheme_stop_indices {
            None => self.string.len(),
            Some(n) => n.len(),
        }
    }

    /// Get a string slice with the first N graphemes. If N is out of bounds,
    /// returns None.
    fn get_n_graphemes(&self, n: usize) -> Option<&str> {
        if n > self.len() {
            None
        } else {
            match &self.grapheme_stop_indices {
                None => Some(&self.string[0..n]),
                Some(v) => {
                    if n == 0 {
                        Some("")
                    } else {
                        Some(&self.string[0..=v[n - 1]])
                    }
                }
            }
        }
    }
}

struct Alignment {
    graphemes: Vec<Graphemes>,
    // longest as in number of graphemes. We cache this for efficiency, it can be
    // computed from the graphemes field easily
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
        let mut graphemes: Vec<Graphemes> = Vec::new();
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
            graphemes.push(Graphemes::new(record.id()));
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
    term_ncols: u16,
    namewidth: u16,
    aln: Alignment,
}

impl View {
    fn new(aln: Alignment) -> View {
        let (ncols, nrows) = terminal::size().unwrap();
        let mut view = View {
            rowstart: 0,
            colstart: 0,
            term_nrows: 0,
            term_ncols: 0,
            namewidth: 0,
            aln,
        };
        // We need to resize before we resize names, because the latter
        // depends on a nonzero terminal size.
        view.resize(ncols, nrows);
        view.resize_names((ncols >> 2) as isize);
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

        // Calculate new namewidth
        if self.namewidth > self.term_ncols - 2 {
            let delta = (self.term_ncols as isize - 2) - self.namewidth as isize;
            self.resize_names(delta);
        }
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

        // Do not exceed longest name shown on screen
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
        let event = event::read().unwrap();

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
            )
            .unwrap();
        }
    }
}

fn draw_ruler<T: Write>(io: &mut T, view: &View) {
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

    queue!(
        io,
        cursor::MoveTo(view.namewidth, 0),
        Print(num_string),
        cursor::MoveTo(view.namewidth, 1),
        Print(tick_string),
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
