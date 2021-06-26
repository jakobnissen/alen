// To do: Auto-detect format and sequence kind
// To do: Option to upper-case

use std::convert::TryInto;
use std::io::{stdin, BufRead, BufReader, Write};
use std::ops::RangeInclusive;
use bio::alphabets;
use bio::io::fasta;
use std::path::Path;
use std::time::Duration;
use std::cmp::min;

use clap;
use unicode_segmentation::UnicodeSegmentation;

use crossterm::{
    cursor,
    event::{self, poll, read, Event, KeyCode, KeyEvent},
    execute, queue,
    style::{self, Print, Color, SetBackgroundColor, SetForegroundColor, ResetColor},
    terminal::{self, disable_raw_mode, enable_raw_mode, ClearType},
};

const HEADER_LINES: usize = 2;
const FOOTER_LINES: usize = 1;

// TODO: Protein/DNA alignment?
// TODO: Remove this debug
#[derive(Debug)]
struct Alignment {
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
}

impl Alignment {
    fn nrows(&self) -> usize {
        self.names.len()
    }

    fn ncols(&self) -> usize {
        self.seqs[0].len()
    }

    fn new<T: BufRead>(file: T) -> Alignment {
        let mut seqs = Vec::new();
        let mut names: Vec<String> = Vec::new();
        let reader = fasta::Reader::new(file);
        // Due to logic later, these MUST be ASCII
        let alphabet = alphabets::Alphabet::new(b"-ACMGRSVTWYHKDBNacmgrsvtwyhkdbn");
        let mut seqlength: Option<usize> = None;
        for result in reader.records() {
            // TODO: Better error message - file nume and record number, perhaps
            let record = result.expect("Error during FASTA parsing");
            let seq: Vec<u8> = record.seq().into();
            if !alphabet.is_word(&seq) {
                panic!("Error: Sequence cannot be understood as DNA") // TODO: More precise error
            }

            // Check identical sequence lengths
            if let Some(len) = seqlength {
                if seq.len() != len {
                    panic!("Error: Sequence lengths uneven") // TODO: More precise error
                }
            } else {
                seqlength = Some(seq.len())
            }
            seqs.push(seq);
            names.push(record.id().to_owned());
        }
        // TODO: Warn if two headers are identical
        if seqlength.map_or(true, |i| i < 1) {
            panic!("Error: Empty alignment") // TODO: More precise error
        }
        Alignment{names, seqs}
    }
}

/// A view object contains all information of what to draw to the screen
struct View {
    rowstart: usize,
    colstart: usize,
    term_nrows: u16, // obtained from terminal
    term_ncols: u16, // obtained from terminal
    namewidth: u16,
    padded_names: Vec<String>,
    aln: Alignment
}

fn calculate_start(
    current: usize,
    delta: isize,
    displaysize: usize,
    n_rows_cols: usize,
) -> usize {
    let last_index = n_rows_cols - displaysize;
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
            padded_names: Vec::new(),
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
        self.namewidth = min(30, ncols >> 2);
        self.update_padded_names();
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
            self.update_padded_names();
            draw_names(io, &self);
        }
        if dx != 0 {
            draw_ruler(io, &self)
        }
        draw_sequences(io, &self);
        io.flush().unwrap();
    }

    /// Update the vector of padded_names, only including the ones that
    /// will be displayed on screen.
    /// Make sure to pad according to unicode graphemes, which I think should
    /// correspond to text width.
    fn update_padded_names(&mut self) {
        let width = self.namewidth as usize;
        if width == 0 {
            self.padded_names.fill("".to_owned());
            return ()
        } else {
            self.padded_names.clear();
        }

        // Else, we can guarantee width is at least 1
        if let Some(range) = self.seq_row_range() {
            for name in self.aln.names[range].iter() {
                let mut grapheme_iter = UnicodeSegmentation::graphemes(name.as_str(), true);
                let mut strvec = (&mut grapheme_iter).take(width - 1).collect::<Vec<_>>();
                if strvec.len() == (width - 1) {
                    strvec.push("…");
                } else if let Some(grapheme) = grapheme_iter.next() {
                    strvec.push(grapheme)
                }
                let mut str = strvec.join("");
                str.push_str(" ".repeat(width - strvec.len()).as_str());
                self.padded_names.push(str);
            }
        }
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
    execute!(io, terminal::EnterAlternateScreen).unwrap();
    
    queue!(
        io,
        style::ResetColor,
        terminal::Clear(ClearType::All),
        cursor::Hide,
        cursor::MoveTo(0, 0)
    ).unwrap();

    draw_all(&mut io, &view);

    //let mut file = std::fs::File::create("/tmp/foo.txt").unwrap();

    io.flush().unwrap();
    loop {
        poll(Duration::from_secs(1_000_000_000)).unwrap();
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
                }
            },
            Event::Resize(ncols, nrows) => {
                view.resize(ncols, nrows);
                draw_all(&mut io, &view);
            }
            _ => ()
        };

        /*
        file.write(
            format!(
        "----------------------\n{}\n{}\n{}\n{}\n{}\n{:?}\n{:?}\n",
        view.rowstart.to_string(), view.colstart.to_string(),
        view.term_ncols.to_string(), view.term_nrows.to_string(), view.namewidth.to_string(),
        view.seq_row_range(), view.seq_col_range()
        ).as_bytes()).unwrap();
        */
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
    if let Some(_range) = view.seq_row_range() {
        for (i, padname) in view.padded_names.iter().enumerate() {
            let termrow = (i + HEADER_LINES) as u16;
            queue!(io, cursor::MoveTo(0, termrow)).unwrap();
            print!("{}│", padname);
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
    if (aln_high - aln_low) + 1 != term_range.len() {
        panic!("{} {} {} {} {} {}", (aln_high - aln_low) + 1, term_range.len(),
        view.colstart, view.seq_ncols_display(), view.term_ncols, view.namewidth);
    }

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
        "q/Esc: Quit   ←→↑↓: Move alignment   Ctrl+←: Move to end"
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
        // We have already checked this is valid ASCII in the Alignment
        // constructor.
        let seq = unsafe {
            std::str::from_utf8_unchecked(&view.aln.seqs[alnrow][col_range.clone()])
        };
        let termrow = (i + HEADER_LINES) as u16;
        queue!(
            io,
            cursor::MoveTo(view.namewidth + 1, termrow),
            style::Print(seq),
        ).unwrap();
    }
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
    let aln = Alignment::new(BufReader::new(buffered_io));
    let mut view = View::new(aln);
    display(&mut view);
}
