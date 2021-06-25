// To do: Auto-detect format and sequence kind
#![allow(dead_code)]
#![allow(unused_variables)]

use std::io::{stdin, BufRead, BufReader, Write};
use std::ops::RangeInclusive;
use bio::alphabets;
use bio::io::fasta;
use std::path::Path;
use std::time::Duration;
use std::cmp::{min, max};

use clap;
use unicode_segmentation::UnicodeSegmentation;

use crossterm::{
    cursor,
    event::{poll, read, Event, KeyCode},
    execute, queue,
    style::{self, Print, Stylize, Color, SetBackgroundColor, SetForegroundColor, ResetColor},
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

fn calculate_start(current: usize, displaysize: usize, n_rows_cols: usize) -> usize {
    let last_index = n_rows_cols - 1;
    // If we would exceed the alignment...
    if current + displaysize > last_index + 1 {
        // If we can display all of them, we are fixed at position 0
        if displaysize >= n_rows_cols {
            return 0
        } else {
            return n_rows_cols - displaysize
        }
    } else {
        return current
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
        view.resize();
        return view
    }

    /// Resize the view to the current terminal window, but do not draw anything
    fn resize(&mut self) {
        // Get terminal size, and set it.
        let (term_col, term_row) = terminal::size().unwrap();
        let (oldcol, oldrow) = (self.colstart, self.rowstart);
        self.term_nrows = term_row;
        self.term_ncols = term_col;
        
        // Calculate new starts (if you zoom out)
        self.rowstart = calculate_start(oldrow, self.nseqs_display(), self.aln.nrows());
        self.colstart = calculate_start(oldcol, self.ncols_display(), self.aln.ncols());

        // Set namewidth and padded names
        // TODO: Better calculation of namewidth, so it doesn't reset if you modify it
        self.namewidth = min(30, term_col >> 2);
        self.update_padded_names();
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
        for name in self.aln.names[self.seq_row_range()].iter() {
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

    fn nseqs_display(&self) -> usize {
        (self.term_nrows as usize) - (HEADER_LINES + FOOTER_LINES)
    }

    fn seq_row_range(&self) -> Option<RangeInclusive<usize>> {
        let nrows = self.nseqs_display();
        if nrows == 0 {
            return None
        } else
        {
            return Some((self.rowstart)..=(
                min(
                self.aln.nrows() - 1, // zero-based indexing
                self.rowstart + nrows - 1
                )
            ))
        }
    }

    fn ncols_display(&self) -> usize {
        (self.term_ncols - self.namewidth) as usize - 1 // one '|' char
    }

    fn seq_col_range(&self) -> Option<RangeInclusive<usize>> {
        let ncols = self.ncols_display();
        if ncols == 0 {
            return None
        } else {
            return Some((self.colstart)..=(
                min(
                self.term_ncols as usize,
                self.colstart + ncols - 1
                )
            ))
        }
    }
}

fn display(view: View) {
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

    io.flush().unwrap();
    loop {
        poll(Duration::from_secs(1_000_000_000)).unwrap();
        let event = read().unwrap();

        // Break on Q or Esc
        if event == Event::Key(KeyCode::Esc.into()) || event == Event::Key(KeyCode::Char('q').into()) {
            break;
        }
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
    if view.term_ncols < 3 || view.term_nrows < 4 {
        draw_easter_egg(io);
    } else {
        //draw_ruler(io, view);
        draw_names(io, view);
        draw_footer(io, view);
        draw_sequences(io, view);
    }
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
    for (i, (alnrow, padname)) in view.seq_row_range()
        .zip(view.padded_names.iter()).enumerate() {
        let termrow = (i + HEADER_LINES) as u16;
        queue!(io, cursor::MoveTo(0, termrow)).unwrap();
        print!("{}│", padname);
    }
}

fn draw_ruler<T: Write>(io: &mut T, view: &View) {
    unimplemented!()
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
    for (i, alnrow) in view.seq_row_range().enumerate() {
        // We have already checked this is valid ASCII in the Alignment
        // constructor.
        let seq = unsafe {
            std::str::from_utf8_unchecked(&view.aln.seqs[alnrow][view.seq_col_range()])
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
    let view = View::new(aln);
    display(view);
}
