// To do: Auto-detect format and sequence kind
#![allow(dead_code)]
#![allow(unused_variables)]

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
    event::{poll, read, Event, KeyCode},
    execute, queue,
    style::{self, Print},
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

impl View {
    fn new(aln: Alignment) -> View {
        let (term_ncols, term_nrows) = terminal::size().unwrap();
        let namewidth = min(30, term_ncols >> 2);

        let mut view = View {
            rowstart: 0,
            colstart: 0,
            term_nrows,
            term_ncols,
            namewidth,
            padded_names: Vec::new(),
            aln
        };
        view.update_padded_names();
        return  view
    }

    /// Resize the view to the current terminal window, but do not draw anything
    fn resize(&mut self) {
        // TODO: Current minimum terminal size
        // 4 rows (2 header + 1 seq + 1 footer)
        // 3 cols (1 name + 1 | + 1 seq)

        // TODO:
        // * Check minimum size, print :( otherwise
        // * Re-calculate row/colstart (if you zoom out)
        // * update_padded_names
        unimplemented!()
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

    fn seq_row_range(&self) -> RangeInclusive<usize> {
        (self.rowstart)..=(
            min(
            self.aln.nrows() - 1, // zero-based indexing
            self.rowstart + self.nseqs_display() - 1
            )
        )
    }

    fn ncols_display(&self) -> usize {
        (self.term_ncols - self.namewidth) as usize - 1 // one '|' char
    }

    fn seq_col_range(&self) -> RangeInclusive<usize> {
        (self.colstart)..=(
            min(
            self.term_ncols as usize,
            self.colstart + self.ncols_display() - 1
            )
        )
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

    println!("Type [Esc] or 'q' to quit.");
    draw_names(&mut io, &view);

    loop {
        poll(Duration::from_secs(1_000_000_000)).unwrap();
        let event = read().unwrap();

        println!("Event::{:?}\r", event);

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
    draw_ruler(io, view);
    draw_names(io, view);
    draw_footer(io, view);
    draw_sequences(io, view);
}

fn draw_names<T: Write>(io: &mut T, view: &View) {    
    for (termrow, (alnrow, padname)) in view.seq_row_range()
        .zip(view.padded_names.iter()).enumerate() {
        
        queue!(io, cursor::MoveTo(0, (termrow + HEADER_LINES) as u16)).unwrap();
        print!("{}│", padname);
    }
    io.flush().unwrap();
}

fn draw_ruler<T: Write>(io: &mut T, view: &View) {
    unimplemented!()
}

fn draw_footer<T: Write>(io: &mut T, view: &View) {
    unimplemented!()
}

fn draw_sequences<T: Write>(io: &mut T, view: &View) {
    unimplemented!()
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
