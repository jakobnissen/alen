// To do: Load in sequences to alignment format
// To do: Figure out how to display to terminal
// To do: Auto-detect format and sequence kind
//

use std::io::{stdin, BufRead, BufReader, Write};
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

struct View {
    rowstart: usize,
    colstart: usize,
    nrows: u16,
    ncols: u16,
    namewidth: u16,
    padded_names: Vec<String>,
    aln: Alignment
}

impl View {
    fn new(aln: Alignment) -> View {
        let (nrows, ncols) = terminal::size().unwrap();
        let namewidth = min(30, nrows >> 2);
        let padded_names = View::padded_names(&aln, namewidth);
        View {rowstart: 1, colstart: 1, nrows, ncols, namewidth, padded_names, aln}
    }

    fn padded_names(aln: &Alignment, width: u16) -> Vec<String> {
        let mut result = Vec::new();
        for name in aln.names.iter() {
            let pvec = UnicodeSegmentation::graphemes(name.as_str(), true)
                .take(width as usize).collect::<Vec<_>>();
            let mut str = pvec.join("");
            str.push_str(" ".repeat((width as usize) - pvec.len()).as_str());
            result.push(str);
        }
        return result
    }

    fn rowmax(&self) -> usize {
        min(self.aln.nrows() - 1, self.rowstart + self.nrows as usize - 1)
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
    draw_names(&io, &view);

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

fn draw_names<T: Write>(mut io: T, view: &View) {
    //let first_tick_col = view.namewidth as u64 + 1 + (view.colstart as u64 % 10);
    queue!(io, cursor::MoveTo(0, 0)).unwrap();
    let mut ruler = String::new(); // ┌─┴
    let mut counter = String::new();
    queue!(
            io,
            cursor::MoveTo(0, 0),
            Print(counter),
            cursor::MoveTo(0, 1),
            Print(ruler),
    ).unwrap();
    
    let ruler_height = 2;
    for (termrow, alnrow) in (view.rowstart..view.rowmax()).enumerate() {
        queue!(io, cursor::MoveTo(0, (termrow + ruler_height) as u16)).unwrap();
        print!("{}|", view.padded_names[alnrow]);
    }
    io.flush().unwrap();
}

fn draw_bottom<T: Write>(mut io: T, view: View) {
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
