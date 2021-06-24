// To do: Load in sequences to alignment format
// To do: Figure out how to display to terminal
// To do: Auto-detect format and sequence kind
//

use bio::alphabets;
use bio::io::fasta;
use std::io::{self, stdin, stdout, BufRead, BufReader};
use std::time::Duration;

use crossterm::{
    event::{poll, read, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode},
    Result,
};

// TODO: Protein/DNA alignment?
#[derive(Debug)]
struct Alignment {
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
}

impl Alignment {
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
        Alignment{names, seqs}
    }
}

fn display(aln: Alignment) {
    loop {
        poll(Duration::from_secs(1_000_000_000)).unwrap();
        let event = read().unwrap();

        println!("Event::{:?}\r", event);

        // Break on Q or Esc
        if event == Event::Key(KeyCode::Esc.into()) || event == Event::Key(KeyCode::Char('q').into()) {
            break;
        }
    }
}

/*
fn main() {
    let buffered = BufReader::new(io::stdin());
    let aln = Alignment::new(buffered);
    println!("{:?}", aln);
}
 */

fn main() -> Result<()> {
    enable_raw_mode().expect("Error enabling raw mode in terminal");
    let aln = Alignment::new(BufReader::new(stdin()));
    display(aln);
    disable_raw_mode()
}