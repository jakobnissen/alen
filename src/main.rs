// To do: Load in sequences to alignment format
// To do: Figure out how to display to terminal
// To do: Auto-detect format and sequence kind
//

use std::io::{stdin, BufRead, BufReader};
use bio::alphabets;
use bio::io::fasta;
use std::path::Path;
use clap::{Arg, App};

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

fn main() {
    let args = App::new("alen")
        .version("0.1")
        .author("Jakob Nybo Nissen <jakobnybonissen@gmail.com>")
        .about("Simple alignment viewer")
        .arg(Arg::with_name("alignment")
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

    let aln = Alignment::new(buffered_io);
    println!("{:?}", aln);
}
