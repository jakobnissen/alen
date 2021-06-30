use crate::constants::{FOOTER_LINES, HEADER_LINES};

use std::cmp::{max, min};
use std::convert::TryInto;
use std::io::BufRead;
use std::ops::RangeInclusive;

use unicode_segmentation::UnicodeSegmentation;

use crossterm::terminal;

use bio::alphabets::Alphabet;
use bio::io::fasta;

/// A string that is separated into its constituent graphemes, used for printing
/// with uniform width.
/// The point of this struct is to avoid doing grapheme computations on ASCII
/// strings, and to only compute grapheme offsets once for each string.

pub struct Graphemes {
    pub string: String,
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
    pub fn len(&self) -> usize {
        match &self.grapheme_stop_indices {
            None => self.string.len(),
            Some(n) => n.len(),
        }
    }

    /// Get a string slice with the first N graphemes. If N is out of bounds,
    /// returns None.
    pub fn get_n_graphemes(&self, n: usize) -> Option<&str> {
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

pub struct Alignment {
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

/// A view object contains all information of what to draw to the screen
pub struct View {
    pub rowstart: usize, // zero-based index
    pub colstart: usize,
    pub term_nrows: u16, // obtained from terminal
    pub term_ncols: u16,
    pub namewidth: u16,
    aln: Alignment,
}

impl View {
    pub fn new(aln: Alignment) -> View {
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

    pub fn from_reader<T: BufRead>(file: T, uppercase: bool, must_aa: bool) -> View {
        View::new(Alignment::new(file, uppercase, must_aa))
    }

    pub fn graphemes(&self) -> &Vec<Graphemes> {
        &self.aln.graphemes
    }

    pub fn seqs(&self) -> &Vec<Vec<u8>> {
        &self.aln.seqs
    }

    pub fn is_aa(&self) -> bool {
        self.aln.is_aa
    }

    /// Resize the view to the current terminal window, but do not draw anything
    pub fn resize(&mut self, ncols: u16, nrows: u16) {
        // Get terminal size, and set it.
        let (oldcol, oldrow) = (self.colstart, self.rowstart);
        self.term_nrows = nrows;
        self.term_ncols = ncols;

        // Calculate new starts (if you zoom out)
        self.rowstart = calculate_start(oldrow, 0, self.seq_nrows_display(), self.nrows());
        self.colstart = calculate_start(oldcol, 0, self.seq_ncols_display(), self.ncols());

        // Calculate new namewidth
        if self.namewidth > self.term_ncols - 2 {
            let delta = (self.term_ncols as isize - 2) - self.namewidth as isize;
            self.resize_names(delta);
        }
    }

    pub fn move_view(&mut self, dy: isize, dx: isize) {
        self.rowstart = calculate_start(
            self.rowstart,
            dy,
            self.seq_nrows_display() as usize,
            self.nrows(),
        );
        self.colstart = calculate_start(
            self.colstart,
            dx,
            self.seq_ncols_display() as usize,
            self.ncols(),
        );
    }

    pub fn resize_names(&mut self, delta: isize) {
        let mut namewidth = (self.namewidth as isize) + delta;
        namewidth = max(0, namewidth); // not negative
        namewidth = min(namewidth, (self.term_ncols as isize).saturating_sub(2)); // do not exceed bounds

        // Do not exceed longest name shown on screen
        namewidth = min(namewidth, self.aln.longest_name as isize);
        self.namewidth = namewidth.try_into().unwrap();
    }

    pub fn nrows(&self) -> usize {
        self.aln.nrows()
    }

    pub fn ncols(&self) -> usize {
        self.aln.ncols()
    }

    pub fn seq_nrows_display(&self) -> usize {
        (self.term_nrows as usize).saturating_sub(HEADER_LINES + FOOTER_LINES)
    }

    pub fn seq_row_range(&self) -> Option<RangeInclusive<usize>> {
        match self.seq_nrows_display() {
            0 => None,
            nrows => Some(self.rowstart..=(min(self.nrows() - 1, self.rowstart + nrows - 1))),
        }
    }

    pub fn seq_ncols_display(&self) -> usize {
        self.term_ncols.saturating_sub(self.namewidth + 1).into() // one '|' char
    }

    pub fn seq_col_range(&self) -> Option<RangeInclusive<usize>> {
        match self.seq_ncols_display() {
            0 => None,
            ncols => Some(self.colstart..=(min(self.ncols() - 1, self.colstart + ncols - 1))),
        }
    }
}
