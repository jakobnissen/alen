use crate::constants::{FOOTER_LINES, HEADER_LINES};

use std::cmp::{max, min};
use std::convert::TryInto;
use std::io::BufRead;
use std::ops::RangeInclusive;

use anyhow::{anyhow, Result};

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
    pub fn new(st: &str) -> Graphemes {
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
                let v: Vec<usize> = UnicodeSegmentation::grapheme_indices(string.as_str(), true)
                    // First is always 0
                    .skip(1)
                    .map(|(index, _grapheme)| index)
                    .collect();
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
            None => self.string.len(), // if is ASCII
            Some(n) => n.len() + 1,    // we skipped first one
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
                        let lastindex = if n == self.len() {
                            self.string.len()
                        } else {
                            v[n]
                        };
                        Some(&self.string[0..lastindex])
                    }
                }
            }
        }
    }
}

// Returns whether it's an AA alphabet, else it default to DNA.
fn verify_alphabet(seqs: &[Vec<u8>], graphemes_vec: &[Graphemes], must_aa: bool) -> Result<bool> {
    let dna_alphabet = Alphabet::new(b"-ACMGRSVTUWYHKDBNacmgrsvtuwyhkdbn");
    let aa_alphabet = Alphabet::new(b"*-ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

    let mut valid_dna = !must_aa;
    for (seq, graphemes) in seqs.iter().zip(graphemes_vec) {
        if valid_dna {
            if dna_alphabet.is_word(seq) {
                continue;
            } else {
                valid_dna = false;
            }
        }
        if !aa_alphabet.is_word(seq) {
            return Err(anyhow!(
                "Sequence \"{}\" cannot be understood as amino acids.",
                graphemes.string
            ));
        }
    }
    Ok(!valid_dna)
}

fn calculate_consensus(seqs: &[Vec<u8>], is_aa: bool) -> Vec<Option<u8>> {
    // We have verified seq is *, - or A-Z, a-z. We uppercase, by masking 3rd bit.
    // * and - are 27th and 28th elements, giving us 28 possible elements total
    let offset = b'A';
    let ncols = seqs[0].len();
    let mut counts = vec![[0u32; 28]; ncols];

    // First loop over sequneces in memory order
    for seq in seqs.iter() {
        for (byte, arr) in seq.iter().zip(counts.iter_mut()) {
            let index = match byte {
                b'*' => 26,
                b'-' => 27,
                // Unset third bit to uppercase ASCII letters
                b => ((b & 0b11011111) - offset) as usize,
            };
            arr[index] += 1
        }
    }

    // We set all ambiguous bases/AAs to 0.
    // These are useless in the consensus
    for arr in counts.iter_mut() {
        for &i in if is_aa {
            b"BJOUXZZZZZZ"
        } else {
            b"MRSVWYHKDBN"
        } {
            arr[(i - offset) as usize] = 0;
        }
    }

    return counts
        .iter()
        .map(|arr| {
            let (mut most_common_byte, count) = arr
                .iter()
                .enumerate()
                .max_by_key(|(_, &x)| x)
                .map(|(i, cnt)| (i as u8 + offset, cnt))
                .unwrap();

            // If the most common is * or -, it becomes \n and \r after unsetting the bit.
            most_common_byte = match most_common_byte {
                b'[' => b'*',
                b'\\' => b'-',
                _ => most_common_byte,
            };

            if *count == 0 {
                None
            } else {
                Some(most_common_byte)
            }
        })
        .collect();
}

fn move_element<T>(v: &mut Vec<T>, from: usize, to: usize) -> Option<()> {
    if from.max(to) >= v.len() {
        return None;
    }
    if from == to {
        return Some(());
    }
    let mut i = from;
    let delta: isize = if to > from { 1 } else { -1 };
    while i != to {
        let i2 = (i as isize + delta) as usize;
        v.swap(i, i2);
        i = i2;
    }
    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_move_element() {
        let mut v = vec![1, 2, 3, 4, 5];
        move_element(&mut v, 1, 3);
        assert_eq!(v, vec![1, 3, 4, 2, 5]);
        move_element(&mut v, 3, 1);
        assert_eq!(v, vec![1, 2, 3, 4, 5]);
        assert_eq!(move_element(&mut v, 0, 5), None);
        assert_eq!(v, vec![1, 2, 3, 4, 5]);
        move_element(&mut v, 3, 0);
        assert_eq!(v, vec![4, 1, 2, 3, 5]);
    }
}

pub struct Alignment {
    graphemes: Vec<Graphemes>,
    // longest as in number of graphemes. We cache this for efficiency, it can be
    // computed from the graphemes field easily
    longest_name: usize,
    seqs: Vec<Vec<u8>>,
    // we calculate this lazily upon demand
    consensus: Option<Vec<Option<u8>>>,
    is_aa: bool,
    ordered: bool,
}

impl Alignment {
    fn nrows(&self) -> usize {
        self.seqs.len()
    }

    fn ncols(&self) -> usize {
        self.seqs[0].len()
    }

    fn new<T: BufRead>(file: T, uppercase: bool, must_aa: bool) -> Result<Alignment> {
        let mut seqs = Vec::new();
        let mut graphemes: Vec<Graphemes> = Vec::new();
        let reader = fasta::Reader::new(file);
        let mut seqlength: Option<usize> = None;
        for result in reader.records() {
            let record = result?;
            let header = Graphemes::new(record.id());
            let seq = record.seq().to_vec();

            // Check identical sequence lengths
            if let Some(len) = seqlength {
                if seq.len() != len {
                    return Err(anyhow!(
                        "Not all input sequences are the same length. \
                    Expected sequence length {}, seq \"{}\" has length {}.",
                        len,
                        &header.string,
                        seq.len()
                    ));
                }
            } else {
                seqlength = Some(seq.len())
            }
            graphemes.push(header);
            seqs.push(seq);
        }

        // Turn uppercase if requested
        if uppercase {
            seqs.iter_mut().for_each(|s| s.make_ascii_uppercase());
        }

        // Verify alphabet
        let is_aa = verify_alphabet(&seqs, &graphemes, must_aa)?;

        if seqlength.map_or(true, |i| i < 1) {
            return Err(anyhow!("Alignment has no seqs, or seqs have length 0."));
        }

        let longest_name = graphemes.iter().map(|v| v.len()).max().unwrap();
        Ok(Alignment {
            graphemes,
            longest_name,
            seqs,
            consensus: None,
            is_aa,
            ordered: false
        })
    }

    /// Reorder the vectors of the alignment such that similar rows are next to each other.
    fn order(&mut self) {
        // If already ordered or 2 or fewer rows, ordering doesn't matter
        if self.ordered || self.nrows() < 3 {
            return;
        }

        // Choices only appear when placing the 3rd seq, so first two are given.
        let mut order: Vec<usize> = vec![0, 1];
        let mut neighbor_distances = vec![jaccard_distance(&self.seqs[0], &self.seqs[1])];
        let mut distances: Vec<usize> = Vec::with_capacity(self.seqs.len());

        for seqindex in 2..self.nrows() {
            // Jaccard distances from seq[seqindex] to all that are already placed.
            distances.clear();
            distances.extend(
                order
                    .iter()
                    .map(|i| jaccard_distance(&self.seqs[seqindex], &self.seqs[*i])),
            );

            // Find minimum distance. First check the "ends", that is, the distance
            // if the new sequence is placed at top or bottom.
            let (mut min_dist, mut min_index) =
                if distances.first().unwrap() < distances.last().unwrap() {
                    (*distances.first().unwrap(), 0)
                } else {
                    (*distances.last().unwrap(), seqindex)
                };

            // Now check the distances if the new sequence is placed in between two existing.
            // This adds 2 new neighbor distances, but removes one.
            for i in 1..order.len() {
                let new_dist = distances[i - 1] + distances[i] - neighbor_distances[i - 1];
                if new_dist < min_dist {
                    min_dist = new_dist;
                    min_index = i;
                }
            }

            order.insert(min_index, seqindex);
            // If min_index is 0, the new sequence is placed at beginning
            if min_index == 0 {
                neighbor_distances.insert(0, min_dist);
            }
            // If it's the last, then we just place it at the end
            else if min_index == seqindex {
                neighbor_distances.push(min_dist);
            }
            // Else, we place it in the middle of two sequences. This removes (overwrites) the
            // neighbor distance of the two seqs that are now separated, but adds two new
            // neighbor pairs - for both the new seq's new neighbors.
            else {
                neighbor_distances[min_index - 1] = distances[min_index - 1];
                neighbor_distances.insert(min_index, distances[min_index]);
            }
        }

        // reorder internal vectors
        if let Some(c) = self.consensus.as_deref_mut() {
            sort_by_indices(c, order.clone())
        };
        sort_by_indices(&mut self.seqs, order.clone());
        sort_by_indices(&mut self.graphemes, order);

        self.ordered = true
    }
}

fn jaccard_distance(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).map(|(i, j)| (i != j) as usize).sum()
}

fn sort_by_indices<T>(data: &mut [T], mut indices: Vec<usize>) {
    for idx in 0..data.len() {
        if indices[idx] != idx {
            let mut current_idx = idx;
            loop {
                let target_idx = indices[current_idx];
                indices[current_idx] = current_idx;
                if indices[target_idx] == target_idx {
                    break;
                }
                data.swap(current_idx, target_idx);
                current_idx = target_idx;
            }
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
    pub namewidth: u16,  // number of graphemes of each name displayed
    pub consensus: bool, // if consensus view is shown
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
            consensus: false,
            aln,
        };
        // We need to resize before we resize names, because the latter
        // depends on a nonzero terminal size.
        view.resize(ncols, nrows);
        view.resize_names((ncols >> 2) as isize);
        view
    }

    pub fn from_reader<T: BufRead>(file: T, uppercase: bool, must_aa: bool) -> Result<View> {
        Ok(View::new(Alignment::new(file, uppercase, must_aa)?))
    }

    pub fn graphemes(&self) -> &Vec<Graphemes> {
        &self.aln.graphemes
    }

    pub fn seqs(&self) -> &Vec<Vec<u8>> {
        &self.aln.seqs
    }

    pub fn consensus(&self) -> Option<&Vec<Option<u8>>> {
        self.aln.consensus.as_ref()
    }

    pub fn calculate_consensus(&mut self) {
        self.aln.consensus = Some(calculate_consensus(&self.aln.seqs, self.is_aa()))
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

    // Returns None if operation failed, Some(()) otherwise
    pub fn move_row(&mut self, from: usize, to: usize) -> Option<()> {
        match move_element(&mut self.aln.graphemes, from, to) {
            Some(_) => {
                // If the first succeeds, the other MUST also succeed,
                // else the fields go out of synch and we must panic
                move_element(&mut self.aln.seqs, from, to).unwrap();
                if let Some(v) = &mut self.aln.consensus {
                    move_element(v, from, to).unwrap();
                }
                Some(())
            }
            None => None,
        }
    }

    pub fn resize_names(&mut self, delta: isize) {
        let mut namewidth = (self.namewidth as isize) + delta;
        namewidth = max(0, namewidth); // not negative
        namewidth = min(namewidth, (self.term_ncols as isize).saturating_sub(2)); // do not exceed bounds

        // Do not exceed longest name shown on screen
        namewidth = min(namewidth, self.aln.longest_name as isize);
        self.namewidth = namewidth.try_into().unwrap();
    }

    pub fn order(&mut self) {
        self.aln.order()
    }

    pub fn nrows(&self) -> usize {
        self.aln.nrows()
    }

    pub fn ncols(&self) -> usize {
        self.aln.ncols()
    }

    pub fn seq_nrows_display(&self) -> usize {
        (self.term_nrows as usize)
            .saturating_sub(HEADER_LINES + FOOTER_LINES + self.consensus as usize)
    }

    /// Index of last seq row
    fn last_seq_row(&self) -> Option<usize> {
        match self.seq_nrows_display() {
            0 => None,
            nrows => Some(min(self.nrows() - 1, self.rowstart + nrows - 1)),
        }
    }

    pub fn seq_row_range(&self) -> Option<RangeInclusive<usize>> {
        Some(self.rowstart..=self.last_seq_row()?)
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

    // Returns None if the row is not drawn on screen, and the u16 otherwise
    pub fn display_row_of_index(&self, index: usize) -> Option<u16> {
        match self.seq_row_range() {
            None => None,
            Some(range) => {
                if range.contains(&index) {
                    Some(
                        (index.checked_sub(self.rowstart).unwrap()
                            + HEADER_LINES
                            + self.consensus as usize) as u16,
                    )
                } else {
                    None
                }
            }
        }
    }
}
