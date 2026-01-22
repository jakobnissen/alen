use crate::constants::{FOOTER_LINES, HEADER_LINES};
use crate::translate::TranslationError;

use std::cell::OnceCell;
use std::cmp::{max, min};
use std::convert::TryInto;
use std::io::BufRead;
use std::num::{NonZeroU8, NonZeroUsize};
use std::ops::RangeInclusive;

use anyhow::{Result, anyhow};

use unicode_segmentation::UnicodeSegmentation;

use crossterm::terminal;

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

// Bitmap of ASCII values 1 << (x - b'E')
// since first letter only in AA alphabet is E.
// Asterisk (*) is also AA only, but we check that separately
const ONLY_AA: u64 = 0x00281cb300281cb3;
const ONLY_AA_OFFSET: u8 = b'E';

// Returns whether it's an AA alphabet, else it default to DNA.
fn verify_alphabet(entries: &[Entry], must_aa: bool) -> Result<bool> {
    let mut must_be_aa = must_aa;
    for entry in entries.iter() {
        for &byte in entry.seq.iter() {
            // Printable ASCII range: ! to ~
            if !(33..=126).contains(&byte) {
                return Err(anyhow!(
                    "Sequence named \"{}\" contains non-ASCII-printable byte 0x{:x}",
                    entry.graphemes.string,
                    byte
                ));
            }
            // It's asterisk, or else it's both at least ONLY_AA_OFFSET (meaning)
            // at least the value of the first only-AA character, and also found in the ONLY_AA bitmap.
            must_be_aa |= (byte == b'*')
                | ((ONLY_AA_OFFSET..=b'z').contains(&byte)
                    & ((ONLY_AA.wrapping_shr(byte.wrapping_sub(ONLY_AA_OFFSET).into()) & 1) == 1));
        }
    }
    Ok(must_be_aa)
}

fn calculate_consensus<'a, T: Iterator<Item = &'a Vec<u8>>>(
    seqs: T,
    ncols: usize,
    is_aa: bool,
) -> Vec<Option<NonZeroU8>> {
    // We have verified seq is between ASCII 33 and 126, inclusive.
    let offset = b'!';
    let nonzero_offset = NonZeroU8::new(offset).unwrap();
    let mut counts = vec![[0u32; 126 - 33 + 1]; ncols];

    // First loop over sequences in memory order
    for seq in seqs {
        for (byte, arr) in seq.iter().zip(counts.iter_mut()) {
            arr[(byte - offset) as usize] += 1
        }
    }

    for arr in counts.iter_mut() {
        // Uppercase all characters
        for i in b'a'..=b'z' {
            arr[(i - offset - 32) as usize] += arr[(i - offset) as usize];
            arr[(i - offset) as usize] = 0
        }

        // We set all ambiguous bases/AAs to 0.
        // These are useless in the consensus
        if is_aa {
            for &i in b"BJOUXZ" {
                arr[(i - offset) as usize] = 0;
            }
        } else {
            for &i in b"MRSVWYHKDBN" {
                arr[(i - offset) as usize] = 0;
            }
        }
    }

    counts
        .iter()
        .map(|arr| {
            let (index, count) = arr.iter().enumerate().max_by_key(|&(_, x)| x).unwrap();

            if *count == 0 {
                None
            } else {
                Some(nonzero_offset.saturating_add(index as u8))
            }
        })
        .collect()
}

fn move_element<T>(v: &mut [T], from: usize, to: usize) -> Option<()> {
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

struct Entry {
    graphemes: Graphemes,
    seq: Vec<u8>,
    original_index: usize,
    // first_nogap..last_nogap into seq gets non-gapped seq
    first_nogap: u32,
    last_nogap: u32,
}

pub struct Alignment {
    entries: Vec<Entry>,
    // longest as in number of graphemes. We cache this for efficiency, it can be
    // computed from the graphemes field easily
    longest_name: usize,
    // we calculate this lazily upon demand (for nucleotide view)
    consensus: OnceCell<Vec<Option<NonZeroU8>>>,
    // consensus for translated protein view (separate cache)
    translated_consensus: OnceCell<Vec<Option<NonZeroU8>>>,
    // When sorting entries by |ent| order[ent.original_index], the rows are ordered.
    // also calculate this lazily
    order: OnceCell<Vec<u32>>,
    is_aa: bool,
    // Translated protein sequences (lazily computed, Result for error handling)
    translated: OnceCell<Result<Vec<Vec<u8>>, TranslationError>>,
}

impl Alignment {
    fn nrows(&self) -> usize {
        self.entries.len()
    }

    fn ncols(&self) -> usize {
        self.entries[0].seq.len()
    }

    fn new<T: BufRead>(file: T, uppercase: bool, must_aa: bool) -> Result<Alignment> {
        let reader = fasta::Reader::new(file);
        let mut seqlength: Option<NonZeroUsize> = None;
        let mut entries = Vec::new();

        for (original_index, result) in reader.records().enumerate() {
            // Turn uppercase
            // Check alphabet
            // Chech sequence length
            let record = result?;
            let graphemes = Graphemes::new(record.id());
            let seq = record.seq().to_vec();
            let this_seq_len = match NonZeroUsize::new(seq.len()) {
                None => {
                    return Err(anyhow!(
                        "Sequence \"{}\" has length zero, which is not allowed.",
                        graphemes.string
                    ));
                }
                Some(len) => len,
            };

            // Check identical sequence lengths
            if let Some(len) = seqlength {
                if len != this_seq_len {
                    return Err(anyhow!(
                        "Sequence \"{}\" has a different length than the previous sequence. \
                        In an alignment, all sequences must have the same length.",
                        graphemes.string
                    ));
                }
            } else {
                seqlength = Some(this_seq_len)
            }

            // start..stop is span of non-deleted symbols
            let (first_nogap, last_nogap) = match seq.iter().position(|&i| i != b'-') {
                None => (0, 0),
                Some(u) => {
                    let last_nogap: u32 =
                        (seq.len() - seq.iter().rev().position(|&i| i != b'-').unwrap()) as u32;
                    (u as u32, last_nogap)
                }
            };

            entries.push(Entry {
                graphemes,
                seq,
                original_index,
                first_nogap,
                last_nogap,
            })
        }

        // Turn uppercase if requested
        if uppercase {
            entries
                .iter_mut()
                .for_each(|s| s.seq.make_ascii_uppercase());
        }

        // Verify alphabet
        let is_aa = verify_alphabet(&entries, must_aa)?;

        let longest_name = entries.iter().map(|v| v.graphemes.len()).max().unwrap();
        Ok(Alignment {
            entries,
            longest_name,
            consensus: OnceCell::new(),
            translated_consensus: OnceCell::new(),
            order: OnceCell::new(),
            is_aa,
            translated: OnceCell::new(),
        })
    }

    /// Reorder the vectors of the alignment such that similar rows are next to each other.
    fn compute_order(&self) -> Vec<u32> {
        // We do this in the caller, to simplify the logic here, so we don't have
        // to handle differing indices and original indices
        assert!(self.entries.is_sorted_by_key(|e| e.original_index));

        // If already ordered or 2 or fewer rows, ordering doesn't matter
        if self.nrows() < 3 {
            return (0..(self.nrows().try_into().unwrap())).collect();
        }

        // Choices only appear when placing the 3rd seq, so first two are given.
        let mut order: Vec<usize> = vec![0, 1];
        let mut neighbor_distances = vec![jaccard_distance(&self.entries[0], &self.entries[1])];
        let mut distances: Vec<f32> = Vec::with_capacity(self.nrows());

        for seqindex in 2..self.nrows() {
            // Jaccard distances from seq[seqindex] to all that are already placed.
            distances.clear();
            distances.extend(
                order
                    .iter()
                    .map(|i| jaccard_distance(&self.entries[seqindex], &self.entries[*i])),
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

        // Now transform order such that if the order begins with [8, 3, 0], then the 8th
        // index is 0, 3th index is 1, 0th index is 2 etc.
        // That means we can sort the entries by looking up directly in the order
        let mut ord: Vec<u32> = vec![u32::MAX; self.nrows()];
        for (i, o) in order.iter().enumerate() {
            // Check each element in order is unique
            assert!(ord[*o] == u32::MAX);
            ord[*o] = i.try_into().unwrap()
        }
        ord
    }

    fn order(&mut self) {
        // Sort by original index. This makes `compute_order` simpler, because
        // we don't have to handle the case when observed and original index
        // are different
        // We sort in O(n) here
        for i in 0..self.entries.len() {
            let mut dst = self.entries[i].original_index;
            while self.entries[i].original_index != i {
                self.entries.swap(i, dst);
                dst = self.entries[i].original_index;
            }
        }

        let order = self.order.get_or_init(|| self.compute_order());

        // Now sort by order in O(n)
        for i in 0..self.entries.len() {
            let mut dst = order[self.entries[i].original_index] as usize;
            while order[self.entries[i].original_index] as usize != i {
                self.entries.swap(i, dst);
                dst = order[self.entries[i].original_index] as usize;
            }
        }
    }
}

fn jaccard_distance(a: &Entry, b: &Entry) -> f32 {
    let start = a.first_nogap.max(b.first_nogap) as usize;
    let stop = a.last_nogap.min(b.last_nogap) as usize;
    if stop <= start {
        return 1.0;
    }
    let d: usize = a.seq[start..stop]
        .iter()
        .zip(b.seq[start..stop].iter())
        .map(|(i, j)| (i != j) as usize)
        .sum();
    ((d as f64) / (stop - start) as f64) as f32
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
    pub namewidth: u16,   // number of graphemes of each name displayed
    pub consensus: bool,  // if consensus view is shown
    pub translated: bool, // if showing translated protein view (DNA only)
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
            translated: false,
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

    pub fn graphemes(&self) -> impl Iterator<Item = &Graphemes> {
        self.aln.entries.iter().map(|e| &e.graphemes)
    }

    pub fn grapheme(&self, n: usize) -> Option<&Graphemes> {
        self.aln.entries.get(n).map(|x| &x.graphemes)
    }

    pub fn seqs(&self) -> impl Iterator<Item = &Vec<u8>> {
        self.aln.entries.iter().map(|e| &e.seq)
    }

    pub fn seq(&self, n: usize) -> Option<&Vec<u8>> {
        self.aln.entries.get(n).map(|x| &x.seq)
    }

    /// Get the consensus sequence for the current view (nucleotide or translated).
    /// Computes lazily on first access for each mode.
    pub fn consensus(&self) -> &Vec<Option<NonZeroU8>> {
        if self.translated {
            self.aln.translated_consensus.get_or_init(|| {
                // Need to get translated sequences first
                let seqs = self.translated_seqs().unwrap();
                calculate_consensus(&mut seqs.iter(), seqs[0].len(), true)
            })
        } else {
            self.aln.consensus.get_or_init(|| {
                calculate_consensus(&mut self.seqs(), self.aln.ncols(), self.aln.is_aa)
            })
        }
    }

    /// Check if consensus has been computed for current view mode
    pub fn consensus_computed(&self) -> bool {
        if self.translated {
            self.aln.translated_consensus.get().is_some()
        } else {
            self.aln.consensus.get().is_some()
        }
    }

    pub fn is_aa(&self) -> bool {
        self.aln.is_aa
    }

    /// Check if source is nucleotide (translation only valid for DNA / RNA)
    pub fn can_translate(&self) -> bool {
        !self.aln.is_aa
    }

    /// Check if translation has been computed
    pub fn translation_computed(&self) -> bool {
        self.aln.translated.get().is_some()
    }

    /// Get translated sequences (computed lazily).
    /// Returns Ok(&sequences) on success, Err(&error) if translation failed.
    pub fn translated_seqs(&self) -> Result<&Vec<Vec<u8>>, &TranslationError> {
        self.aln
            .translated
            .get_or_init(|| {
                self.aln
                    .entries
                    .iter()
                    .map(|e| crate::translate::translate_sequence(&e.seq))
                    .collect::<Result<Vec<_>, _>>()
            })
            .as_ref()
    }

    /// Get a single translated sequence by index.
    /// Panics on out-of-bounds (programmer error).
    pub fn translated_seq(&self, n: usize) -> Option<&Vec<u8>> {
        let seqs = self.translated_seqs().ok()?;
        Some(&seqs[n])
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
        self.rowstart = calculate_start(self.rowstart, dy, self.seq_nrows_display(), self.nrows());
        self.colstart = calculate_start(self.colstart, dx, self.seq_ncols_display(), self.ncols());
    }

    // Returns None if operation failed, Some(()) otherwise
    pub fn move_row(&mut self, from: usize, to: usize) -> Option<()> {
        move_element(&mut self.aln.entries, from, to)
    }

    pub fn resize_names(&mut self, delta: isize) {
        let mut namewidth = (self.namewidth as isize) + delta;
        namewidth = max(0, namewidth); // not negative
        namewidth = min(namewidth, (self.term_ncols as isize).saturating_sub(2)); // do not exceed bounds

        // Do not exceed longest name shown on screen
        namewidth = min(namewidth, self.aln.longest_name as isize);
        self.namewidth = namewidth.try_into().unwrap();
    }

    pub fn order_original(&mut self) {
        self.aln
            .entries
            .sort_unstable_by(|a, b| a.original_index.cmp(&b.original_index))
    }

    pub fn order(&mut self) {
        self.aln.order()
    }

    pub fn nrows(&self) -> usize {
        self.aln.nrows()
    }

    pub fn ncols(&self) -> usize {
        if self.translated {
            self.aln.ncols().div_ceil(3)
        } else {
            self.aln.ncols()
        }
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
