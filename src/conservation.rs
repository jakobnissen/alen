use crate::alphabet::residue_index;

const MAX_ALPHABET_SIZE: usize = 20;

fn alphabet_size(is_aa: bool) -> usize {
    if is_aa { 20 } else { 4 }
}

pub fn max_conservation(is_aa: bool) -> f64 {
    (alphabet_size(is_aa) as f64).log2()
}

fn compute_sequence_conservation(
    counts: &[u32],
    total_non_gap: usize,
    num_sequences: usize,
    is_aa: bool,
) -> f64 {
    // Sequence conservation in bits, scaled by occupancy and clamped at zero.
    if total_non_gap == 0 {
        return 0.0;
    }

    let alphabet_size = alphabet_size(is_aa);
    let total_non_gap = total_non_gap as f64;
    let num_sequences = num_sequences as f64;
    let max_bits = max_conservation(is_aa);

    let entropy = counts
        .iter()
        .take(alphabet_size)
        .filter(|&&count| count > 0)
        .map(|&count| {
            let p = count as f64 / total_non_gap;
            -p * p.log2()
        })
        .sum::<f64>();

    // Small-sample correction so short columns do not overstate conservation.
    let en = (alphabet_size as f64 - 1.0) / (2.0 * total_non_gap * std::f64::consts::LN_2);

    let information = (max_bits - (entropy + en)).max(0.0);
    // Reduce sequence conservation for sparse, partly gapped columns.
    let occupancy = total_non_gap / num_sequences;
    occupancy * information
}

pub fn compute_conservation_profile<'a, I>(seqs: I, ncols: usize, is_aa: bool) -> Vec<f64>
where
    I: IntoIterator<Item = &'a [u8]>,
{
    let mut counts = vec![[0u32; MAX_ALPHABET_SIZE]; ncols];
    let mut total_non_gap = vec![0usize; ncols];
    let mut num_sequences = 0usize;

    for seq in seqs {
        num_sequences += 1;
        for (&byte, (column_counts, column_total)) in seq
            .iter()
            .zip(counts.iter_mut().zip(total_non_gap.iter_mut()))
        {
            if let Some(index) = residue_index(byte, is_aa) {
                column_counts[index] += 1;
                *column_total += 1;
            }
        }
    }

    counts
        .iter()
        .zip(total_non_gap.iter())
        .map(|(column_counts, &column_total)| {
            compute_sequence_conservation(column_counts, column_total, num_sequences, is_aa)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-12;

    #[test]
    fn test_compute_sequence_conservation() {
        let conservation = compute_sequence_conservation(&[4, 0, 0, 0], 4, 4, false);
        let expected = 2.0 - 3.0 / (8.0 * std::f64::consts::LN_2);
        assert!((conservation - expected).abs() < EPSILON);

        let conservation = compute_sequence_conservation(&[1, 1, 0, 0], 2, 2, false);
        assert!((conservation - 0.0).abs() < EPSILON);

        let conservation = compute_sequence_conservation(&[0, 0, 0, 0], 0, 4, false);
        assert!((conservation - 0.0).abs() < EPSILON);
    }

    #[test]
    fn test_compute_conservation_profile_nucleotide_rules() {
        let seqs = [
            b"A".as_slice(),
            b"-".as_slice(),
            b"N".as_slice(),
            b"a".as_slice(),
        ];
        let profile = compute_conservation_profile(seqs.iter().copied(), 1, false);
        let expected = 0.5 * (2.0 - 3.0 / (4.0 * std::f64::consts::LN_2));
        assert!((profile[0] - expected).abs() < EPSILON);

        let seqs = [b"-".as_slice(), b".".as_slice(), b"N".as_slice()];
        let profile = compute_conservation_profile(seqs.iter().copied(), 1, false);
        assert!((profile[0] - 0.0).abs() < EPSILON);

        let mixed_case = [b"aU".as_slice(), b"At".as_slice(), b"AT".as_slice()];
        let uppercase = [b"AT".as_slice(), b"AT".as_slice(), b"AT".as_slice()];

        let mixed_case_profile = compute_conservation_profile(mixed_case.iter().copied(), 2, false);
        let uppercase_profile = compute_conservation_profile(uppercase.iter().copied(), 2, false);

        assert_eq!(mixed_case_profile, uppercase_profile);
    }

    #[test]
    fn test_compute_conservation_profile_amino_acid_mode() {
        let seqs = [
            b"ME".as_slice(),
            b"ME".as_slice(),
            b"MM".as_slice(),
            b"MM".as_slice(),
        ];
        let profile = compute_conservation_profile(seqs.iter().copied(), 2, true);
        let expected = max_conservation(true) - 19.0 / (8.0 * std::f64::consts::LN_2);

        assert!((profile[0] - expected).abs() < EPSILON);
        assert!((profile[1] - 0.0).abs() < EPSILON);
    }
}
