use crate::alphabet::{
    AMINO_ACID_ALPHABET_SIZE, Alphabet, NUCLEOTIDE_ALPHABET_SIZE, amino_acid_index,
    nucleotide_index,
};

pub fn max_conservation(alphabet: Alphabet) -> f64 {
    (alphabet.size() as f64).log2()
}

fn compute_sequence_conservation<const N: usize>(
    counts: &[u32; N],
    total_non_gap: usize,
    num_sequences: usize,
) -> f64 {
    // Sequence conservation in bits, scaled by occupancy and clamped at zero.
    if total_non_gap == 0 {
        return 0.0;
    }

    let total_non_gap = total_non_gap as f64;
    let num_sequences = num_sequences as f64;
    let max_bits = (N as f64).log2();

    let entropy = counts
        .iter()
        .filter(|&&count| count > 0)
        .map(|&count| {
            let p = count as f64 / total_non_gap;
            -p * p.log2()
        })
        .sum::<f64>();

    // Small-sample correction so short columns do not overstate conservation.
    let en = (N as f64 - 1.0) / (2.0 * total_non_gap * std::f64::consts::LN_2);

    let information = (max_bits - (entropy + en)).max(0.0);
    // Reduce sequence conservation for sparse, partly gapped columns.
    let occupancy = total_non_gap / num_sequences;
    occupancy * information
}

fn compute_conservation_profile_inner<'a, I, const N: usize, F>(
    seqs: I,
    ncols: usize,
    residue_index: F,
) -> Vec<f64>
where
    I: IntoIterator<Item = &'a [u8]>,
    F: Fn(u8) -> Option<usize> + Copy,
{
    let mut counts = vec![[0u32; N]; ncols];
    let mut total_non_gap = vec![0usize; ncols];
    let mut num_sequences = 0usize;

    for seq in seqs {
        num_sequences += 1;
        for (&byte, (column_counts, column_total)) in seq
            .iter()
            .zip(counts.iter_mut().zip(total_non_gap.iter_mut()))
        {
            if let Some(index) = residue_index(byte) {
                column_counts[index] += 1;
                *column_total += 1;
            }
        }
    }

    counts
        .iter()
        .zip(total_non_gap.iter())
        .map(|(column_counts, &column_total)| {
            compute_sequence_conservation(column_counts, column_total, num_sequences)
        })
        .collect()
}

pub fn compute_conservation_profile<'a, I>(seqs: I, ncols: usize, alphabet: Alphabet) -> Vec<f64>
where
    I: IntoIterator<Item = &'a [u8]>,
{
    match alphabet {
        Alphabet::AminoAcid => {
            compute_conservation_profile_inner::<_, AMINO_ACID_ALPHABET_SIZE, _>(
                seqs,
                ncols,
                amino_acid_index,
            )
        }
        Alphabet::Nucleotide => {
            compute_conservation_profile_inner::<_, NUCLEOTIDE_ALPHABET_SIZE, _>(
                seqs,
                ncols,
                nucleotide_index,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Alphabet;

    const EPSILON: f64 = 1e-12;

    #[test]
    fn test_compute_sequence_conservation() {
        let conservation = compute_sequence_conservation(&[4, 0, 0, 0], 4, 4);
        let expected = 2.0 - 3.0 / (8.0 * std::f64::consts::LN_2);
        assert!((conservation - expected).abs() < EPSILON);

        let conservation = compute_sequence_conservation(&[1, 1, 0, 0], 2, 2);
        assert!((conservation - 0.0).abs() < EPSILON);

        let conservation = compute_sequence_conservation(&[0, 0, 0, 0], 0, 4);
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
        let profile = compute_conservation_profile(seqs.iter().copied(), 1, Alphabet::Nucleotide);
        let expected = 0.5 * (2.0 - 3.0 / (4.0 * std::f64::consts::LN_2));
        assert!((profile[0] - expected).abs() < EPSILON);

        let seqs = [b"-".as_slice(), b".".as_slice(), b"N".as_slice()];
        let profile = compute_conservation_profile(seqs.iter().copied(), 1, Alphabet::Nucleotide);
        assert!((profile[0] - 0.0).abs() < EPSILON);

        let mixed_case = [b"aU".as_slice(), b"At".as_slice(), b"AT".as_slice()];
        let uppercase = [b"AT".as_slice(), b"AT".as_slice(), b"AT".as_slice()];

        let mixed_case_profile =
            compute_conservation_profile(mixed_case.iter().copied(), 2, Alphabet::Nucleotide);
        let uppercase_profile =
            compute_conservation_profile(uppercase.iter().copied(), 2, Alphabet::Nucleotide);

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
        let profile = compute_conservation_profile(seqs.iter().copied(), 2, Alphabet::AminoAcid);
        let expected =
            max_conservation(Alphabet::AminoAcid) - 19.0 / (8.0 * std::f64::consts::LN_2);

        assert!((profile[0] - expected).abs() < EPSILON);
        assert!((profile[1] - 0.0).abs() < EPSILON);
    }
}
