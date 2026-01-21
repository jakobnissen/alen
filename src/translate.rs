/// Universal genetic code translation table.
/// Codon index: (base1 * 16) + (base2 * 4) + base3
/// where A=0, C=1, G=2, T/U=3
const CODON_TABLE: [u8; 64] = [
    b'K', b'N', b'K', b'N', // AAA, AAC, AAG, AAT
    b'T', b'T', b'T', b'T', // ACA, ACC, ACG, ACT
    b'R', b'S', b'R', b'S', // AGA, AGC, AGG, AGT
    b'I', b'I', b'M', b'I', // ATA, ATC, ATG, ATT
    b'Q', b'H', b'Q', b'H', // CAA, CAC, CAG, CAT
    b'P', b'P', b'P', b'P', // CCA, CCC, CCG, CCT
    b'R', b'R', b'R', b'R', // CGA, CGC, CGG, CGT
    b'L', b'L', b'L', b'L', // CTA, CTC, CTG, CTT
    b'E', b'D', b'E', b'D', // GAA, GAC, GAG, GAT
    b'A', b'A', b'A', b'A', // GCA, GCC, GCG, GCT
    b'G', b'G', b'G', b'G', // GGA, GGC, GGG, GGT
    b'V', b'V', b'V', b'V', // GTA, GTC, GTG, GTT
    b'*', b'Y', b'*', b'Y', // TAA, TAC, TAG, TAT (stop codons = *)
    b'S', b'S', b'S', b'S', // TCA, TCC, TCG, TCT
    b'*', b'C', b'W', b'C', // TGA, TGC, TGG, TGT (TGA = stop)
    b'L', b'F', b'L', b'F', // TTA, TTC, TTG, TTT
];

fn base_to_index(base: u8) -> Option<usize> {
    match base | 0x20 {
        // lowercase conversion
        b'a' => Some(0),
        b'c' => Some(1),
        b'g' => Some(2),
        b't' | b'u' => Some(3),
        _ => None,
    }
}

/// Translate a 3-base codon to an amino acid.
/// Returns '-' if codon contains any gap, 'X' for ambiguous/unknown bases.
pub fn translate_codon(codon: &[u8]) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }

    // Check for gaps
    if codon.iter().any(|&b| b == b'-') {
        return b'-';
    }

    let b1 = match base_to_index(codon[0]) {
        Some(i) => i,
        None => return b'X',
    };
    let b2 = match base_to_index(codon[1]) {
        Some(i) => i,
        None => return b'X',
    };
    let b3 = match base_to_index(codon[2]) {
        Some(i) => i,
        None => return b'X',
    };

    CODON_TABLE[b1 * 16 + b2 * 4 + b3]
}

/// Translate a DNA sequence to protein.
/// Assumes sequence is in-frame (starts at codon position 0).
/// Truncates to complete codons.
pub fn translate_sequence(seq: &[u8]) -> Vec<u8> {
    let protein_len = seq.len() / 3;
    let mut protein = Vec::with_capacity(protein_len);

    for i in 0..protein_len {
        let codon = &seq[i * 3..i * 3 + 3];
        protein.push(translate_codon(codon));
    }

    protein
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_codon() {
        // Standard codons
        assert_eq!(translate_codon(b"ATG"), b'M'); // Start codon
        assert_eq!(translate_codon(b"TAA"), b'*'); // Stop
        assert_eq!(translate_codon(b"TAG"), b'*'); // Stop
        assert_eq!(translate_codon(b"TGA"), b'*'); // Stop
        assert_eq!(translate_codon(b"TTT"), b'F'); // Phe
        assert_eq!(translate_codon(b"GGG"), b'G'); // Gly

        // Case insensitive
        assert_eq!(translate_codon(b"atg"), b'M');
        assert_eq!(translate_codon(b"AtG"), b'M');

        // U instead of T (RNA)
        assert_eq!(translate_codon(b"AUG"), b'M');

        // Gaps
        assert_eq!(translate_codon(b"A-G"), b'-');
        assert_eq!(translate_codon(b"---"), b'-');

        // Ambiguous bases
        assert_eq!(translate_codon(b"ATN"), b'X');
        assert_eq!(translate_codon(b"NNN"), b'X');
    }

    #[test]
    fn test_translate_sequence() {
        assert_eq!(translate_sequence(b"ATGTTT"), b"MF");
        assert_eq!(translate_sequence(b"ATGTTTGGG"), b"MFG");
        // Partial codon at end is truncated
        assert_eq!(translate_sequence(b"ATGTTTT"), b"MF");
        assert_eq!(translate_sequence(b"ATGTTTTT"), b"MF");
        // Gaps
        assert_eq!(translate_sequence(b"ATG---TTT"), b"M-F");
    }
}
