use std::fmt;

/// Error type for translation failures
#[derive(Debug, Clone)]
pub enum TranslationError {
    /// A codon contains a mix of gaps and non-gaps (e.g., "A-G")
    /// This is invalid because the ungapped sequence would translate differently
    MixedGapCodon { position: usize },
}

impl fmt::Display for TranslationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TranslationError::MixedGapCodon { position } => {
                write!(
                    f,
                    "Mixed gap codon at position {} (gaps must align with codon boundaries)",
                    position + 1 // 1-indexed for display
                )
            }
        }
    }
}

impl std::error::Error for TranslationError {}

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
/// Returns '-' if codon is all gaps, 'X' for ambiguous/unknown bases.
/// Returns error if codon has mixed gaps (some gaps, some non-gaps).
pub fn translate_codon(codon: [u8; 3]) -> Result<u8, TranslationError> {
    // Check for gaps - must be all gaps or no gaps
    let gap_count = codon.iter().filter(|&&b| b == b'-').count();
    match gap_count {
        0 => {} // No gaps, continue to translation
        3 => return Ok(b'-'), // All gaps
        _ => return Err(TranslationError::MixedGapCodon { position: 0 }), // Mixed gaps - error
    }

    let b1 = match base_to_index(codon[0]) {
        Some(i) => i,
        None => return Ok(b'X'),
    };
    let b2 = match base_to_index(codon[1]) {
        Some(i) => i,
        None => return Ok(b'X'),
    };
    let b3 = match base_to_index(codon[2]) {
        Some(i) => i,
        None => return Ok(b'X'),
    };

    Ok(CODON_TABLE[b1 * 16 + b2 * 4 + b3])
}

/// Translate a DNA sequence to protein.
/// Assumes sequence is in-frame (starts at codon position 0).
/// Partial codons at the end are translated to 'X'.
/// Returns error if any codon contains mixed gaps.
pub fn translate_sequence(seq: &[u8]) -> Result<Vec<u8>, TranslationError> {
    let protein_len = seq.len().div_ceil(3);
    let mut protein = Vec::with_capacity(protein_len);

    // TODO: Use array_chunks when stabilized
    for (i, chunk) in seq.chunks(3).enumerate() {
        if chunk.len() == 3 {
            let codon: [u8; 3] = chunk.try_into().unwrap();
            protein.push(translate_codon(codon).map_err(|e| match e {
                TranslationError::MixedGapCodon { .. } => {
                    TranslationError::MixedGapCodon { position: i * 3 }
                }
            })?);
        } else {
            // Partial codon at end -> X
            protein.push(b'X');
        }
    }

    Ok(protein)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_codon() {
        // Standard codons
        assert_eq!(translate_codon(*b"ATG").unwrap(), b'M'); // Start codon
        assert_eq!(translate_codon(*b"TAA").unwrap(), b'*'); // Stop
        assert_eq!(translate_codon(*b"TAG").unwrap(), b'*'); // Stop
        assert_eq!(translate_codon(*b"TGA").unwrap(), b'*'); // Stop
        assert_eq!(translate_codon(*b"TTT").unwrap(), b'F'); // Phe
        assert_eq!(translate_codon(*b"GGG").unwrap(), b'G'); // Gly

        // Case insensitive
        assert_eq!(translate_codon(*b"atg").unwrap(), b'M');
        assert_eq!(translate_codon(*b"AtG").unwrap(), b'M');

        // U instead of T (RNA)
        assert_eq!(translate_codon(*b"AUG").unwrap(), b'M');

        // All gaps
        assert_eq!(translate_codon(*b"---").unwrap(), b'-');

        // Ambiguous bases
        assert_eq!(translate_codon(*b"ATN").unwrap(), b'X');
        assert_eq!(translate_codon(*b"NNN").unwrap(), b'X');
    }

    #[test]
    fn test_translate_codon_mixed_gaps() {
        // Mixed gaps should return error
        assert!(translate_codon(*b"A-G").is_err());
        assert!(translate_codon(*b"AT-").is_err());
        assert!(translate_codon(*b"-TG").is_err());
        assert!(translate_codon(*b"--G").is_err());
        assert!(translate_codon(*b"A--").is_err());
    }

    #[test]
    fn test_translate_sequence() {
        assert_eq!(translate_sequence(b"ATGTTT").unwrap(), b"MF");
        assert_eq!(translate_sequence(b"ATGTTTGGG").unwrap(), b"MFG");
        // Partial codon at end produces X
        assert_eq!(translate_sequence(b"ATGTTTT").unwrap(), b"MFX");
        assert_eq!(translate_sequence(b"ATGTTTTT").unwrap(), b"MFX");
        // All-gap codons work
        assert_eq!(translate_sequence(b"ATG---TTT").unwrap(), b"M-F");
    }

    #[test]
    fn test_translate_sequence_mixed_gaps() {
        // Mixed gaps in sequence should return error with correct position
        let result = translate_sequence(b"ATGA-GTTT");
        assert!(result.is_err());
        if let Err(TranslationError::MixedGapCodon { position }) = result {
            assert_eq!(position, 3); // Error at codon starting at position 3
        } else {
            panic!("Expected MixedGapCodon error");
        }
    }
}
