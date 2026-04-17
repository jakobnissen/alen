pub const PRINTABLE_ASCII_OFFSET: u8 = b'!';
pub const PRINTABLE_ASCII_LEN: usize = (b'~' - b'!' + 1) as usize;
// Keep U out of this set so RNA alignments do not auto-detect as amino acids.
const AMINO_ACID_DETECTION_BYTES: &[u8; 11] = b"*EFIJLOPQXZ";

pub const fn is_printable_ascii(byte: u8) -> bool {
    byte >= b'!' && byte <= b'~'
}

pub fn indicates_amino_acid(byte: u8) -> bool {
    AMINO_ACID_DETECTION_BYTES.contains(&byte.to_ascii_uppercase())
}

pub fn nucleotide_index(byte: u8) -> Option<usize> {
    match byte.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

fn amino_acid_index(byte: u8) -> Option<usize> {
    match byte.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'D' => Some(2),
        b'E' => Some(3),
        b'F' => Some(4),
        b'G' => Some(5),
        b'H' => Some(6),
        b'I' => Some(7),
        b'K' => Some(8),
        b'L' => Some(9),
        b'M' => Some(10),
        b'N' => Some(11),
        b'P' => Some(12),
        b'Q' => Some(13),
        b'R' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'V' => Some(17),
        b'W' => Some(18),
        b'Y' => Some(19),
        _ => None,
    }
}

pub fn residue_index(byte: u8, is_aa: bool) -> Option<usize> {
    if is_aa {
        amino_acid_index(byte)
    } else {
        nucleotide_index(byte)
    }
}
