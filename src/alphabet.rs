use std::num::NonZeroU8;

pub const PRINTABLE_ASCII_OFFSET: u8 = b'!';
pub const PRINTABLE_ASCII_LEN: usize = (b'~' - b'!' + 1) as usize;

const AMINO_ACID_DETECTION_LUT: [bool; 256] = build_bool_lut(b"*EFIJLOPQXZ");
const AMINO_ACID_INDEX_LUT: [Option<NonZeroU8>; 256] = build_index_lut(b"ACDEFGHIKLMNPQRSTVWY");
const NUCLEOTIDE_INDEX_LUT: [Option<NonZeroU8>; 256] = {
    let mut lut = build_index_lut(b"ACGT");
    // Include 'U' for RNA sequences, mapping it to the same index as 'T'.
    set_index(&mut lut, b'U', 3);
    lut
};

const fn ascii_lowercase(byte: u8) -> u8 {
    if byte >= b'A' && byte <= b'Z' {
        byte + (b'a' - b'A')
    } else {
        byte
    }
}

const fn encode_index(index: u8) -> NonZeroU8 {
    match NonZeroU8::new(index.wrapping_add(1)) {
        Some(value) => value,
        None => panic!("index must fit in NonZeroU8"),
    }
}

const fn set_index(lut: &mut [Option<NonZeroU8>; 256], byte: u8, index: u8) {
    let encoded = encode_index(index);
    lut[byte as usize] = Some(encoded);
    lut[ascii_lowercase(byte) as usize] = Some(encoded);
}

const fn build_index_lut(bytes: &[u8]) -> [Option<NonZeroU8>; 256] {
    let mut lut = [None; 256];
    let mut i = 0;

    while i < bytes.len() {
        set_index(&mut lut, bytes[i], i as u8);
        i += 1;
    }

    lut
}

const fn build_bool_lut(entries: &[u8]) -> [bool; 256] {
    let mut lut = [false; 256];
    let mut i = 0;

    while i < entries.len() {
        let byte = entries[i];
        lut[byte as usize] = true;
        lut[ascii_lowercase(byte) as usize] = true;
        i += 1;
    }

    lut
}

pub fn indicates_amino_acid(byte: u8) -> bool {
    lookup_flag(&AMINO_ACID_DETECTION_LUT, byte)
}

fn lookup_index(lut: &[Option<NonZeroU8>; 256], byte: u8) -> Option<usize> {
    lut[byte as usize].map(|index| usize::from(index.get()) - 1)
}

fn lookup_flag(lut: &[bool; 256], byte: u8) -> bool {
    lut[byte as usize]
}

pub fn nucleotide_index(byte: u8) -> Option<usize> {
    lookup_index(&NUCLEOTIDE_INDEX_LUT, byte)
}

fn amino_acid_index(byte: u8) -> Option<usize> {
    lookup_index(&AMINO_ACID_INDEX_LUT, byte)
}

pub fn residue_index(byte: u8, is_aa: bool) -> Option<usize> {
    if is_aa {
        amino_acid_index(byte)
    } else {
        nucleotide_index(byte)
    }
}
