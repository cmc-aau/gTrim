use anyhow::{bail, Result};
use log::warn;
use std::collections::HashSet;

use ndarray::Array2;

pub const GAP_CHARS: [u8; 2] = [b'-', b'.'];
pub const DNA_BASES: [u8; 4] = *b"ACGT";
pub const AMINO_ACIDS: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

pub struct SequenceTypeDetector;

impl SequenceTypeDetector {
    /// Detect whether sequences are DNA or amino acids
    /// Returns (is_dna, char_set)
    pub fn detect(alignment_array: &Array2<u8>) -> Result<(bool, Vec<u8>)> {
        let chars_to_remove: HashSet<u8> = b"Nn-?Xx.".iter().copied().collect();

        let dna_set: HashSet<u8> = DNA_BASES.iter().copied().collect();
        let amino_acid_set: HashSet<u8> = AMINO_ACIDS.iter().copied().collect();

        let mut seq_set: HashSet<u8> = HashSet::new();

        // Collect unique characters from first 100,000 elements of the array
        for (idx, &ch) in alignment_array.iter().enumerate() {
            seq_set.insert(ch.to_ascii_uppercase());

            if idx > 100000 { break }
        }

        // Remove ambiguous/gap characters
        seq_set.retain(|c| !chars_to_remove.contains(c));

        if seq_set.is_empty() {
            warn!("⚠️ Could not detect sequence type: all residues were ambiguous or removed.");
            bail!("No valid sequence characters detected.");
        }

        if seq_set.is_subset(&dna_set) {
            Ok((true, DNA_BASES.to_vec()))
        } else if seq_set.is_subset(&amino_acid_set) {
            Ok((false, AMINO_ACIDS.to_vec()))
        } else {
            let mut sorted_chars: Vec<char> = seq_set.iter().map(|&c| c as char).collect();
            sorted_chars.sort();
            warn!("⚠️  Sequence contains unexpected characters: {:?}", sorted_chars);
            warn!("Defaulting to amino acid interpretation.");
            Ok((false, AMINO_ACIDS.to_vec()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;

    #[test]
    fn test_dna_detection() {
        let alignment_array = arr2(&[*b"AA-G", *b"A--G", *b"CTAG"]);

        let (is_dna, _char_set) = SequenceTypeDetector::detect(&alignment_array).unwrap();
        assert!(is_dna);
    }

    #[test]
    fn test_amino_acid_detection() {
        let alignment_array = arr2(&[*b"MK-W", *b"L--T"]);

        let (is_dna, _char_set) = SequenceTypeDetector::detect(&alignment_array).unwrap();
        assert!(!is_dna);
    }
}
