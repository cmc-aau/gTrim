use anyhow::{Context, Result};
use ndarray::Array2;
use needletail::parse_fastx_file;

/// Reads and processes alignment files
pub struct AlignmentReader {
    file: String,
}

impl AlignmentReader {
    pub fn new(file: String) -> Self {
        Self { file }
    }

    /// Main method to read alignment and create array representation
    pub fn run(&self) -> Result<(Vec<String>, Array2<char>)> {
        let (seq_ids, alignment_array) = self.read_alignment()?;
        Ok((seq_ids, alignment_array))
    }

    /// Read the alignment file using needletail
    fn read_alignment(&self) -> Result<(Vec<String>, Array2<char>)> {
        // read sequence IDs and residues
        let mut seq_ids = Vec::new();
        let mut reader = parse_fastx_file(&self.file)
            .with_context(|| format!("Failed to open alignment file: {}", self.file))?;

        let mut seqs = Vec::new();
        while let Some(record) = reader.next() {
            let rec = record.with_context(|| "Failed to parse FASTA record")?;
            let id = String::from_utf8_lossy(rec.id()).to_string();
            let seq: Vec<char> = rec.seq().iter().map(|&b| (b as char).to_ascii_uppercase()).collect();
            seq_ids.push(id);
            seqs.push(seq);
        }

        if seq_ids.is_empty() {
            anyhow::bail!("No sequences found in alignment file");
        }

        // Verify all sequences have the same length
        let seq_len = seqs[0].len();
        for (idx, seq) in seqs.iter().enumerate() {
            if seq.len() != seq_len {
                anyhow::bail!(
                    "Sequence {} has length {} but expected {}. All sequences must have the same length.",
                    idx,
                    seq.len(),
                    seq_len
                );
            }
        }

        let alignment_array = Array2::from_shape_fn((seq_ids.len(), seq_len), |(seq_idx, pos)| {
            seqs[seq_idx][pos]
        });

        Ok((seq_ids, alignment_array))
    }
}

// Read binary mask (0 = false; 1 = true).
pub fn read_mask(mask_file: &str) -> Result<Vec<bool>> {
    let mask_string = std::fs::read_to_string(mask_file)
        .with_context(|| format!("Failed to read mask file: {}", mask_file))?;

    let mask_line = mask_string.lines().next().unwrap_or_default();

    let mut mask = Vec::new();
    for char_code in mask_line.chars() {
        match char_code {
            '0' => mask.push(false),
            '1' => mask.push(true),
            _ => anyhow::bail!("Invalid character '{}' found in mask file. Only '0' and '1' are allowed.", char_code),
        }
    }

    Ok(mask)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignment_reader_basic() {
        // This test would need a test file
        // For now, just test the struct creation
        let reader = AlignmentReader::new("test.fasta".to_string());
        assert_eq!(reader.file, "test.fasta");
    }
}
