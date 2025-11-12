use anyhow::{Context, Result};
use ndarray::Array2;
use needletail::parse_fastx_file;

/// Represents a sequence record with ID and sequence data
#[derive(Debug, Clone)]
pub struct SeqRecord {
    pub id: String,
    pub seq: Vec<char>,
}

/// Reads and processes alignment files
pub struct AlignmentReader {
    file: String,
}

impl AlignmentReader {
    pub fn new(file: String) -> Self {
        Self { file }
    }

    /// Main method to read alignment and create array representation
    pub fn run(&self) -> Result<(Vec<SeqRecord>, Array2<char>)> {
        let alignment = self.read_alignment()?;
        let alignment_array = self.make_alignment_array(&alignment)?;
        Ok((alignment, alignment_array))
    }

    /// Read the alignment file using needletail
    fn read_alignment(&self) -> Result<Vec<SeqRecord>> {
        let mut records = Vec::new();
        let mut reader = parse_fastx_file(&self.file)
            .with_context(|| format!("Failed to open alignment file: {}", self.file))?;

        while let Some(record) = reader.next() {
            let rec = record.with_context(|| "Failed to parse FASTA record")?;
            let id = String::from_utf8_lossy(rec.id()).to_string();
            let seq: Vec<char> = rec.seq().iter().map(|&b| (b as char).to_ascii_uppercase()).collect();
            records.push(SeqRecord { id, seq });
        }

        if records.is_empty() {
            anyhow::bail!("No sequences found in alignment file");
        }

        Ok(records)
    }

    /// Convert alignment to 2D numpy-like array
    fn make_alignment_array(&self, alignment: &[SeqRecord]) -> Result<Array2<char>> {
        if alignment.is_empty() {
            anyhow::bail!("Cannot create array from empty alignment");
        }

        let num_seqs = alignment.len();
        let seq_len = alignment[0].seq.len();

        // Verify all sequences have the same length
        for (i, record) in alignment.iter().enumerate() {
            if record.seq.len() != seq_len {
                anyhow::bail!(
                    "Sequence {} has length {} but expected {}. All sequences must have the same length.",
                    i,
                    record.seq.len(),
                    seq_len
                );
            }
        }

        let alignment_array = Array2::from_shape_fn((num_seqs, seq_len), |(seq_idx, pos)| {
            alignment[seq_idx].seq[pos]
        });

        Ok(alignment_array)
    }
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
