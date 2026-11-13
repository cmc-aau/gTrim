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

    /// Read the alignment into sequence IDs and an array of uppercased residues as ASCII bytes
    /// (rows = sequences, columns = alignment positions).
    pub fn run(&self) -> Result<(Vec<String>, Array2<u8>)> {
        let mut reader = parse_fastx_file(&self.file)
            .with_context(|| format!("Failed to open alignment file: {}", self.file))?;

        let mut builder = AlignmentBuilder::default();
        while let Some(record) = reader.next() {
            let rec = record.with_context(|| "Failed to parse FASTA record")?;
            builder.push(String::from_utf8_lossy(rec.id()).into_owned(), &rec.seq())?;
        }

        builder.finish()
    }
}

/// Collects sequences into a single row-major buffer that becomes the backing storage of the
/// alignment array, so residues are stored once at one byte each.
#[derive(Default)]
struct AlignmentBuilder {
    seq_ids: Vec<String>,
    residues: Vec<u8>,
    seq_len: Option<usize>,
}

impl AlignmentBuilder {
    fn push(&mut self, id: String, seq: &[u8]) -> Result<()> {
        let expected_len = *self.seq_len.get_or_insert(seq.len());
        if seq.len() != expected_len {
            anyhow::bail!(
                "Sequence {} has length {} but expected {}. All sequences must have the same length.",
                self.seq_ids.len(),
                seq.len(),
                expected_len
            );
        }

        self.residues.extend(seq.iter().map(u8::to_ascii_uppercase));
        self.seq_ids.push(id);

        Ok(())
    }

    fn finish(mut self) -> Result<(Vec<String>, Array2<u8>)> {
        let Some(seq_len) = self.seq_len else {
            anyhow::bail!("No sequences found in alignment file");
        };

        self.residues.shrink_to_fit();
        let alignment_array = Array2::from_shape_vec((self.seq_ids.len(), seq_len), self.residues)
            .context("Failed to build alignment array")?;

        Ok((self.seq_ids, alignment_array))
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
    use ndarray::arr2;

    fn read(records: &[(&str, &str)]) -> Result<(Vec<String>, Array2<u8>)> {
        let mut builder = AlignmentBuilder::default();
        for (id, seq) in records {
            builder.push(id.to_string(), seq.as_bytes())?;
        }
        builder.finish()
    }

    #[test]
    fn test_builder() {
        let (ids, alignment) = read(&[("s1", "ac-g"), ("s2", "A.tG"), ("s3", "xKn*")]).unwrap();
        assert_eq!(ids, vec!["s1", "s2", "s3"]);
        assert_eq!(alignment, arr2(&[*b"AC-G", *b"A.TG", *b"XKN*"]));
    }

    #[test]
    fn test_builder_unequal_lengths() {
        let err = read(&[("s1", "ACGT"), ("s2", "ACG")]).unwrap_err();
        assert!(err.to_string().contains("Sequence 1 has length 3 but expected 4"));
    }

    #[test]
    fn test_builder_empty() {
        assert!(read(&[]).is_err());
    }
}
