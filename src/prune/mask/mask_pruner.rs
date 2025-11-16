use anyhow::Result;
use log::info;

use crate::sequence_type_detector::SequenceTypeDetector;
use crate::alignment_reader::{AlignmentReader, read_mask};
use crate::alignment_writer::{get_output_prefix, check_compressed_input, write_filtered_alignment};

pub struct MaskPruner {
    alignment_file: String,
    mask_file: String,
    out_dir: String,
}

impl MaskPruner {
    pub fn new(
        alignment_file: String,
        mask_file: String,
        out_dir: String,
    ) -> Self {
        Self {
            alignment_file,
            mask_file,
            out_dir,
        }
    }

    pub fn run(&mut self) -> Result<()> {
        info!("Starting mask workflow.");

        info!("Reading alignment from file: {}", self.alignment_file);
        let reader = AlignmentReader::new(self.alignment_file.clone());
        let (seq_ids, alignment_array) = reader.run()?;
        info!("Alignment contains {} sequences and {} columns.", seq_ids.len(), alignment_array.ncols());

        let (is_dna, _char_set) = SequenceTypeDetector::detect(&alignment_array)?;
        info!("Sequence type detected: {}", if is_dna { "DNA" } else { "amino acid" });

        info!("Reading mask from file: {}", self.mask_file);
        let mask = read_mask(&self.mask_file)?;
        info!("Mask specifies {} columns to retain and {} columns to filter.", 
            mask.iter().filter(|&state| *state).count(),  
            mask.iter().filter(|&state| !state).count());

        // check that mask has same length as alignment
        if mask.len() != alignment_array.ncols() {
            anyhow::bail!("Mask length {} does not match alignment length {}.", mask.len(), alignment_array.ncols());
        }

        // write results to output files
        let output_prefix = get_output_prefix(&self.alignment_file, &self.out_dir)?;
        let output_suffix = "_masked.fasta";

        info!("Writing filtered alignment to file.");
        let output_alignment_file = format!("{}{}", output_prefix, output_suffix);
        let compressed = check_compressed_input(&self.alignment_file);
        write_filtered_alignment(&seq_ids, &alignment_array, &mask, &output_alignment_file, compressed)?;

        Ok(())
    }
}
