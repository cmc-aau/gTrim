use std::fs::File;
use std::io::Write;

use anyhow::{Context, Result};
use log::info;
use ndarray::Array2;

use crate::sequence_type_detector::{GAP_CHARS, SequenceTypeDetector};
use crate::alignment_reader::AlignmentReader;
use crate::column_counts::map_column_counts;
use crate::alignment_writer::{get_output_prefix, check_compressed_input, write_filtered_alignment, write_mask};

#[derive(Default, Clone)]
struct PerColStats {
    num_gaps: usize,
    gap_perc: f64,
}
pub struct GapPruner {
    alignment_file: String,
    gap_threshold: f64,
    max_residue_pruned: usize,
    out_dir: String,
}

impl GapPruner {
    pub fn new(
        alignment_file: String,
        gap_threshold: f64,
        max_residue_pruned: usize,
        out_dir: String,
    ) -> Self {
        Self {
            alignment_file,
            gap_threshold,
            max_residue_pruned,
            out_dir,
        }
    }

    pub fn run(&mut self) -> Result<()> {
        info!("Starting gap pruning workflow.");

        info!("Reading alignment from file: {}", self.alignment_file);
        let reader = AlignmentReader::new(self.alignment_file.clone());
        let (seq_ids, alignment_array) = reader.run()?;
        info!("Alignment contains {} sequences and {} columns.", seq_ids.len(), alignment_array.ncols());

        let (is_dna, _char_set) = SequenceTypeDetector::detect(&alignment_array)?;
        info!("Sequence type detected: {}", if is_dna { "DNA" } else { "amino acid" });

        // determine columns exceeding gap threshold that should be filtered
        info!("Pruning gappy columns.");
        let col_stats = self.calculate_gap_stats(&alignment_array);
        let mask = self.filter_gappy_columns(&col_stats, self.gap_threshold);
        info!("Pruned {} columns.", mask.iter().filter(|&state| !state).count());

        // write results to output files
        let output_prefix = get_output_prefix(&self.alignment_file, &self.out_dir)?;
        let output_suffix = format!("_gt{}_m{}_pruned.fasta", self.gap_threshold, self.max_residue_pruned);

        let output_stats_file = output_suffix.replace(".fasta", ".tsv");
        let output_stats_file = format!("{}{}", output_prefix, output_stats_file);
        self.write_gap_stats(&col_stats, &mask, &output_stats_file)?;

        let output_mask_file = output_suffix.replace(".fasta", ".mask");
        let output_mask_file = format!("{}{}", output_prefix, output_mask_file);
        write_mask(&mask, &output_mask_file)?;

        info!("Writing filtered alignment to file.");
        let output_alignment_file = format!("{}{}", output_prefix, output_suffix);
        let compressed = check_compressed_input(&self.alignment_file);
        write_filtered_alignment(&seq_ids, &alignment_array, &mask, &output_alignment_file, compressed)?;

        Ok(())
    }

    /// Determine gap statistics for each column in alignment.
    fn calculate_gap_stats(&self, alignment_array: &Array2<u8>) -> Vec<PerColStats> {
        let num_rows = alignment_array.nrows() as f64;
        map_column_counts(alignment_array, |counts| {
            let num_gaps: usize = GAP_CHARS.iter().map(|&g| counts[g as usize] as usize).sum();
            let gap_perc = (num_gaps as f64 / num_rows) * 100.0;
            PerColStats { num_gaps, gap_perc }
        })
    }

    /// Filter columns that exceed the gap threshold.
    fn filter_gappy_columns(&self, gap_stats: &[PerColStats], gap_threshold: f64) -> Vec<bool> {
        // Build boolean mask: true = keep column, false = remove
        gap_stats
            .iter()
            .map(|s| s.gap_perc <= gap_threshold)
            .collect()
    }

    /// Write gap statistics.
    fn write_gap_stats(&self, col_stats: &[PerColStats], mask: &[bool], out_file: &str) -> Result<()> {
        let mut file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file))?;

        writeln!(file, "column_idx\tnum_gaps\tgap_perc\tretained")?;
        for (col_idx, cur_col_stats) in col_stats.iter().enumerate() {
            writeln!(file, "{}\t{}\t{}\t{}", 
                            col_idx, 
                            cur_col_stats.num_gaps, 
                            cur_col_stats.gap_perc,
                            mask[col_idx])?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;

    fn create_pruner() -> GapPruner {
        GapPruner::new("".to_string(), 0.0, 0, "".to_string())
    }

    #[test]
    fn test_calculate_gap_statistics() {
        let alignment_array = arr2(&[*b"AA-G", *b"A--G", *b"CTAG", *b"-TCG"]);
        let pruner = create_pruner();
        let stats = pruner.calculate_gap_stats(&alignment_array);

        assert_eq!(stats.len(), 4);

        assert_eq!(stats[0].num_gaps, 1);
        assert_eq!(stats[0].gap_perc, 25.0);

        assert_eq!(stats[1].num_gaps, 1);
        assert_eq!(stats[1].gap_perc, 25.0);

        assert_eq!(stats[2].num_gaps, 2);
        assert_eq!(stats[2].gap_perc, 50.0);

        assert_eq!(stats[3].num_gaps, 0);
        assert_eq!(stats[3].gap_perc, 0.0);
    }

    #[test]
    fn test_filter_gappy_columns() {
        let stats = vec![
            PerColStats { num_gaps: 1, gap_perc: 25.0 },
            PerColStats { num_gaps: 2, gap_perc: 50.0 },
            PerColStats { num_gaps: 3, gap_perc: 75.0 },
            PerColStats { num_gaps: 0, gap_perc: 0.0 },
        ];
        let pruner = create_pruner();

        let mask1 = pruner.filter_gappy_columns(&stats, 50.0);
        assert_eq!(mask1, vec![true, true, false, true]);

        let mask2 = pruner.filter_gappy_columns(&stats, 24.0);
        assert_eq!(mask2, vec![false, false, false, true]);

        let mask3 = pruner.filter_gappy_columns(&stats, 100.0);
        assert_eq!(mask3, vec![true, true, true, true]);
        
        let mask4 = pruner.filter_gappy_columns(&stats, 0.0);
        assert_eq!(mask4, vec![false, false, false, true]);
    }
}