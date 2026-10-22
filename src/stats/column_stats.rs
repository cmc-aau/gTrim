use std::fs::File;
use std::io::Write;
use std::collections::{HashSet, HashMap};

use anyhow::{Context, Result};
use log::info;
use ndarray::Array2;

use crate::sequence_type_detector::{GAP_CHARS, SequenceTypeDetector};
use crate::alignment_reader::AlignmentReader;
use crate::alignment_writer::get_output_prefix;

#[derive(Default, Clone)]
struct PerColStats {
    gap_count: usize,
    residue_counts: HashMap<char, usize>,
    column_len: usize,
}
pub struct ColumnStats {
    alignment_file: String,
    out_dir: String,
}

impl ColumnStats {
    pub fn new(
        alignment_file: String,
        out_dir: String,
    ) -> Self {
        Self {
            alignment_file,
            out_dir,
        }
    }

    pub fn run(&mut self) -> Result<()> {
        info!("Starting column stats workflow.");

        info!("Reading alignment from file: {}", self.alignment_file);
        let reader = AlignmentReader::new(self.alignment_file.clone());
        let (seq_ids, alignment_array) = reader.run()?;
        info!("Alignment contains {} sequences and {} columns.", seq_ids.len(), alignment_array.ncols());

        let (is_dna, _char_set) = SequenceTypeDetector::detect(&alignment_array)?;
        info!("Sequence type detected: {}", if is_dna { "DNA" } else { "amino acid" });

        // calculate columns statistics
        let col_stats = self.calculate_columns_stats(&alignment_array);

        // write results to output files
        let output_prefix = get_output_prefix(&self.alignment_file, &self.out_dir)?;
        let output_suffix = "_column_stats.tsv";

        let output_stats_file = format!("{}{}", output_prefix, output_suffix);
        self.write_column_stats(&col_stats, &output_stats_file)?;

        Ok(())
    }

    /// Determine gap statistics for each column in alignment.
    fn calculate_columns_stats(&self, alignment_array: &Array2<char>) -> Vec<PerColStats> {
        let gaps: HashSet<char> = HashSet::from_iter(GAP_CHARS.iter().copied());
        let mut gap_stats = Vec::with_capacity(alignment_array.ncols());
        for col_idx in 0..alignment_array.ncols() {
            let column = alignment_array.column(col_idx);
            
            let mut gap_count = 0;
            let mut residue_counts = HashMap::new();
            for ch in column.iter() {
                if gaps.contains(ch) {
                    gap_count += 1;
                } else {
                    *residue_counts.entry(*ch).or_insert(0) += 1;
                }
            }

            gap_stats.push(PerColStats { 
                gap_count, 
                residue_counts,
                column_len: column.len(),
            });
        }
        
        gap_stats
    }

    /// Write column statistics.
    fn write_column_stats(&self, col_stats: &[PerColStats], out_file: &str) -> Result<()> {
        let mut file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file))?;

        writeln!(file, "column_idx\tnum_gaps\tgap_perc\tnum_residues\tresidue_perc\tper_residue_counts\tper_residue_perc")?;
        for (col_idx, cur_col_stats) in col_stats.iter().enumerate() {
            let gap_perc = (cur_col_stats.gap_count as f64 / cur_col_stats.column_len as f64) * 100.0;
            let num_residues = cur_col_stats.column_len - cur_col_stats.gap_count;
            let residue_perc = (num_residues as f64 / cur_col_stats.column_len as f64) * 100.0;

            let mut sorted_residues: Vec<(&char, &usize)> = cur_col_stats.residue_counts.iter().collect();
            sorted_residues.sort_by(|a, b| b.1.cmp(a.1));

            let residue_counts_str = sorted_residues
                .iter()
                .map(|(res, count)| format!("{}:{}", res, count))
                .collect::<Vec<String>>()
                .join(",");

            let residue_perc_str = sorted_residues
                .into_iter()
                .map(|(res, count)| format!("{}:{:.2}", res, (*count as f64 / cur_col_stats.column_len as f64) * 100.0))
                .collect::<Vec<String>>()
                .join(",");

            writeln!(file, "{}\t{}\t{:.2}\t{}\t{:.2}\t{}\t{}", 
                            col_idx, 
                            cur_col_stats.gap_count, 
                            gap_perc,
                            num_residues,
                            residue_perc,
                            residue_counts_str,
                            residue_perc_str)?;
        }

        Ok(())
    }
}
