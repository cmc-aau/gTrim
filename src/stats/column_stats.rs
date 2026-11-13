use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use log::info;
use ndarray::Array2;

use crate::sequence_type_detector::{GAP_CHARS, SequenceTypeDetector};
use crate::alignment_reader::AlignmentReader;
use crate::column_counts::map_column_counts;
use crate::alignment_writer::get_output_prefix;

#[derive(Default, Clone)]
struct PerColStats {
    gap_count: usize,
    /// Non-gap residues and their counts, most common first (ties in character code order).
    residue_counts: Vec<(u8, usize)>,
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

    /// Determine gap and residue statistics for each column in alignment.
    fn calculate_columns_stats(&self, alignment_array: &Array2<u8>) -> Vec<PerColStats> {
        let column_len = alignment_array.nrows();
        map_column_counts(alignment_array, |counts| {
            let mut gap_count = 0;
            let mut residue_counts = Vec::new();
            for (code, &n) in counts.iter().enumerate() {
                let code = code as u8;
                if GAP_CHARS.contains(&code) {
                    gap_count += n as usize;
                } else if n > 0 {
                    residue_counts.push((code, n as usize));
                }
            }
            residue_counts.sort_by(|a, b| b.1.cmp(&a.1));

            PerColStats {
                gap_count,
                residue_counts,
                column_len,
            }
        })
    }

    /// Write column statistics.
    fn write_column_stats(&self, col_stats: &[PerColStats], out_file: &str) -> Result<()> {
        let file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file))?;
        let mut file = BufWriter::new(file);

        writeln!(file, "column_idx\tnum_gaps\tgap_perc\tnum_residues\tresidue_perc\tper_residue_counts\tper_residue_perc")?;
        for (col_idx, cur_col_stats) in col_stats.iter().enumerate() {
            let gap_perc = (cur_col_stats.gap_count as f64 / cur_col_stats.column_len as f64) * 100.0;
            let num_residues = cur_col_stats.column_len - cur_col_stats.gap_count;
            let residue_perc = (num_residues as f64 / cur_col_stats.column_len as f64) * 100.0;

            let residue_counts_str = cur_col_stats
                .residue_counts
                .iter()
                .map(|&(res, count)| format!("{}:{}", res as char, count))
                .collect::<Vec<String>>()
                .join(",");

            let residue_perc_str = cur_col_stats
                .residue_counts
                .iter()
                .map(|&(res, count)| format!("{}:{:.2}", res as char, (count as f64 / cur_col_stats.column_len as f64) * 100.0))
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

        file.flush()?;

        Ok(())
    }
}
