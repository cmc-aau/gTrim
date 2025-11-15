use std::fs::File;
use std::io::Write;
use std::collections::HashSet;

use anyhow::{Context, Result};
use log::info;
use ndarray::Array2;
use rayon::iter::{ParallelIterator, IntoParallelIterator};

use crate::sequence_type_detector::{GAP_CHARS, SequenceTypeDetector};
use crate::alignment_reader::{AlignmentReader, SeqRecord};

pub struct GapPruner {
    file: String,
    gap_threshold: f64,
    max_residue_pruned: usize,
    out_dir: String,
}

#[derive(Default, Clone)]
struct PerColStats {
    num_gaps: usize,
    gap_perc: f64,
}

impl GapPruner {
    pub fn new(
        file: String,
        gap_threshold: f64,
        max_residue_pruned: usize,
        out_dir: String,
    ) -> Self {
        Self {
            file,
            gap_threshold,
            max_residue_pruned,
            out_dir,
        }
    }

    pub fn run(&mut self) -> Result<()> {
        info!("Starting gap pruning workflow.");

        info!("Reading alignment from file: {}", self.file);
        let reader = AlignmentReader::new(self.file.clone());
        let (alignment, alignment_array) = reader.run()?;
        info!("Alignment contains {} sequences and {} columns.", alignment.len(), alignment_array.ncols());

        let (is_dna, _char_set) = SequenceTypeDetector::detect(&alignment)?;
        info!("Sequence type detected: {}", if is_dna { "DNA" } else { "amino acid" });

        // determine columns exceeding gap threshold that should be filtered
        info!("Pruning gappy columns.");
        let col_stats = self.calculate_gap_statistics(&alignment_array);
        let mask = self.filter_gappy_columns(&col_stats, self.gap_threshold);
        info!("Pruned {} columns.", mask.iter().filter(|&state| !state).count());

        // Write outputs
        let suffix = format!("_gt{}_m{}_pruned.fasta", self.gap_threshold, self.max_residue_pruned);

        // get the filename from the input path
        let mut filename = std::path::Path::new(&self.file)
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| anyhow::anyhow!("Input path '{}' does not contain a valid filename", self.file))?;

        // create owned stem and strip known extensions
        filename = filename
            .trim_end_matches(".gz")
            .trim_end_matches(".faa")
            .trim_end_matches(".fna")
            .trim_end_matches(".fasta");

        // join output directory with the sanitized stem and return as String
        let base_path = std::path::Path::new(&self.out_dir)
            .join(filename)
            .to_string_lossy()
            .into_owned();

        let output_stats_file = suffix.replace(".fasta", ".tsv");
        let output_stats_file = format!("{}{}", base_path, output_stats_file);
        self.write_gap_stats(&col_stats, &mask, &output_stats_file)?;

        let output_mask_file = suffix.replace(".fasta", ".mask");
        let output_mask_file = format!("{}{}", base_path, output_mask_file);
        self.write_mask(&mask, &output_mask_file)?;

        info!("Writing filtered alignment to file.");
        let output_alignment_file = format!("{}{}", base_path, suffix);
        self.write_filtered_alignment(&alignment, &alignment_array, &mask, &output_alignment_file)?;

        Ok(())
    }

    /// Determine gap statistics for each column in alignment.
    fn calculate_gap_statistics(&self, alignment_array: &Array2<char>) -> Vec<PerColStats> {
        let gaps: HashSet<char> = HashSet::from_iter(GAP_CHARS.iter().copied());
        let num_rows = alignment_array.nrows() as f64;
        let mut gap_stats = Vec::with_capacity(alignment_array.ncols());
        for col_idx in 0..alignment_array.ncols() {
            let column = alignment_array.column(col_idx);
            let num_gaps = column.iter().filter(|&&c| gaps.contains(&c)).count();
            let gap_perc = (num_gaps as f64 / num_rows) * 100.0;
            gap_stats.push(PerColStats { num_gaps, gap_perc });
        }
        
        gap_stats
    }

    /// Filter columns that exceed the gap threshold.
    fn filter_gappy_columns(&self, gap_stats: &[PerColStats], gap_threshold: f64) -> Vec<bool> {
        // Build boolean mask: true = keep column, false = remove
        gap_stats
            .iter()
            .map(|s| s.gap_perc <= gap_threshold)
            .collect()
    }

    /// Write filtered alignment to FASTA file.
    pub fn write_filtered_alignment(
        &self, 
        records: &[SeqRecord],
        alignment_array: &Array2<char>,
        mask: &[bool],
        out_file: &str,
    ) -> Result<()> {
        let mut file = File::create(out_file)
            .with_context(|| format!("Failed to create file: {}", out_file))?;

        for (row_idx, record) in records.iter().enumerate() {
            writeln!(file, ">{}", record.id)?;
            let seq: String = alignment_array
                .row(row_idx)
                .iter()
                .copied()
                .enumerate()
                .filter_map(|(i, ch)| if mask[i] { Some(ch) } else { None })
                .collect();
            writeln!(file, "{}", seq)?;
        }

        Ok(())
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

    /// Write out mask (true = 1; false = 0).
    fn write_mask(&self, mask: &[bool], out_file: &str) -> Result<()> {
        // Convert bool mask to a string of '1' and '0' characters
        let mask_string: String = mask.iter().map(|&b| if b { '1' } else { '0' }).collect();

        let mut file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file))?;
        writeln!(file, "{}", mask_string)?;

        Ok(())
    }
}