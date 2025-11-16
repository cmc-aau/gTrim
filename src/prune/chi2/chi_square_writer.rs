use anyhow::{Context, Result};
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::alignment_writer;

/// Write alignment to FASTA file
pub fn write_alignment(
    seq_ids: &[String],
    alignment_array: &Array2<char>,
    output_file: &str,
    compress: bool,
) -> Result<()> {
    let mask = vec![true; alignment_array.ncols()];
    alignment_writer::write_filtered_alignment(seq_ids, alignment_array, &mask, output_file, compress)
}

/// Structure to hold pruning information for a column
#[derive(Debug)]
pub struct PruneInfo {
    pub iteration: usize,
    pub original_index: usize,
    pub global_chi_score: f64,
    pub initial_significant_chi_score: f64,
    pub chi_score_difference: f64,
    pub significant_taxa: usize,
}

/// Write pruning dictionary to TSV file
pub fn write_pruned_dict_to_tsv(
    prune_dict: &HashMap<usize, PruneInfo>,
    file_name: &str,
) -> Result<()> {
    let mut file = File::create(file_name)
        .with_context(|| format!("Failed to create file: {}", file_name))?;

    // Write header
    writeln!(
        file,
        "Iteration\tOriginal Index\tGlobal ChiScore\tInitial Significant ChiScore\tChiScore Difference\tSignificant taxa"
    )?;

    // Sort by iteration first, then by chi_score_difference (ascending - most negative first)
    let mut items: Vec<_> = prune_dict.iter().collect();
    items.sort_by(|(_, a), (_, b)| {
        a.iteration.cmp(&b.iteration)
            .then_with(|| a.chi_score_difference.partial_cmp(&b.chi_score_difference).unwrap())
    });

    for (_, info) in items {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}",
            info.iteration,
            info.original_index,
            info.global_chi_score,
            info.initial_significant_chi_score,
            info.chi_score_difference,
            info.significant_taxa
        )?;
    }

    Ok(())
}

/// Write binary mask where 1 indicate a column is retained and 0 indicates it was pruned.
pub fn write_mask(
    prune_dict: &HashMap<usize, PruneInfo>,
    num_cols: usize,
    out_file: &str,
) -> Result<()> {
    let mut mask = vec![true; num_cols];

    for info in prune_dict.values() {
        if info.original_index < num_cols {
            mask[info.original_index] = false;
        }
    }

    alignment_writer::write_mask(&mask, out_file)
}

/// Score dictionary structure for JSON output
#[derive(Debug, Serialize, Deserialize)]
pub struct ScoreDict {
    pub before_permuted: Vec<f64>,
    pub before_real: Vec<f64>,
    pub after_real: Vec<f64>,
}

/// Write score dictionary to JSON file
pub fn write_score_dict_to_json(score_dict: &ScoreDict, file_name: &str) -> Result<()> {
    let json = serde_json::to_string_pretty(score_dict)
        .context("Failed to serialize score dictionary to JSON")?;

    let mut file = File::create(file_name)
        .with_context(|| format!("Failed to create file: {}", file_name))?;

    file.write_all(json.as_bytes())?;
    Ok(())
}

/// Structure for per-row score information
#[derive(Debug)]
pub struct RowScoreInfo {
    pub empirical_pvalue: f64,
    pub zscore: f64,
}

/// Write score dictionary to TSV file, ordered by descending absolute z-score
pub fn write_score_dict_to_tsv(
    row_dict: &HashMap<String, RowScoreInfo>,
    file_name: &str,
) -> Result<()> {
    let mut file = File::create(file_name)
        .with_context(|| format!("Failed to create file: {}", file_name))?;

    // Write header
    writeln!(file, "Row\tEmpirical-Pvalue\tZ-Score")?;

    // Sort by absolute z-score in descending order
    let mut items: Vec<_> = row_dict.iter().collect();
    items.sort_by(|a, b| {
        b.1.zscore
            .abs()
            .partial_cmp(&a.1.zscore.abs())
            .unwrap()
    });

    for (row_name, info) in items {
        writeln!(
            file,
            "{}\t{}\t{}",
            row_name, info.empirical_pvalue, info.zscore
        )?;
    }

    Ok(())
}

/// Make a dictionary of chi-squared scores for each row in the alignment
pub fn make_score_dict(
    per_row_chi2: &[f64],
    permutated_per_row_chi2: &[f64],
    empirical_pvalues: &[f64],
    seq_ids: &[String],
) -> HashMap<String, RowScoreInfo> {
    let row_names: Vec<String> = seq_ids.to_vec();
    let mean_perm_chi2: f64 = permutated_per_row_chi2.iter().sum::<f64>() / permutated_per_row_chi2.len() as f64;

    let variance: f64 = permutated_per_row_chi2
        .iter()
        .map(|&x| (x - mean_perm_chi2).powi(2))
        .sum::<f64>() / permutated_per_row_chi2.len() as f64;
    let sd_perm_chi2 = variance.sqrt();

    let mut row_dict = HashMap::new();

    for i in 0..row_names.len() {
        let zscore = (per_row_chi2[i] - mean_perm_chi2) / sd_perm_chi2;
        row_dict.insert(
            row_names[i].clone(),
            RowScoreInfo {
                empirical_pvalue: empirical_pvalues[i],
                zscore,
            },
        );
    }

    row_dict
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_score_dict() {
        let per_row_chi2 = vec![10.0, 15.0];
        let permutated = vec![8.0, 9.0, 10.0, 11.0];
        let empirical_pvalues = vec![0.5, 0.25];
        let seq_ids = vec![
            "seq1".to_string(),
            "seq2".to_string(),
        ];

        let score_dict = make_score_dict(&per_row_chi2, &permutated, &empirical_pvalues, &seq_ids);

        assert_eq!(score_dict.len(), 2);
        assert!(score_dict.contains_key("seq1"));
        assert!(score_dict.contains_key("seq2"));
    }
}
