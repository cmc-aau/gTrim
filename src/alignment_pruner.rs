use std::collections::{HashMap, HashSet};

use anyhow::Result;
use log::info;
use ndarray::Array2;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::iter::{ParallelIterator, IntoParallelIterator};

use crate::alignment_reader::AlignmentReader;
use crate::chi_square_calculator::ChiSquareCalculator;
use crate::permutation_test::PermutationTest;
use crate::sequence_type_detector::SequenceTypeDetector;
use crate::utils::{
    make_score_dict, write_alignment, write_pruned_dict_to_tsv, write_score_dict_to_json,
    write_score_dict_to_tsv, PruneInfo, ScoreDict,
};

pub struct AlignmentPruner {
    file: String,
    max_residue_pruned: usize,
    permutations: usize,
    residues_to_prune_perc: f64,
    initial_residues_to_prune: usize,
    residues_to_prune: usize,
    touchdown: bool,
    threads_permutations: usize,
    alignment_size: usize,
    progress: Option<ProgressBar>,
}

struct PerRowStats {
    count_rows: Array2<f64>,
    expected_observed: Array2<f64>,
    per_row_chi2: Vec<f64>,
    median: f64,
    mean: f64,
    q95: f64,
    sum: f64,
}

impl AlignmentPruner {
    pub fn new(
        file: String,
        max_residue_pruned: usize,
        permutations: usize,
        residues_to_prune_perc: f64,
        touchdown: bool,
        threads_permutations: usize,
    ) -> Self {
        let pb = ProgressBar::new_spinner();
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap()
                .tick_chars("⣾⣷⣯⣟⡿⢿⣻⣽"),
        );

        Self {
            file,
            max_residue_pruned,
            permutations,
            residues_to_prune_perc,
            initial_residues_to_prune: 0,
            residues_to_prune: 0,
            touchdown,
            threads_permutations,
            alignment_size: 0,
            progress: Some(pb),
        }
    }

    pub fn run(&mut self) -> Result<()> {
        info!("Starting pruning workflow.");
        
        info!("Reading alignment from file: {}", self.file);
        let reader = AlignmentReader::new(self.file.clone());
        let (alignment, alignment_array) = reader.run()?;
        info!("Alignment contains {} sequences and {} columns.", alignment.len(), alignment_array.ncols());

        let (is_dna, char_set) = SequenceTypeDetector::detect(&alignment)?;
        info!("Sequence type detected: {}", if is_dna { "DNA" } else { "amino acid" });

        // determine initial residues to prune per iteration
        self.initial_residues_to_prune = ((self.residues_to_prune_perc / 100.0) * alignment_array.ncols() as f64).floor() as usize;
        if self.initial_residues_to_prune == 0 {
            self.initial_residues_to_prune = 1;
        }
        self.residues_to_prune = self.initial_residues_to_prune;
        info!("Initially columns to prune per iteration: {}", self.initial_residues_to_prune);

        let permutation_test = PermutationTest::new(self.threads_permutations, self.permutations);
        let chi_square_calculator = ChiSquareCalculator::new(char_set);

        info!("Running {} permutations:", self.permutations);
        let (sums, _maxes, _upper_box_threshold, upper_threshold, permutated_per_row_chi2) =
            permutation_test.run(&alignment_array, &chi_square_calculator);

        info!("Starting recursive pruning:");
        let (pruned_alignment_array, prune_dict, score_dict) = self.recursive_prune(
            alignment_array,
            &chi_square_calculator,
            &permutation_test,
            &sums,
            upper_threshold,
            &permutated_per_row_chi2,
        )?;

        // Write outputs
        let suffix = if self.touchdown {
            format!("_p{}_s{}_m{}_touchdown_pruned.fasta", self.permutations, self.initial_residues_to_prune, self.max_residue_pruned)
        } else {
            format!("_p{}_s{}_m{}_pruned.fasta", self.permutations, self.initial_residues_to_prune, self.max_residue_pruned)
        };

        let base_path = self.file.trim_end_matches(".faa.gz").trim_end_matches(".faa")
            .trim_end_matches(".fna.gz").trim_end_matches(".fna")
            .trim_end_matches(".fasta").trim_end_matches(".fasta.gz");
        let output_alignment_file = format!("{}{}", base_path, suffix);
        let output_tsv_file = suffix.replace(".fasta", ".tsv");
        let output_tsv_file = format!("{}{}", base_path, output_tsv_file);

        write_alignment(&alignment, &pruned_alignment_array, &output_alignment_file)?;
        write_pruned_dict_to_tsv(&prune_dict, &output_tsv_file)?;

        // Calculate and write empirical p-values
        let empirical_pvalues = permutation_test.calc_empirical_pvalue(
            &score_dict.after_real,
            &permutated_per_row_chi2,
        );

        let row_empirical_pvalue_dict = make_score_dict(
            &score_dict.after_real,
            &permutated_per_row_chi2,
            &empirical_pvalues,
            &alignment,
        );

        let output_score_tsv_file = suffix.replace(".fasta", "_scores.tsv");
        let output_score_tsv_file = format!("{}{}", base_path, output_score_tsv_file);
        write_score_dict_to_tsv(&row_empirical_pvalue_dict, &output_score_tsv_file)?;

        let output_json_file = suffix.replace(".fasta", "_score_dict.json");
        let output_json_file = format!("{}{}", base_path, output_json_file);
        write_score_dict_to_json(&score_dict, &output_json_file)?;

        info!("Taxa p-values and z-scores printed to {}", output_score_tsv_file);
        info!("Chiscore values from permutation testing printed to {}", output_json_file);

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn recursive_prune(
        &mut self,
        mut alignment_array: Array2<char>,
        chi_square_calculator: &ChiSquareCalculator,
        permutation_test: &PermutationTest,
        sums: &[f64],
        upper_threshold: f64,
        permutated_per_row_chi2: &[f64],
    ) -> Result<(Array2<char>, HashMap<usize, PruneInfo>, ScoreDict)> {
        let mut prune_dict = HashMap::new();
        let mut score_dict = ScoreDict {
            before_permuted: permutated_per_row_chi2.to_vec(),
            before_real: Vec::new(),
            after_real: Vec::new(),
        };

        let mut iteration = 0;
        let mut removed_columns_count = 0;
        self.alignment_size = alignment_array.ncols();
        let mut original_indices: Vec<usize> = (0..alignment_array.ncols()).collect();

        let mean_perm_chi2: f64 =
            permutated_per_row_chi2.iter().sum::<f64>() / permutated_per_row_chi2.len() as f64;
        let variance: f64 = permutated_per_row_chi2
            .iter()
            .map(|&x| (x - mean_perm_chi2).powi(2))
            .sum::<f64>()
            / permutated_per_row_chi2.len() as f64;
        let sd_perm_chi2 = variance.sqrt();

        let mut residues_to_prune_indices: Vec<(usize, f64)> = Vec::with_capacity(self.residues_to_prune);
        let mut initial_global_chi2 = 0.0;

        loop {
            if removed_columns_count >= self.max_residue_pruned {
                break;
            }

            let stats = self.calculate_per_row_stats(&alignment_array, chi_square_calculator);

            let (alignment_empirical_p, significant_count) = self.calc_empirical_pvals(
                &stats,
                sums,
                permutated_per_row_chi2,
                permutation_test,
            );

            if iteration == 0 {
                score_dict.before_real = stats.per_row_chi2.clone();
            } else {
                // Record pruned columns
                for &(col, chi2_diff) in &residues_to_prune_indices {
                    let pruned_col = original_indices[col];
                    prune_dict.insert(
                        pruned_col,
                        PruneInfo {
                            iteration,
                            original_index: pruned_col,
                            global_chi_score: stats.sum,
                            initial_significant_chi_score: initial_global_chi2,
                            chi_score_difference: chi2_diff,
                            significant_taxa: significant_count,
                        },
                    );
                    removed_columns_count += 1;
                }

                // Update original indices - remove pruned columns
                let cols_to_remove: HashSet<usize> = residues_to_prune_indices.iter().map(|(c, _)| *c).collect();
                let mut new_indices = Vec::with_capacity(original_indices.len() - cols_to_remove.len());
                for (j, &idx) in original_indices.iter().enumerate() {
                    if !cols_to_remove.contains(&j) {
                        new_indices.push(idx);
                    }
                }
                original_indices = new_indices;
            }

            let (should_stop, stop_reason) = self.check_stopping_criteria(
                &stats,
                permutated_per_row_chi2,
                alignment_empirical_p,
                significant_count,
            );

            if let Some(pb) = &self.progress {
                pb.set_message(format!(
                    "Iteration: {} | Columns removed: {} | Total columns removed: {} ({:.2}%) | Biased taxa: {} | Mean Z-score: {:.2} | q95 Z-score: {:.2} | Alignment p-value: {:.2}",
                    iteration,
                    self.residues_to_prune,
                    removed_columns_count,
                    (removed_columns_count as f64 / self.alignment_size as f64) * 100.0,
                    significant_count,
                    (stats.mean - mean_perm_chi2) / sd_perm_chi2,
                    (stats.q95 - upper_threshold) / sd_perm_chi2,
                    alignment_empirical_p
                ));
                pb.tick();
            }

            if should_stop {
                // Finish the progress bar
                if let Some(pb) = &self.progress {
                    pb.finish_with_message("Pruning complete.");
                }
                self.progress = None;

                info!("Pruning stopped because of {}.", stop_reason);
                break;
            }

            // Calculate chi2 differences for pruning
            let chi2_differences = self.prune_quartic(
                &stats.count_rows,
                &stats.expected_observed,
                &alignment_array,
                chi_square_calculator,
            );

            initial_global_chi2 = chi_square_calculator
                .calculate_quartic_row_global_chi2(&stats.expected_observed, &stats.count_rows);

            // Sort columns by chi2 difference and get top N
            // Use index as secondary sort key for deterministic ordering
            let mut diff_vec: Vec<(usize, f64)> = chi2_differences.into_iter().collect();
            diff_vec.sort_by(|a, b| {
                a.1.partial_cmp(&b.1)
                    .unwrap()
                    .then_with(|| a.0.cmp(&b.0))
            });

            let num_to_remove = self.residues_to_prune.min(diff_vec.len());
            residues_to_prune_indices = diff_vec[diff_vec.len() - num_to_remove..].to_vec();
            residues_to_prune_indices.sort_by(|a, b| a.0.cmp(&b.0));

            // Remove all columns in a single pass
            let cols_to_remove: Vec<usize> = residues_to_prune_indices.iter().map(|(idx, _)| *idx).collect();
            alignment_array = remove_columns(&alignment_array, &cols_to_remove);

            iteration += 1;
        }

        score_dict.after_real = self
            .calculate_per_row_stats(&alignment_array, chi_square_calculator)
            .per_row_chi2;

        self.residues_to_prune = self.initial_residues_to_prune;
        Ok((alignment_array, prune_dict, score_dict))
    }

    fn prune_quartic(
        &self,
        count_rows_array: &Array2<f64>,
        expected_observed: &Array2<f64>,
        alignment_array: &Array2<char>,
        chi_square_calculator: &ChiSquareCalculator,
    ) -> HashMap<usize, f64> {
        let initial_global_chi2 = chi_square_calculator
            .calculate_quartic_row_global_chi2(expected_observed, count_rows_array);

        chi_square_calculator.calculate_quartic_chi2_difference(
            count_rows_array,
            alignment_array,
            initial_global_chi2,
        )
    }

    fn check_stopping_criteria(
        &mut self,
        stats: &PerRowStats,
        permuted_chi2: &[f64],
        alignment_empirical_p: f64,
        significant_count: usize,
    ) -> (bool, String) {
        if self.touchdown {
            let mut sorted_perm = permuted_chi2.to_vec();
            sorted_perm.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let p99_idx = ((99.0 / 100.0) * (sorted_perm.len() - 1) as f64) as usize;
            let p99 = sorted_perm[p99_idx];

            if stats.median <= p99 {
                let new_residues_to_prune = self.alignment_size / 1000;
                if new_residues_to_prune < self.residues_to_prune {
                    self.residues_to_prune = new_residues_to_prune;
                }

                if self.residues_to_prune == 0 {
                    self.residues_to_prune = 1;
                }
            }
        }

        if alignment_empirical_p >= 0.95 {
            return (true, "alignment p-value".to_string());
        }

        if significant_count == 0 {
            return (true, "no significant taxa".to_string());
        }

        (false, String::new())
    }

    fn calc_empirical_pvals(
        &self,
        stats: &PerRowStats,
        permuted_sums: &[f64],
        permuted_chi2: &[f64],
        permutation_test: &PermutationTest,
    ) -> (f64, usize) {
        let alignment_empirical_p =
            permutation_test.calc_empirical_pvalue(&[stats.sum], permuted_sums)[0];

        let empirical_pvals =
            permutation_test.calc_empirical_pvalue(&stats.per_row_chi2, permuted_chi2);

        let significant_count = empirical_pvals.iter().filter(|&&p| p <= 0.05).count();

        (alignment_empirical_p, significant_count)
    }

    fn calculate_per_row_stats(
        &self,
        alignment_array: &Array2<char>,
        chi_square_calculator: &ChiSquareCalculator,
    ) -> PerRowStats {
        let count_rows_array = chi_square_calculator.calculate_row_counts_parallel(alignment_array);
        let expected_observed = chi_square_calculator.calculate_expected_observed(&count_rows_array);
        let per_row_chi2 = chi_square_calculator.calculate_row_chi2(&expected_observed, &count_rows_array);

        let mut sorted_chi2 = per_row_chi2.clone();
        sorted_chi2.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median_idx = (0.5 * sorted_chi2.len() as f64) as usize;
        let q95_idx = (0.95 * sorted_chi2.len() as f64) as usize;

        PerRowStats {
            count_rows: count_rows_array,
            expected_observed,
            per_row_chi2: per_row_chi2.clone(),
            median: sorted_chi2[median_idx],
            mean: per_row_chi2.iter().sum::<f64>() / per_row_chi2.len() as f64,
            q95: sorted_chi2[q95_idx],
            sum: per_row_chi2.iter().sum(),
        }
    }
}

/// Remove multiple columns from a 2D array efficiently in parallel
/// Uses bitset for O(1) lookup and parallelizes across rows
fn remove_columns(array: &Array2<char>, cols_to_remove: &[usize]) -> Array2<char> {
    let nrows = array.nrows();
    let ncols = array.ncols();
    let new_ncols = ncols - cols_to_remove.len();

    // Build a bitset: mark columns to keep
    let mut keep_cols = vec![true; ncols];
    for &col_idx in cols_to_remove {
        keep_cols[col_idx] = false;
    }

    // Parallelize row-wise filtering
    let rows: Vec<Vec<char>> = (0..nrows)
        .into_par_iter()
        .map(|row_idx| {
            let mut row_data = Vec::with_capacity(new_ncols);
            for (col_idx, &keep) in keep_cols.iter().enumerate() {
                if keep {
                    row_data.push(array[[row_idx, col_idx]]);
                }
            }
            row_data
        })
        .collect();

    // Flatten and construct result array
    let new_data: Vec<char> = rows.into_iter().flat_map(|row| row.into_iter()).collect();
    Array2::from_shape_vec((nrows, new_ncols), new_data).unwrap()
}
