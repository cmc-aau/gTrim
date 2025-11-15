use std::sync::Arc;

use ndarray::Array2;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::progress::progress_bar;
use super::chi_square_calculator::ChiSquareCalculator;

pub struct PermutationTest {
    pool: rayon::ThreadPool,
    permutations: usize,
}

impl PermutationTest {
    pub fn new(threads: usize, permutations: usize) -> Self {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to create thread pool");

        Self {
            pool,
            permutations,
        }
    }

    /// Calculate empirical p-values
    /// If per_row_chi2 is a vector (per taxon), calculate p-value for each taxon
    /// If per_row_chi2 is a single value (alignment sum), calculate single p-value
    pub fn calc_empirical_pvalue(
        &self,
        per_row_chi2: &[f64],
        permutated_per_row_chi2: &[f64],
    ) -> Vec<f64> {
        let mut sorted_perms: Vec<f64> = permutated_per_row_chi2.to_vec();
        sorted_perms.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n_perms = sorted_perms.len();
        let p_len = per_row_chi2.len();

        self.pool.install(|| {
            per_row_chi2
                .par_iter()
                .map(|&chi2_value| {
                    // find first index >= chi2_value via binary search style:
                    // partition_point is ideal: first index where x >= chi2_value
                    let idx = sorted_perms.partition_point(|&x| x < chi2_value);
                    let count = n_perms - idx;

                    if p_len > 1 {
                        (count as f64 / n_perms as f64) * p_len as f64
                    } else {
                        count as f64 / n_perms as f64
                    }
                })
                .collect::<Vec<f64>>()
        })
    }

    /// Permute alignment columns and calculate chi-squared values
    fn permute_and_calculate_chi2(
        &self,
        alignment_array: &Array2<char>,
        chi_square_calculator: &ChiSquareCalculator,
    ) -> Vec<Vec<f64>> {
        let num_rows = alignment_array.nrows();
        let num_cols = alignment_array.ncols();

        // Progress bar for permutations
        let pb = progress_bar(self.permutations as u64);
        let pb = Arc::new(pb);

        // Perform permutations in parallel
        let results: Vec<Vec<f64>> = self.pool.install(|| {
            (0..self.permutations)
                .into_par_iter()
                .map(|i| {
                    // Use fixed seed for reproducibility
                    let seed = 12345u64 + i as u64;
                    let mut rng = ChaCha8Rng::seed_from_u64(seed);

                    // Reuse a single permuted array buffer (allocated once per thread)
                    let mut permuted_array = alignment_array.clone();

                    // Reuse a single column buffer
                    let mut column = Vec::with_capacity(num_rows);

                    for col_idx in 0..num_cols {
                        column.clear();
                        column.extend(permuted_array.column(col_idx).iter().copied());
                        column.shuffle(&mut rng);
                        for (row_idx, &value) in column.iter().enumerate() {
                            permuted_array[[row_idx, col_idx]] = value;
                        }
                    }

                    // Calculate chi2 for permuted array
                    let count_rows_array = chi_square_calculator.calculate_row_counts(&permuted_array);
                    let expected_observed = chi_square_calculator.calculate_expected_observed(&count_rows_array);
                    let per_row_chi2 = chi_square_calculator.calculate_row_chi2(&expected_observed, &count_rows_array);

                    // Increment progress bar
                    let pb = pb.clone();
                    pb.inc(1);

                    per_row_chi2
                })
                .collect()
        });

        // Finish progress bar
        Arc::try_unwrap(pb)
            .map(|pb| pb.finish())
            .ok();

        results
    }

    /// Run the permutation test and get chi-squared score percentiles and distribution
    /// Returns: (sums, maxes, upper_box_threshold, upper_threshold, permutated_per_row_chi2)
    pub fn run(
        &self,
        alignment_array: &Array2<char>,
        chi_square_calculator: &ChiSquareCalculator,
    ) -> (Vec<f64>, Vec<f64>, f64, f64, Vec<f64>) {
        let permutated_per_row_chi2 = self.permute_and_calculate_chi2(alignment_array, chi_square_calculator);

        // Calculate maxes, sums, and flatten in a single pass
        let (mut maxes, mut sums, mut flattened) = (Vec::new(), Vec::new(), Vec::new());
        for row in permutated_per_row_chi2 {
            let max = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let sum: f64 = row.iter().sum();
            maxes.push(max);
            sums.push(sum);
            flattened.extend(row);
        }

        // Calculate thresholds
        let upper_box_threshold = percentile(&flattened, 75.0);
        let upper_threshold = percentile(&flattened, 95.0);

        (sums, maxes, upper_box_threshold, upper_threshold, flattened)
    }
}

/// Calculate percentile of a vector using linear interpolation (matching NumPy's default)
fn percentile(data: &[f64], percentile: f64) -> f64 {
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let idx = percentile / 100.0 * (sorted.len() - 1) as f64;
    let lower = idx.floor() as usize;
    let upper = idx.ceil() as usize;
    let fraction = idx - lower as f64;

    if lower == upper {
        sorted[lower]
    } else {
        sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_percentile() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let p50 = percentile(&data, 50.0);
        assert!((p50 - 5.5).abs() < 0.1);
    }

    #[test]
    fn test_empirical_pvalue() {
        let pt = PermutationTest::new(1, 100);
        let per_row_chi2 = vec![10.0];
        let permutated = vec![5.0, 8.0, 12.0, 15.0]; // 2 out of 4 are >= 10.0
        let pvals = pt.calc_empirical_pvalue(&per_row_chi2, &permutated);
        assert_eq!(pvals.len(), 1);
        assert!((pvals[0] - 0.5).abs() < 0.01); // 2/4 = 0.5
    }
}
