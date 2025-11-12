use std::collections::HashMap;

use ndarray::{Array1, Array2, Axis};
use rayon::prelude::*;


pub struct ChiSquareCalculator {
    char_set: Vec<char>,
    char_to_idx: HashMap<char, usize>,
}

impl ChiSquareCalculator {
    pub fn new(char_set: Vec<char>) -> Self {
        // Build character to index mapping for fast lookups
        let char_to_idx: HashMap<char, usize> = char_set
            .iter()
            .enumerate()
            .map(|(idx, &ch)| (ch, idx))
            .collect();

        Self {
            char_set,
            char_to_idx,
        }
    }

    /// Count occurrences of all characters in the char_set for each sequence (row)
    /// Returns a 2D array where rows are characters and columns are sequences
    pub fn calculate_row_counts(&self, alignment_array: &Array2<char>) -> Array2<f64> {
        let num_chars = self.char_set.len();
        let num_seqs = alignment_array.nrows();
        
        let mut counts_vec: Vec<Vec<f64>> = Vec::with_capacity(num_seqs);
        for seq_idx in 0..num_seqs {
            let row = alignment_array.row(seq_idx);
            let mut seq_counts = vec![0.0; num_chars];
            for (char_idx, &ch) in self.char_set.iter().enumerate() {
                seq_counts[char_idx] = row.iter().filter(|&&c| c == ch).count() as f64;
            }
            counts_vec.push(seq_counts);
        }

        // Build the count rows array from the results
        let mut count_rows = Array2::zeros((num_chars, num_seqs));
        for (seq_idx, seq_counts) in counts_vec.into_iter().enumerate() {
            for (char_idx, count) in seq_counts.into_iter().enumerate() {
                count_rows[[char_idx, seq_idx]] = count;
            }
        }

        // Avoid division by zero - if any zero exists, add 1 to all values
        if count_rows.iter().any(|&x| x == 0.0) {
            count_rows.mapv_inplace(|x| x + 1.0);
        }

        count_rows
    }

    /// Count occurrences of all characters in the char_set for each sequence (row) in parallel
    /// Returns a 2D array where rows are characters and columns are sequences
    pub fn calculate_row_counts_parallel(&self, alignment_array: &Array2<char>) -> Array2<f64> {
        let num_chars = self.char_set.len();
        let num_seqs = alignment_array.nrows();

        // Compute counts in parallel using global threadpool - each sequence processed independently
        let counts_vec: Vec<Vec<f64>> = (0..num_seqs)
            .into_par_iter()
            .map(|seq_idx| {
                let row = alignment_array.row(seq_idx);
                let mut seq_counts = vec![0.0; num_chars];
                for (char_idx, &ch) in self.char_set.iter().enumerate() {
                    seq_counts[char_idx] = row.iter().filter(|&&c| c == ch).count() as f64;
                }
                seq_counts
            })
            .collect();

        // Build the count rows array from the parallel results
        let mut count_rows = Array2::zeros((num_chars, num_seqs));
        for (seq_idx, seq_counts) in counts_vec.into_iter().enumerate() {
            for (char_idx, count) in seq_counts.into_iter().enumerate() {
                count_rows[[char_idx, seq_idx]] = count;
            }
        }

        // Avoid division by zero - if any zero exists, add 1 to all values
        if count_rows.iter().any(|&x| x == 0.0) {
            count_rows.mapv_inplace(|x| x + 1.0);
        }

        count_rows
    }

    /// Fully vectorized calculation of expected and observed frequencies for each sequence
    pub fn calculate_expected_observed(&self, count_rows_array: &Array2<f64>) -> Array2<f64> {
        // Sum across characters for each sequence (axis 0)
        let total_chars: Array1<f64> = count_rows_array.sum_axis(Axis(0));
        let total_chars_sum: f64 = total_chars.sum();

        // Calculate global frequencies for each character
        let row_sums: Array1<f64> = count_rows_array.sum_axis(Axis(1));
        let global_frequencies = &row_sums / total_chars_sum;

        // Calculate expected frequencies: outer product of global_frequencies and sequence totals
        let num_chars = count_rows_array.nrows();
        let num_seqs = count_rows_array.ncols();
        let mut expected_frequencies = Array2::zeros((num_chars, num_seqs));

        for i in 0..num_chars {
            for j in 0..num_seqs {
                expected_frequencies[[i, j]] = global_frequencies[i] * total_chars[j];
            }
        }

        // Avoid division by zero - if any zero exists, add 1 to all values
        if expected_frequencies.iter().any(|&x| x == 0.0) {
            expected_frequencies.mapv_inplace(|x| x + 1.0);
        }

        expected_frequencies
    }

    /// Calculate the chi-squared score for all taxa (sequences)
    pub fn calculate_row_chi2(&self, expected_values: &Array2<f64>, count_rows_array: &Array2<f64>) -> Vec<f64> {
        let diff = count_rows_array - expected_values;
        let chi2_values = &diff * &diff / expected_values;

        // Sum across characters (axis 0) to get chi2 for each sequence
        chi2_values.sum_axis(Axis(0)).to_vec()
    }

    /// Calculate the global quartic chi-squared score (sum of squared chi-squares)
    pub fn calculate_quartic_row_global_chi2(&self, expected_values: &Array2<f64>, count_rows_array: &Array2<f64>) -> f64 {
        let chi2_sum = self.calculate_row_chi2(expected_values, count_rows_array);
        chi2_sum.iter().map(|&x| x * x).sum()
    }

    /// Calculation of quartic chi2 for a single column
    /// Optimized for large arrays (200k+ sequences) with cache-friendly access
    fn calculate_chi2_for_column(
        &self,
        count_rows_array: &Array2<f64>,
        alignment_array: &Array2<char>,
        col_idx: usize,
    ) -> f64 {
        let num_chars = self.char_set.len();
        let num_seqs = count_rows_array.ncols();

        // Extract column view once - critical for large arrays (200k rows × 40k cols)
        // This avoids strided access of alignment_array[[seq_idx, col_idx]]
        let column_view = alignment_array.column(col_idx);

        // Pre-compute character indices using column view
        let mut seq_char_indices: Vec<Option<usize>> = Vec::with_capacity(num_seqs);
        for &char_at_pos in column_view.iter() {
            seq_char_indices.push(self.char_to_idx.get(&char_at_pos).copied());
        }

        // Pre-allocate all buffers
        let mut seq_totals = vec![0.0; num_seqs];
        let mut char_totals = vec![0.0; num_chars];
        let mut adjusted_values = Vec::with_capacity(num_chars * num_seqs);

        // For count_rows_array, direct indexing is fine since it's (num_chars, num_seqs)
        // Access pattern [[char_idx, seq_idx]] for fixed char_idx is sequential
        for char_idx in 0..num_chars {
            for seq_idx in 0..num_seqs {
                let original_count = count_rows_array[[char_idx, seq_idx]];
                let matches = seq_char_indices[seq_idx] == Some(char_idx);
                let adjusted = original_count - 1.0 - if matches { 1.0 } else { 0.0 };

                adjusted_values.push(adjusted);
                seq_totals[seq_idx] += adjusted;
                char_totals[char_idx] += adjusted;
            }
        }

        // Calculate total and global frequencies
        let total_chars: f64 = char_totals.iter().sum();
        let global_frequencies: Vec<f64> = char_totals
            .iter()
            .map(|&x| x / total_chars)
            .collect();

        // Calculate quartic chi2
        let mut quartic_chi2 = 0.0;

        for seq_idx in 0..num_seqs {
            let mut seq_chi2 = 0.0;
            let seq_total = seq_totals[seq_idx];

            for char_idx in 0..num_chars {
                let observed = adjusted_values[char_idx * num_seqs + seq_idx];
                let mut expected = global_frequencies[char_idx] * seq_total;

                // Handle zero division case
                if expected == 0.0 {
                    expected = 1.0;
                }

                let diff = observed - expected;
                seq_chi2 += (diff * diff) / expected;
            }

            quartic_chi2 += seq_chi2 * seq_chi2;
        }

        quartic_chi2
    }

    /// Calculate the chi-squared difference for each column
    /// based on its impact on the squared chi-squared score of rows.
    pub fn calculate_quartic_chi2_difference(
        &self,
        count_rows_array: &Array2<f64>,
        alignment_array: &Array2<char>,
        initial_global_chi2: f64,
    ) -> HashMap<usize, f64> {
        let num_cols = alignment_array.ncols();

        // Use global thread pool
        let results: Vec<(usize, f64)> = (0..num_cols)
            .into_par_iter()
            .map(|col_idx| {
                let chi2_without_col = self.calculate_chi2_for_column(
                    count_rows_array,
                    alignment_array,
                    col_idx
                );
                let difference = initial_global_chi2 - chi2_without_col;
                (col_idx, difference)
            })
            .collect();

        results.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_row_counts() {
        let char_set = vec!['A', 'C', 'G', 'T'];
        let calculator = ChiSquareCalculator::new(char_set);

        // Create a simple 5x4 array (5 sequences, 4 columns)
        let alignment = Array2::from_shape_vec(
            (5, 4),
            vec![
                'A', 'C', 'G', 'T',
                'A', 'C', 'G', 'T',
                'A', 'C', 'G', 'T',
                'A', 'C', 'G', 'T',
                'A', 'A', 'A', 'A',
            ],
        ).unwrap();

        let counts = calculator.calculate_row_counts(&alignment);

        // After zero-avoidance (+1 to all values because sequence 4 has zeros)
        // Sequence 0 has 1 A, 1 C, 1 G, 1 T → becomes 2, 2, 2, 2
        assert_eq!(counts[[0, 0]], 2.0); // 'A' in sequence 0
        assert_eq!(counts[[1, 0]], 2.0); // 'C' in sequence 0

        // Sequence 4 has 4 As, 0 Cs, 0 Gs, 0 Ts → becomes 5, 1, 1, 1
        assert_eq!(counts[[0, 4]], 5.0); // 'A' in sequence 4
        assert_eq!(counts[[1, 4]], 1.0); // 'C' in sequence 4

        // perform same test with parallel implementation
        let counts = calculator.calculate_row_counts_parallel(&alignment);

        // After zero-avoidance (+1 to all values because sequence 4 has zeros)
        // Sequence 0 has 1 A, 1 C, 1 G, 1 T → becomes 2, 2, 2, 2
        assert_eq!(counts[[0, 0]], 2.0); // 'A' in sequence 0
        assert_eq!(counts[[1, 0]], 2.0); // 'C' in sequence 0

        // Sequence 4 has 4 As, 0 Cs, 0 Gs, 0 Ts → becomes 5, 1, 1, 1
        assert_eq!(counts[[0, 4]], 5.0); // 'A' in sequence 4
        assert_eq!(counts[[1, 4]], 1.0); // 'C' in sequence 4
    }
}
