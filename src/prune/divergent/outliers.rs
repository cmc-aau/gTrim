use std::cmp::Ordering;
use std::collections::BTreeMap;

use clap::ValueEnum;

/// Methods for calculating the expected value of each locus comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum OutlierMethod {
    /// Running median over the N points nearest in x, made non-decreasing in x by isotonic regression. The expected
    /// value should not fall as relatedness rises, so local dips in the running median are treated as noise.
    MonotoneMedian,
    /// Running median over the N points nearest in x.
    Median,
}

/// Median of a set of values (mean of the two middle values for an even count).
fn median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        0.5 * (values[n / 2 - 1] + values[n / 2])
    }
}

/// Running median of y where the window for each point is the `window` points nearest to it in x.
/// Points must be sorted by x. If there are fewer than `window` points, all points are used.
///
/// Points with the same x are equally near, so the window is extended to include every point tied
/// with the x values at its ends. Points with the same x therefore get the same expected value, and
/// the result does not depend on the order of tied points.
pub fn running_median(xs: &[f64], ys: &[f64], window: usize) -> Vec<f64> {
    let n = xs.len();
    if n == 0 {
        return Vec::new();
    }

    // the N nearest points to x[p] form a contiguous window [lo, lo + w) in sorted order,
    // and lo is non-decreasing as p increases
    let w = window.clamp(1, n);
    let mut lo = 0;
    let mut buf = Vec::with_capacity(w);
    let mut expected = Vec::with_capacity(n);
    let mut p = 0;
    while p < n {
        // points tied with x[p] share the same window
        let mut run_end = p + 1;
        while run_end < n && xs[run_end] == xs[p] {
            run_end += 1;
        }

        while lo + w < n && (xs[lo + w] - xs[p]) < (xs[p] - xs[lo]) {
            lo += 1;
        }

        let mut start = lo;
        while start > 0 && xs[start - 1] == xs[start] {
            start -= 1;
        }
        let mut end = lo + w;
        while end < n && xs[end] == xs[end - 1] {
            end += 1;
        }

        buf.clear();
        buf.extend_from_slice(&ys[start..end]);
        let m = median(&mut buf);
        expected.extend(std::iter::repeat_n(m, run_end - p));
        p = run_end;
    }

    expected
}

/// Least-squares non-decreasing fit to `values` (pool adjacent violators). Points must be sorted by x. Points with
/// the same x start in the same block, so they always get the same fitted value.
pub fn non_decreasing_fit(xs: &[f64], values: &[f64]) -> Vec<f64> {
    let n = values.len();
    // blocks of (sum, count); each block's mean is at least the previous block's mean
    let mut blocks: Vec<(f64, usize)> = Vec::new();
    let mut p = 0;
    while p < n {
        let mut end = p + 1;
        while end < n && xs[end] == xs[p] {
            end += 1;
        }
        let mut block = (values[p..end].iter().sum::<f64>(), end - p);
        while let Some(&(sum, count)) = blocks.last() {
            if sum / count as f64 > block.0 / block.1 as f64 {
                blocks.pop();
                block = (block.0 + sum, block.1 + count);
            } else {
                break;
            }
        }
        blocks.push(block);
        p = end;
    }

    blocks.iter().flat_map(|&(sum, count)| std::iter::repeat_n(sum / count as f64, count)).collect()
}

/// Expected y value for each point under the specified outlier method. Points must be sorted by x.
pub fn expected_values(method: OutlierMethod, xs: &[f64], ys: &[f64], window: usize) -> Vec<f64> {
    match method {
        OutlierMethod::MonotoneMedian => non_decreasing_fit(xs, &running_median(xs, ys, window)),
        OutlierMethod::Median => running_median(xs, ys, window),
    }
}

/// Outlier counts for a focal genome and locus.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct VoteCounts {
    pub comparisons: usize,
    pub outliers: usize,
}

impl VoteCounts {
    /// Locus is masked if the majority (>50%) of comparisons are outliers.
    pub fn is_masked(&self) -> bool {
        2 * self.outliers > self.comparisons
    }
}

/// Tally outlier votes for each (genome, locus) pair.
pub fn tally_votes(
    keys: impl Iterator<Item = (usize, usize)>,
    outlier_flags: impl Iterator<Item = bool>,
) -> BTreeMap<(usize, usize), VoteCounts> {
    let mut votes: BTreeMap<(usize, usize), VoteCounts> = BTreeMap::new();
    for (key, is_outlier) in keys.zip(outlier_flags) {
        let v = votes.entry(key).or_default();
        v.comparisons += 1;
        v.outliers += is_outlier as usize;
    }

    votes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median() {
        assert_eq!(median(&mut [3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&mut [4.0, 1.0, 3.0, 2.0]), 2.5);
    }

    #[test]
    fn test_running_median_window_nearest_in_x() {
        // the last point is far from the others, so it is excluded from the window of
        // point 3 (points 1, 2, 3) and its own window is points 2, 3, 4
        let xs = [0.0, 1.0, 2.0, 3.0, 100.0];
        let ys = [10.0, 20.0, 30.0, 40.0, 50.0];
        let expected = running_median(&xs, &ys, 3);
        assert_eq!(expected, vec![20.0, 20.0, 30.0, 30.0, 40.0]);
    }

    #[test]
    fn test_running_median_ties() {
        // a run of tied x longer than the window uses the whole run, so the
        // expected value doesn't depend on the order of the tied points
        let xs = [1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 9.0];
        let expected = running_median(&xs, &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 100.0], 3);
        assert_eq!(&expected[1..6], &[3.0; 5]);
        let reordered = running_median(&xs, &[0.0, 5.0, 4.0, 3.0, 2.0, 1.0, 100.0], 3);
        assert_eq!(expected, reordered);

        // a window ending inside a run of ties is extended to the whole run
        // (the window of x = 1 is x = 0, 1, 2, and all three points at x = 2 are included)
        let xs = [0.0, 1.0, 2.0, 2.0, 2.0];
        assert_eq!(running_median(&xs, &[10.0, 20.0, 30.0, 40.0, 50.0], 3)[1], 30.0);
    }

    #[test]
    fn test_running_median_small_input() {
        let xs = [0.0, 1.0];
        let ys = [10.0, 20.0];
        assert_eq!(running_median(&xs, &ys, 40), vec![15.0, 15.0]);
        assert!(running_median(&[], &[], 40).is_empty());
    }

    #[test]
    fn test_non_decreasing_fit() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0];
        // already non-decreasing: unchanged
        assert_eq!(non_decreasing_fit(&xs, &[1.0, 2.0, 2.0, 5.0, 9.0]), vec![1.0, 2.0, 2.0, 5.0, 9.0]);
        // a dip is pooled with the values before it
        assert_eq!(non_decreasing_fit(&xs, &[1.0, 6.0, 2.0, 7.0, 8.0]), vec![1.0, 4.0, 4.0, 7.0, 8.0]);
        // pooling cascades backwards
        assert_eq!(non_decreasing_fit(&xs, &[3.0, 6.0, 5.0, 1.0, 10.0]), vec![3.0, 4.0, 4.0, 4.0, 10.0]);
        // decreasing input becomes its mean
        assert_eq!(non_decreasing_fit(&[1.0, 2.0, 3.0], &[9.0, 6.0, 3.0]), vec![6.0, 6.0, 6.0]);
        assert!(non_decreasing_fit(&[], &[]).is_empty());
    }

    #[test]
    fn test_non_decreasing_fit_ties() {
        // points with the same x always get the same value, even if a block boundary would otherwise fall between them
        let fit = non_decreasing_fit(&[1.0, 2.0, 2.0, 3.0], &[5.0, 4.0, 8.0, 9.0]);
        assert_eq!(fit, vec![5.0, 6.0, 6.0, 9.0]);
        let fit = non_decreasing_fit(&[1.0, 2.0, 2.0, 3.0], &[7.0, 4.0, 6.0, 9.0]);
        assert!((fit[0] - 17.0 / 3.0).abs() < 1e-12 && fit[0] == fit[1] && fit[1] == fit[2] && fit[3] == 9.0);
    }

    #[test]
    fn test_expected_values_methods() {
        let xs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = [10.0, 30.0, 20.0, 40.0, 50.0];
        let median = expected_values(OutlierMethod::Median, &xs, &ys, 1);
        assert_eq!(median, ys.to_vec());
        let monotone = expected_values(OutlierMethod::MonotoneMedian, &xs, &ys, 1);
        assert_eq!(monotone, vec![10.0, 25.0, 25.0, 40.0, 50.0]);
        assert!(monotone.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn test_majority_vote() {
        assert!(VoteCounts { comparisons: 3, outliers: 2 }.is_masked());
        assert!(!VoteCounts { comparisons: 3, outliers: 1 }.is_masked());
        assert!(!VoteCounts { comparisons: 2, outliers: 1 }.is_masked());
        assert!(VoteCounts { comparisons: 1, outliers: 1 }.is_masked());
    }

    #[test]
    fn test_tally_votes() {
        let keys = [(0, 0), (0, 0), (0, 0), (1, 0)];
        let flags = [true, true, false, false];
        let votes = tally_votes(keys.into_iter(), flags.into_iter());
        assert_eq!(votes[&(0, 0)], VoteCounts { comparisons: 3, outliers: 2 });
        assert_eq!(votes[&(1, 0)], VoteCounts { comparisons: 1, outliers: 0 });
    }
}
