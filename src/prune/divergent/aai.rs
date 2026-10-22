use std::cmp::Ordering;

use rayon::prelude::*;

use crate::progress::progress_bar;
use super::locus_set::{LocusSet, GAP_CODE};

/// Comparison of a locus between a focal genome and one of its closest neighbours.
#[derive(Debug, Clone)]
pub struct LocusComparison {
    pub locus_idx: usize,
    pub genome_idx: usize,
    pub neighbour_idx: usize,
    /// Rank of the neighbour among neighbours with this locus (1 = closest).
    pub neighbour_rank: usize,
    /// Genome-wide AAI (%) excluding the focal locus.
    pub similarity: f64,
    /// Percent identity of the focal locus.
    pub locus_identity: f64,
    /// Number of aligned positions for the focal locus.
    pub locus_aligned: u32,
}

/// Genome-wide AAI (%) between two genomes, pooled over all loci. None if they have no aligned positions.
pub fn genome_aai(locus_set: &LocusSet, genome_a: usize, genome_b: usize) -> Option<f64> {
    let (identical, aligned) = identity_counts(locus_set.genome_slice(genome_a), locus_set.genome_slice(genome_b));
    (aligned > 0).then(|| 100.0 * identical as f64 / aligned as f64)
}

/// Count identical and aligned positions between two encoded sequences.
/// A position is aligned if both sequences have a valid residue.
#[inline]
pub fn identity_counts(a: &[u8], b: &[u8]) -> (u32, u32) {
    let mut identical = 0u32;
    let mut aligned = 0u32;
    for (&x, &y) in a.iter().zip(b) {
        let both = (x != GAP_CODE) & (y != GAP_CODE);
        aligned += both as u32;
        identical += (both & (x == y)) as u32;
    }

    (identical, aligned)
}

/// Check if the number of aligned positions covers the required fraction of the shorter sequence.
#[inline]
pub fn sufficient_overlap(aligned: u32, nongap_a: u32, nongap_b: u32, min_aligned_frac: f64) -> bool {
    aligned > 0 && aligned as f64 >= min_aligned_frac * nongap_a.min(nongap_b) as f64
}

/// Minimum overlap for a neighbour to be compared to a focal genome at a locus.
#[derive(Debug, Clone, Copy)]
pub struct OverlapFilter {
    /// Minimum fraction of the shorter sequence that must be aligned, at the locus and across all other loci.
    pub min_aligned_frac: f64,
    /// Minimum number of aligned positions at the locus. Identity over a few positions is too noisy to test.
    pub min_aligned_columns: u32,
}

/// Identify closest neighbours and compare each locus of a genome to its closest neighbours.
///
/// Two sets of neighbours are identified for each locus of a genome:
///   - lower: closest genomes by AAI calculated with the locus left out, used to identify loci
///     with unexpectedly low identity
///   - upper: closest genomes by identity of the locus, used to identify loci with unexpectedly
///     high identity relative to genome-wide AAI (e.g. contamination or lateral transfer), and to
///     report the genome whose copy of the locus is most similar
pub struct NeighbourComparer<'a> {
    locus_set: &'a LocusSet,
    neighbours: usize,
    upper_neighbours: usize,
    overlap: OverlapFilter,
}

impl<'a> NeighbourComparer<'a> {
    /// Create comparer that selects `neighbours` by AAI and `upper_neighbours` by locus identity.
    pub fn new(locus_set: &'a LocusSet, neighbours: usize, upper_neighbours: usize, overlap: OverlapFilter) -> Self {
        Self {
            locus_set,
            neighbours,
            upper_neighbours,
            overlap,
        }
    }

    /// Compare loci of all genomes to their closest neighbours. Returns comparisons to neighbours
    /// selected by AAI (lower) and by locus identity (upper).
    pub fn run(&self) -> (Vec<LocusComparison>, Vec<LocusComparison>) {
        let num_genomes = self.locus_set.num_genomes();
        let pb = progress_bar(num_genomes as u64);

        let comparisons: Vec<(Vec<LocusComparison>, Vec<LocusComparison>)> = (0..num_genomes)
            .into_par_iter()
            .map(|genome_idx| {
                let cmps = self.compare_genome(genome_idx);
                pb.inc(1);
                cmps
            })
            .collect();

        pb.finish();

        let (lower, upper): (Vec<_>, Vec<_>) = comparisons.into_iter().unzip();
        (lower.into_iter().flatten().collect(), upper.into_iter().flatten().collect())
    }

    /// Compare each locus of the focal genome to its closest neighbours with the locus, where
    /// neighbours are ranked by AAI calculated with the focal locus left out (lower) and by
    /// identity of the locus (upper).
    fn compare_genome(&self, genome_idx: usize) -> (Vec<LocusComparison>, Vec<LocusComparison>) {
        let ls = self.locus_set;
        let num_genomes = ls.num_genomes();
        let num_loci = ls.loci.len();

        // per-locus (identical, aligned) counts to all other genomes, indexed by [j * num_loci + locus_idx],
        // along with counts summed across all loci
        let mut locus_counts = vec![(0u32, 0u32); num_genomes * num_loci];
        let mut total_counts = vec![(0u32, 0u32); num_genomes];
        for j in (0..num_genomes).filter(|&j| j != genome_idx) {
            for locus_idx in 0..num_loci {
                if !ls.present[[genome_idx, locus_idx]] || !ls.present[[j, locus_idx]] {
                    continue;
                }

                let (identical, aligned) = identity_counts(
                    ls.locus_slice(genome_idx, locus_idx),
                    ls.locus_slice(j, locus_idx),
                );
                locus_counts[j * num_loci + locus_idx] = (identical, aligned);
                total_counts[j].0 += identical;
                total_counts[j].1 += aligned;
            }
        }

        let mut lower = Vec::new();
        let mut upper = Vec::new();
        let mut candidates: Vec<Candidate> = Vec::with_capacity(num_genomes);
        let mut by_identity: Vec<Candidate> = Vec::with_capacity(num_genomes);
        for locus_idx in 0..num_loci {
            if !ls.present[[genome_idx, locus_idx]] {
                continue;
            }

            let focal_nongap = ls.nongap[[genome_idx, locus_idx]];
            let focal_nongap_loo = ls.nongap_total[genome_idx] - focal_nongap;

            // calculate AAI excluding the focal locus to all genomes with the locus and sufficient overlap
            candidates.clear();
            for j in (0..num_genomes).filter(|&j| j != genome_idx && ls.present[[j, locus_idx]]) {
                let (locus_identical, locus_aligned) = locus_counts[j * num_loci + locus_idx];
                let nongap = ls.nongap[[j, locus_idx]];
                if locus_aligned < self.overlap.min_aligned_columns
                    || !sufficient_overlap(locus_aligned, focal_nongap, nongap, self.overlap.min_aligned_frac)
                {
                    continue;
                }

                let loo_identical = total_counts[j].0 - locus_identical;
                let loo_aligned = total_counts[j].1 - locus_aligned;
                let nongap_loo = ls.nongap_total[j] - nongap;
                if !sufficient_overlap(loo_aligned, focal_nongap_loo, nongap_loo, self.overlap.min_aligned_frac) {
                    continue;
                }

                candidates.push(Candidate {
                    genome_idx: j,
                    aai_loo: 100.0 * loo_identical as f64 / loo_aligned as f64,
                    locus_identity: 100.0 * locus_identical as f64 / locus_aligned as f64,
                    locus_aligned,
                });
            }

            // select closest neighbours by decreasing locus identity, breaking ties by decreasing AAI
            by_identity.clear();
            by_identity.extend_from_slice(&candidates);
            let cmp = |a: &Candidate, b: &Candidate| {
                b.locus_identity
                    .partial_cmp(&a.locus_identity)
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| b.aai_loo.partial_cmp(&a.aai_loo).unwrap_or(Ordering::Equal))
                    .then_with(|| a.genome_idx.cmp(&b.genome_idx))
            };
            select_top(&mut by_identity, self.upper_neighbours, cmp);
            upper.extend(by_identity.iter().enumerate().map(|(rank, c)| c.comparison(locus_idx, genome_idx, rank)));

            // select closest neighbours by decreasing AAI, using genome index for deterministic ordering
            let cmp = |a: &Candidate, b: &Candidate| {
                b.aai_loo
                    .partial_cmp(&a.aai_loo)
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| a.genome_idx.cmp(&b.genome_idx))
            };
            select_top(&mut candidates, self.neighbours, cmp);
            lower.extend(candidates.iter().enumerate().map(|(rank, c)| c.comparison(locus_idx, genome_idx, rank)));
        }

        (lower, upper)
    }
}

/// Candidate neighbour of a focal genome for a single locus.
#[derive(Clone, Copy)]
struct Candidate {
    genome_idx: usize,
    aai_loo: f64,
    locus_identity: f64,
    locus_aligned: u32,
}

impl Candidate {
    fn comparison(&self, locus_idx: usize, focal_idx: usize, rank: usize) -> LocusComparison {
        LocusComparison {
            locus_idx,
            genome_idx: focal_idx,
            neighbour_idx: self.genome_idx,
            neighbour_rank: rank + 1,
            similarity: self.aai_loo,
            locus_identity: self.locus_identity,
            locus_aligned: self.locus_aligned,
        }
    }
}

/// Reduce candidates to the top `n` under the specified ordering, sorted by that ordering.
fn select_top<F: Fn(&Candidate, &Candidate) -> Ordering>(candidates: &mut Vec<Candidate>, n: usize, cmp: F) {
    if candidates.len() > n {
        candidates.select_nth_unstable_by(n - 1, &cmp);
        candidates.truncate(n);
    }
    candidates.sort_by(&cmp);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_counts() {
        let a = [b'A', b'C', GAP_CODE, b'D', b'E'];
        let b = [b'A', b'D', b'C', GAP_CODE, b'E'];

        // aligned positions: 0, 1, 4; identical: 0, 4
        assert_eq!(identity_counts(&a, &b), (2, 3));
    }

    #[test]
    fn test_select_top_by_identity_breaks_ties_by_aai() {
        let c = |genome_idx, aai_loo, locus_identity| Candidate { genome_idx, aai_loo, locus_identity, locus_aligned: 100 };
        let mut candidates = vec![c(0, 60.0, 100.0), c(1, 95.0, 100.0), c(2, 99.0, 90.0), c(3, 80.0, 100.0)];
        let cmp = |a: &Candidate, b: &Candidate| {
            b.locus_identity
                .partial_cmp(&a.locus_identity)
                .unwrap_or(Ordering::Equal)
                .then_with(|| b.aai_loo.partial_cmp(&a.aai_loo).unwrap_or(Ordering::Equal))
                .then_with(|| a.genome_idx.cmp(&b.genome_idx))
        };
        select_top(&mut candidates, 2, cmp);

        // genomes tied at 100% identity are ranked by AAI, so the focal genome's own lineage is preferred
        let selected: Vec<usize> = candidates.iter().map(|c| c.genome_idx).collect();
        assert_eq!(selected, vec![1, 3]);
    }

    #[test]
    fn test_sufficient_overlap() {
        assert!(sufficient_overlap(50, 100, 200, 0.5));
        assert!(!sufficient_overlap(49, 100, 200, 0.5));
        assert!(!sufficient_overlap(0, 0, 0, 0.0));
    }
}
