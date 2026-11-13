use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result, bail};
use log::{info, warn};
use ndarray::{Array2, s};
use rand::SeedableRng;
use rand::seq::index::sample;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::alignment_reader::AlignmentReader;
use crate::alignment_writer::{check_compressed_input, get_output_prefix, write_filtered_alignment, write_mask};
use crate::prune::trim::marker_info::{MarkerInfo, read_marker_info};
use crate::sequence_type_detector::SequenceTypeDetector;

/// Characters treated as gaps (matches GenomeTreeTk `msa_filter`).
pub const TRIM_GAP_CHARS: [char; 4] = ['-', '.', '_', '*'];

/// Number of adjacent columns summarised together, so each sequence is read in contiguous blocks.
const COL_BLOCK: usize = 64;

fn is_gap(c: char) -> bool {
    TRIM_GAP_CHARS.contains(&c)
}

/// Outcome of filtering a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnClass {
    Retained,
    Gappy,
    Conserved,
    Unconserved,
}

/// Gap count and most common residue of a column.
#[derive(Debug, Clone, PartialEq)]
pub struct ColumnSummary {
    pub num_gaps: usize,
    pub top_residue: Option<char>,
    pub top_residue_count: usize,
}

/// Criteria for retaining a column.
#[derive(Debug, Clone)]
pub struct ColumnFilter {
    max_gap_frac: f64,
    min_consensus_frac: f64,
    max_consensus_frac: f64,
}

impl ColumnFilter {
    /// All arguments are percentages in [0, 100].
    pub fn new(min_perc_aa: f64, min_consensus: f64, max_consensus: f64) -> Self {
        Self {
            max_gap_frac: (100.0 - min_perc_aa) / 100.0,
            min_consensus_frac: min_consensus / 100.0,
            max_consensus_frac: max_consensus / 100.0,
        }
    }

    /// Classify a column. A column is gappy if more than (100 - min_perc_aa)% of sequences have a gap
    /// (or every sequence does). Otherwise the consensus is the frequency of the most common residue
    /// among non-gap characters: below min_consensus is unconserved, at or above max_consensus is
    /// conserved, and anything in between is retained.
    pub fn classify(&self, num_seqs: usize, summary: &ColumnSummary) -> ColumnClass {
        let num_residues = num_seqs - summary.num_gaps;
        if num_residues == 0 || summary.num_gaps as f64 / num_seqs as f64 > self.max_gap_frac {
            return ColumnClass::Gappy;
        }

        let consensus = summary.top_residue_count as f64 / num_residues as f64;
        if consensus < self.min_consensus_frac {
            ColumnClass::Unconserved
        } else if consensus >= self.max_consensus_frac {
            ColumnClass::Conserved
        } else {
            ColumnClass::Retained
        }
    }
}

/// Summarise each column of the alignment, in parallel over blocks of columns.
pub fn summarise_columns(alignment: &Array2<char>) -> Vec<ColumnSummary> {
    let ncols = alignment.ncols();
    let blocks: Vec<usize> = (0..ncols).step_by(COL_BLOCK).collect();

    blocks
        .par_iter()
        .map(|&start| summarise_block(alignment, start, (start + COL_BLOCK).min(ncols)))
        .collect::<Vec<_>>()
        .into_iter()
        .flatten()
        .collect()
}

fn summarise_block(alignment: &Array2<char>, start: usize, end: usize) -> Vec<ColumnSummary> {
    let width = end - start;
    let mut num_gaps = vec![0usize; width];
    let mut ascii_counts = vec![[0u32; 128]; width];
    let mut other_counts: HashMap<(usize, char), usize> = HashMap::new();

    for row in alignment.rows() {
        for (j, &c) in row.slice(s![start..end]).iter().enumerate() {
            if is_gap(c) {
                num_gaps[j] += 1;
            } else if c.is_ascii() {
                ascii_counts[j][c as usize] += 1;
            } else {
                *other_counts.entry((j, c)).or_default() += 1;
            }
        }
    }

    (0..width)
        .map(|j| {
            // most common residue; ties are broken by the lowest character code
            let mut top_residue = None;
            let mut top_residue_count = 0;
            let ascii = ascii_counts[j].iter().enumerate().map(|(code, &n)| (code as u8 as char, n as usize));
            let other = other_counts.iter().filter(|((col, _), _)| *col == j).map(|((_, c), &n)| (*c, n));
            for (c, n) in ascii.chain(other) {
                if n > top_residue_count || (n == top_residue_count && n > 0 && Some(c) < top_residue) {
                    top_residue = Some(c);
                    top_residue_count = n;
                }
            }

            ColumnSummary {
                num_gaps: num_gaps[j],
                top_residue,
                top_residue_count,
            }
        })
        .collect()
}

/// Smallest number of columns to sample from each marker that gives an alignment of at least
/// `target_len` columns. If the target cannot be reached, every retained column is used.
pub fn columns_per_marker(retained_counts: &[usize], target_len: usize) -> usize {
    let max_retained = retained_counts.iter().copied().max().unwrap_or(0);
    (1..=max_retained)
        .find(|&k| retained_counts.iter().map(|&n| n.min(k)).sum::<usize>() >= target_len)
        .unwrap_or(max_retained)
}

/// Filtering results for a single marker.
#[derive(Debug, Clone, Default)]
struct MarkerFilterInfo {
    retained_cols: Vec<usize>,
    sampled_cols: usize,
    gappy: usize,
    conserved: usize,
    unconserved: usize,
}

/// Trim a concatenated marker gene alignment by removing gappy, highly conserved and poorly
/// conserved columns, then randomly subsampling the same number of columns from each marker.
pub struct MsaTrimmer {
    alignment_file: String,
    marker_info_file: String,
    msa_length: usize,
    min_perc_aa: f64,
    min_consensus: f64,
    max_consensus: f64,
    seed: u64,
    out_dir: String,
}

impl MsaTrimmer {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        alignment_file: String,
        marker_info_file: String,
        msa_length: usize,
        min_perc_aa: f64,
        min_consensus: f64,
        max_consensus: f64,
        seed: u64,
        out_dir: String,
    ) -> Self {
        Self {
            alignment_file,
            marker_info_file,
            msa_length,
            min_perc_aa,
            min_consensus,
            max_consensus,
            seed,
            out_dir,
        }
    }

    pub fn run(&mut self) -> Result<()> {
        info!("Starting MSA trimming workflow.");

        if self.min_consensus >= self.max_consensus {
            bail!(
                "Minimum consensus ({}) must be less than maximum consensus ({}).",
                self.min_consensus,
                self.max_consensus
            );
        }

        info!("Reading marker info file: {}", self.marker_info_file);
        let markers = read_marker_info(&self.marker_info_file)?;
        let concat_len: usize = markers.iter().map(|m| m.length).sum();
        info!("Read {} markers with a concatenated length of {} columns.", markers.len(), concat_len);

        info!("Reading alignment from file: {}", self.alignment_file);
        let reader = AlignmentReader::new(self.alignment_file.clone());
        let (seq_ids, alignment) = reader.run()?;
        info!("Alignment contains {} sequences and {} columns.", seq_ids.len(), alignment.ncols());

        if alignment.ncols() != concat_len {
            bail!(
                "Alignment length ({}) does not equal the summed length of the markers in the marker info file ({}).",
                alignment.ncols(),
                concat_len
            );
        }

        let (is_dna, alphabet) = SequenceTypeDetector::detect(&alignment)?;
        info!("Sequence type detected: {}", if is_dna { "DNA" } else { "amino acid" });

        // filter columns
        info!(
            "Filtering columns with < {}% residues, < {}% consensus, or >= {}% consensus.",
            self.min_perc_aa, self.min_consensus, self.max_consensus
        );
        let summaries = summarise_columns(&alignment);
        let filter = ColumnFilter::new(self.min_perc_aa, self.min_consensus, self.max_consensus);
        let classes: Vec<ColumnClass> = summaries.iter().map(|s| filter.classify(alignment.nrows(), s)).collect();
        self.report_nonstandard_consensus(&summaries, &classes, &alphabet);

        let mut filter_info = Self::filter_markers(&markers, &classes);
        let retained_counts: Vec<usize> = filter_info.iter().map(|f| f.retained_cols.len()).collect();
        let total_retained: usize = retained_counts.iter().sum();
        info!("Gappy columns: {}", filter_info.iter().map(|f| f.gappy).sum::<usize>());
        info!("Conserved columns: {}", filter_info.iter().map(|f| f.conserved).sum::<usize>());
        info!("Unconserved columns: {}", filter_info.iter().map(|f| f.unconserved).sum::<usize>());
        info!(
            "Retained {} columns (min: {}, max: {} per marker) after filtering.",
            total_retained,
            retained_counts.iter().min().unwrap_or(&0),
            retained_counts.iter().max().unwrap_or(&0)
        );

        if total_retained == 0 {
            bail!("No columns passed the filtering criteria.");
        }

        // randomly subsample columns from each marker
        info!("Randomly subsampling columns to obtain an MSA of ~{} columns (seed = {}).", self.msa_length, self.seed);
        let cols_per_marker = columns_per_marker(&retained_counts, self.msa_length);
        let sampled_cols = self.subsample_columns(&mut filter_info, cols_per_marker);
        info!("Sampling up to {} columns from each marker.", cols_per_marker);
        info!(
            "{} markers had fewer than {} columns available for sampling.",
            retained_counts.iter().filter(|&&n| n < cols_per_marker).count(),
            cols_per_marker
        );
        info!("Trimmed MSA has {} columns.", sampled_cols.len());

        let filtered_mask: Vec<bool> = classes.iter().map(|&c| c == ColumnClass::Retained).collect();
        let mut final_mask = vec![false; alignment.ncols()];
        for &col_idx in &sampled_cols {
            final_mask[col_idx] = true;
        }

        // write results
        let output_prefix = get_output_prefix(&self.alignment_file, &self.out_dir)?;
        let output_suffix = format!(
            "_l{}_aa{}_c{}-{}_trimmed.fasta",
            self.msa_length, self.min_perc_aa, self.min_consensus, self.max_consensus
        );
        let output_file = |ext: &str| format!("{}{}", output_prefix, output_suffix.replace(".fasta", ext));

        write_mask(&filtered_mask, &output_file("_filtered.mask"))?;
        write_mask(&final_mask, &output_file(".mask"))?;
        self.write_marker_stats(&markers, &filter_info, &output_file("_markers.tsv"))?;
        self.write_genome_stats(&seq_ids, &alignment, &sampled_cols, &output_file("_genomes.tsv"))?;

        info!("Writing trimmed alignment to file.");
        let output_alignment_file = format!("{}{}", output_prefix, output_suffix);
        let compressed = check_compressed_input(&self.alignment_file);
        write_filtered_alignment(&seq_ids, &alignment, &final_mask, &output_alignment_file, compressed)?;

        Ok(())
    }

    /// Assign column classes to markers.
    fn filter_markers(markers: &[MarkerInfo], classes: &[ColumnClass]) -> Vec<MarkerFilterInfo> {
        let mut start = 0;
        markers
            .iter()
            .map(|marker| {
                let mut info = MarkerFilterInfo::default();
                for (col_idx, class) in classes.iter().enumerate().skip(start).take(marker.length) {
                    match class {
                        ColumnClass::Retained => info.retained_cols.push(col_idx),
                        ColumnClass::Gappy => info.gappy += 1,
                        ColumnClass::Conserved => info.conserved += 1,
                        ColumnClass::Unconserved => info.unconserved += 1,
                    }
                }
                start += marker.length;
                info
            })
            .collect()
    }

    /// Randomly sample up to `cols_per_marker` retained columns from each marker. Returns the
    /// sampled column indices in alignment order.
    fn subsample_columns(&self, filter_info: &mut [MarkerFilterInfo], cols_per_marker: usize) -> Vec<usize> {
        let mut rng = ChaCha8Rng::seed_from_u64(self.seed);
        let mut sampled = BTreeSet::new();
        for info in filter_info.iter_mut() {
            let num_retained = info.retained_cols.len();
            let amount = cols_per_marker.min(num_retained);
            for idx in sample(&mut rng, num_retained, amount).iter() {
                sampled.insert(info.retained_cols[idx]);
            }
            info.sampled_cols = amount;
        }

        sampled.into_iter().collect()
    }

    /// Warn about columns whose most common residue is not a standard residue.
    fn report_nonstandard_consensus(&self, summaries: &[ColumnSummary], classes: &[ColumnClass], alphabet: &[char]) {
        let mut nonstandard: HashMap<char, usize> = HashMap::new();
        for (summary, class) in summaries.iter().zip(classes) {
            if *class == ColumnClass::Gappy {
                continue;
            }
            if let Some(c) = summary.top_residue.filter(|c| !alphabet.contains(c)) {
                *nonstandard.entry(c).or_default() += 1;
            }
        }

        if !nonstandard.is_empty() {
            let mut counts: Vec<(char, usize)> = nonstandard.into_iter().collect();
            counts.sort();
            let total: usize = counts.iter().map(|(_, n)| n).sum();
            let counts: Vec<String> = counts.iter().map(|(c, n)| format!("{}:{}", c, n)).collect();
            warn!(
                "Most common residue was not in the standard alphabet in {} columns ({}).",
                total,
                counts.join(", ")
            );
        }
    }

    /// Write per-marker filtering statistics.
    fn write_marker_stats(&self, markers: &[MarkerInfo], filter_info: &[MarkerFilterInfo], out_file: &str) -> Result<()> {
        let file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file))?;
        let mut writer = BufWriter::new(file);

        writeln!(
            writer,
            "marker_id\tdescription\tlength\tretained_columns\tretained_columns_perc\tsampled_columns\tgappy_columns\tconserved_columns\tunconserved_columns"
        )?;
        for (marker, info) in markers.iter().zip(filter_info) {
            let retained_perc = if marker.length > 0 {
                100.0 * info.retained_cols.len() as f64 / marker.length as f64
            } else {
                0.0
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
                marker.id,
                marker.description,
                marker.length,
                info.retained_cols.len(),
                retained_perc,
                info.sampled_cols,
                info.gappy,
                info.conserved,
                info.unconserved
            )?;
        }
        writer.flush()?;

        Ok(())
    }

    /// Write the number of residues each sequence has in the trimmed alignment.
    fn write_genome_stats(&self, seq_ids: &[String], alignment: &Array2<char>, sampled_cols: &[usize], out_file: &str) -> Result<()> {
        let residue_counts: Vec<usize> = (0..alignment.nrows())
            .into_par_iter()
            .map(|row_idx| {
                let row = alignment.row(row_idx);
                sampled_cols.iter().filter(|&&col_idx| row[col_idx].is_ascii_alphabetic()).count()
            })
            .collect();

        let file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file))?;
        let mut writer = BufWriter::new(file);

        let msa_len = sampled_cols.len();
        writeln!(writer, "genome_id\tmsa_length\tresidues\tresidues_perc")?;
        for (seq_id, &num_residues) in seq_ids.iter().zip(&residue_counts) {
            writeln!(
                writer,
                "{}\t{}\t{}\t{:.2}",
                seq_id,
                msa_len,
                num_residues,
                100.0 * num_residues as f64 / msa_len as f64
            )?;
        }
        writer.flush()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;

    fn summary(num_gaps: usize, top_residue_count: usize) -> ColumnSummary {
        ColumnSummary {
            num_gaps,
            top_residue: Some('A'),
            top_residue_count,
        }
    }

    #[test]
    fn test_summarise_columns() {
        let alignment = arr2(&[
            ['A', '-', 'C', '*'],
            ['A', '.', 'D', '_'],
            ['C', 'K', 'D', '-'],
            ['A', 'K', 'C', 'X'],
        ]);
        let summaries = summarise_columns(&alignment);

        assert_eq!(summaries[0], ColumnSummary { num_gaps: 0, top_residue: Some('A'), top_residue_count: 3 });
        assert_eq!(summaries[1], ColumnSummary { num_gaps: 2, top_residue: Some('K'), top_residue_count: 2 });
        // tie between C and D is broken by the lowest character code
        assert_eq!(summaries[2], ColumnSummary { num_gaps: 0, top_residue: Some('C'), top_residue_count: 2 });
        assert_eq!(summaries[3], ColumnSummary { num_gaps: 3, top_residue: Some('X'), top_residue_count: 1 });
    }

    #[test]
    fn test_summarise_columns_multiple_blocks() {
        let ncols = 2 * COL_BLOCK + 5;
        let alignment = Array2::from_shape_fn((3, ncols), |(r, c)| if r == 0 && c % 2 == 0 { '-' } else { 'L' });
        let summaries = summarise_columns(&alignment);

        assert_eq!(summaries.len(), ncols);
        for (c, s) in summaries.iter().enumerate() {
            let expected_gaps = if c % 2 == 0 { 1 } else { 0 };
            assert_eq!(s.num_gaps, expected_gaps);
            assert_eq!(s.top_residue_count, 3 - expected_gaps);
        }
    }

    #[test]
    fn test_classify() {
        // defaults: >= 50% residues, consensus in [25%, 95%)
        let filter = ColumnFilter::new(50.0, 25.0, 95.0);

        assert_eq!(filter.classify(10, &summary(5, 3)), ColumnClass::Retained); // 50% gaps is allowed
        assert_eq!(filter.classify(10, &summary(6, 3)), ColumnClass::Gappy);
        assert_eq!(filter.classify(10, &summary(10, 0)), ColumnClass::Gappy);
        assert_eq!(filter.classify(8, &summary(0, 1)), ColumnClass::Unconserved); // 12.5%
        assert_eq!(filter.classify(8, &summary(0, 2)), ColumnClass::Retained); // 25% is retained
        assert_eq!(filter.classify(20, &summary(0, 18)), ColumnClass::Retained); // 90%
        assert_eq!(filter.classify(20, &summary(0, 19)), ColumnClass::Conserved); // 95% is conserved
        assert_eq!(filter.classify(20, &summary(0, 20)), ColumnClass::Conserved);

        // consensus is calculated over non-gap residues
        assert_eq!(filter.classify(20, &summary(10, 10)), ColumnClass::Conserved);

        // all-gap column is gappy even when gaps are allowed
        let no_gap_filter = ColumnFilter::new(0.0, 25.0, 95.0);
        assert_eq!(no_gap_filter.classify(4, &summary(4, 0)), ColumnClass::Gappy);
    }

    #[test]
    fn test_columns_per_marker() {
        assert_eq!(columns_per_marker(&[10, 10, 10], 3), 1);
        assert_eq!(columns_per_marker(&[10, 10, 10], 4), 2);
        assert_eq!(columns_per_marker(&[10, 10, 10], 30), 10);
        assert_eq!(columns_per_marker(&[2, 10, 10], 12), 5);
        // target cannot be reached: use every retained column
        assert_eq!(columns_per_marker(&[2, 10, 10], 100), 10);
        assert_eq!(columns_per_marker(&[0, 0], 100), 0);
    }

    #[test]
    fn test_filter_and_subsample() {
        let markers = vec![
            MarkerInfo { id: "M1".into(), description: "".into(), length: 4 },
            MarkerInfo { id: "M2".into(), description: "".into(), length: 3 },
        ];
        use ColumnClass::*;
        let classes = vec![Retained, Gappy, Retained, Retained, Conserved, Retained, Unconserved];
        let mut filter_info = MsaTrimmer::filter_markers(&markers, &classes);

        assert_eq!(filter_info[0].retained_cols, vec![0, 2, 3]);
        assert_eq!(filter_info[0].gappy, 1);
        assert_eq!(filter_info[1].retained_cols, vec![5]);
        assert_eq!((filter_info[1].conserved, filter_info[1].unconserved), (1, 1));

        let trimmer = MsaTrimmer::new(String::new(), String::new(), 3, 50.0, 25.0, 95.0, 1, String::new());
        let sampled = trimmer.subsample_columns(&mut filter_info, 2);
        assert_eq!(sampled.len(), 3);
        assert!(sampled.contains(&5));
        assert!(sampled.windows(2).all(|w| w[0] < w[1]));
        assert_eq!((filter_info[0].sampled_cols, filter_info[1].sampled_cols), (2, 1));

        // sampling is deterministic for a given seed
        assert_eq!(sampled, trimmer.subsample_columns(&mut filter_info, 2));
    }
}
