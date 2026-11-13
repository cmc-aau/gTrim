use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{Context, Result, bail};
use log::{info, warn};
use ndarray::Array2;
use rand::SeedableRng;
use rand::seq::index::sample;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::alignment_reader::AlignmentReader;
use crate::column_counts::map_column_counts;
use crate::alignment_writer::{check_compressed_input, get_output_prefix, write_filtered_alignment, write_mask};
use crate::prune::trim::column_mask::read_column_mask;
use crate::prune::trim::marker_info::{MarkerInfo, read_marker_info};
use crate::sequence_type_detector::SequenceTypeDetector;

/// Characters treated as gaps (matches GenomeTreeTk `msa_filter`).
pub const TRIM_GAP_CHARS: [u8; 4] = *b"-._*";

fn is_gap(c: u8) -> bool {
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
    pub top_residue: Option<u8>,
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

/// Summarise each column of the alignment.
pub fn summarise_columns(alignment: &Array2<u8>) -> Vec<ColumnSummary> {
    map_column_counts(alignment, |counts| {
        let num_gaps = TRIM_GAP_CHARS.iter().map(|&g| counts[g as usize] as usize).sum();

        // most common residue; ties are broken by the lowest character code
        let mut top_residue = None;
        let mut top_residue_count = 0;
        for (code, &n) in counts.iter().enumerate() {
            let n = n as usize;
            if n > top_residue_count && !is_gap(code as u8) {
                top_residue = Some(code as u8);
                top_residue_count = n;
            }
        }

        ColumnSummary {
            num_gaps,
            top_residue,
            top_residue_count,
        }
    })
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
    /// Columns removed by the input mask before trimming.
    masked: usize,
    /// Indices of retained columns in the input alignment.
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
    mask_file: Option<String>,
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
        mask_file: Option<String>,
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
            mask_file,
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

        // columns of the full concatenated alignment present in the input alignment
        let full_mask = self.read_input_mask(&markers, concat_len)?;
        let input_cols: Vec<usize> = (0..concat_len).filter(|&i| full_mask[i]).collect();
        let marker_lens = Self::masked_marker_lengths(&markers, &full_mask);

        info!("Reading alignment from file: {}", self.alignment_file);
        let reader = AlignmentReader::new(self.alignment_file.clone());
        let (seq_ids, alignment) = reader.run()?;
        info!("Alignment contains {} sequences and {} columns.", seq_ids.len(), alignment.ncols());

        if self.mask_file.is_none() && alignment.ncols() != concat_len {
            bail!(
                "Alignment length ({}) does not equal the summed length of the markers in the marker info file ({}).",
                alignment.ncols(),
                concat_len
            );
        } else if alignment.ncols() != input_cols.len() {
            bail!(
                "Alignment length ({}) does not equal the number of columns retained by the mask file ({} of {}). The alignment must contain only the retained columns of the full concatenated alignment.",
                alignment.ncols(),
                input_cols.len(),
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

        let mut filter_info = Self::filter_markers(&marker_lens, &classes);
        for ((info, marker), &len) in filter_info.iter_mut().zip(&markers).zip(&marker_lens) {
            info.masked = marker.length - len;
        }
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
        let cols_per_marker = columns_per_marker(&retained_counts, self.msa_length);
        let sampled_cols = self.subsample_columns(&mut filter_info, cols_per_marker);
        if total_retained <= self.msa_length {
            info!(
                "MSA length ({}) is not less than the number of retained columns; retaining all {} columns.",
                self.msa_length, total_retained
            );
        } else {
            info!("Randomly subsampling columns to obtain an MSA of ~{} columns (seed = {}).", self.msa_length, self.seed);
            info!("Sampling up to {} columns from each marker.", cols_per_marker);
            info!(
                "{} markers had fewer than {} columns available for sampling.",
                retained_counts.iter().filter(|&&n| n < cols_per_marker).count(),
                cols_per_marker
            );
        }
        info!("Trimmed MSA has {} columns.", sampled_cols.len());

        // masks over the input alignment
        let mut final_mask = vec![false; alignment.ncols()];
        for &col_idx in &sampled_cols {
            final_mask[col_idx] = true;
        }

        // masks over the full concatenated alignment, which are written to file
        let mut full_filtered_mask = vec![false; concat_len];
        for (col_idx, &class) in classes.iter().enumerate() {
            if class == ColumnClass::Retained {
                full_filtered_mask[input_cols[col_idx]] = true;
            }
        }
        let mut full_final_mask = vec![false; concat_len];
        for &col_idx in &sampled_cols {
            full_final_mask[input_cols[col_idx]] = true;
        }

        // write results
        let output_prefix = format!("{}.trimmed", get_output_prefix(&self.alignment_file, &self.out_dir)?);
        let output_file = |ext: &str| format!("{}{}", output_prefix, ext);

        write_mask(&full_filtered_mask, &output_file(".filtered.mask"))?;
        write_mask(&full_final_mask, &output_file(".mask"))?;
        self.write_marker_stats(&markers, &filter_info, &output_file(".markers.tsv"))?;
        self.write_genome_stats(&seq_ids, &alignment, &sampled_cols, &output_file(".genomes.tsv"))?;

        info!("Writing trimmed alignment to file.");
        let output_alignment_file = output_file(".fasta");
        let compressed = check_compressed_input(&self.alignment_file);
        write_filtered_alignment(&seq_ids, &alignment, &final_mask, &output_alignment_file, compressed)?;

        Ok(())
    }

    /// Read the mask file, if any, and return which columns of the full concatenated alignment are
    /// in the input alignment.
    fn read_input_mask(&self, markers: &[MarkerInfo], concat_len: usize) -> Result<Vec<bool>> {
        let Some(mask_file) = &self.mask_file else {
            return Ok(vec![true; concat_len]);
        };

        info!("Reading mask file: {}", mask_file);
        let mask = read_column_mask(mask_file, concat_len)?;
        let marker_lens = Self::masked_marker_lengths(markers, &mask.retained);
        let num_retained: usize = marker_lens.iter().sum();
        let markers_masked = markers.iter().zip(&marker_lens).filter(|(m, len)| m.length != **len).count();
        info!(
            "Mask file ({} format) retains {} of {} columns; {} columns were removed from {} markers.",
            mask.format,
            num_retained,
            concat_len,
            concat_len - num_retained,
            markers_masked
        );

        let empty: Vec<&str> = markers
            .iter()
            .zip(&marker_lens)
            .filter(|(m, len)| m.length > 0 && **len == 0)
            .map(|(m, _)| m.id.as_str())
            .collect();
        if !empty.is_empty() {
            warn!("{} markers have no columns left after masking: {}", empty.len(), empty.join(", "));
        }

        Ok(mask.retained)
    }

    /// Number of columns of each marker retained by a mask over the full concatenated alignment.
    fn masked_marker_lengths(markers: &[MarkerInfo], full_mask: &[bool]) -> Vec<usize> {
        let mut start = 0;
        markers
            .iter()
            .map(|marker| {
                let len = full_mask[start..start + marker.length].iter().filter(|&&r| r).count();
                start += marker.length;
                len
            })
            .collect()
    }

    /// Assign column classes to markers, given the length of each marker in the input alignment.
    fn filter_markers(marker_lens: &[usize], classes: &[ColumnClass]) -> Vec<MarkerFilterInfo> {
        let mut start = 0;
        marker_lens
            .iter()
            .map(|&len| {
                let mut info = MarkerFilterInfo::default();
                for (col_idx, class) in classes.iter().enumerate().skip(start).take(len) {
                    match class {
                        ColumnClass::Retained => info.retained_cols.push(col_idx),
                        ColumnClass::Gappy => info.gappy += 1,
                        ColumnClass::Conserved => info.conserved += 1,
                        ColumnClass::Unconserved => info.unconserved += 1,
                    }
                }
                start += len;
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
    fn report_nonstandard_consensus(&self, summaries: &[ColumnSummary], classes: &[ColumnClass], alphabet: &[u8]) {
        let mut nonstandard: HashMap<u8, usize> = HashMap::new();
        for (summary, class) in summaries.iter().zip(classes) {
            if *class == ColumnClass::Gappy {
                continue;
            }
            if let Some(c) = summary.top_residue.filter(|c| !alphabet.contains(c)) {
                *nonstandard.entry(c).or_default() += 1;
            }
        }

        if !nonstandard.is_empty() {
            let mut counts: Vec<(u8, usize)> = nonstandard.into_iter().collect();
            counts.sort();
            let total: usize = counts.iter().map(|(_, n)| n).sum();
            let counts: Vec<String> = counts.iter().map(|&(c, n)| format!("{}:{}", c as char, n)).collect();
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
            "marker_id\tdescription\tlength\tmasked_columns\tretained_columns\tretained_columns_perc\tsampled_columns\tgappy_columns\tconserved_columns\tunconserved_columns"
        )?;
        for (marker, info) in markers.iter().zip(filter_info) {
            let retained_perc = if marker.length > 0 {
                100.0 * info.retained_cols.len() as f64 / marker.length as f64
            } else {
                0.0
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
                marker.id,
                marker.description,
                marker.length,
                info.masked,
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
    fn write_genome_stats(&self, seq_ids: &[String], alignment: &Array2<u8>, sampled_cols: &[usize], out_file: &str) -> Result<()> {
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
            top_residue: Some(b'A'),
            top_residue_count,
        }
    }

    #[test]
    fn test_summarise_columns() {
        let alignment = arr2(&[*b"A-C*", *b"A.D_", *b"CKD-", *b"AKCX"]);
        let summaries = summarise_columns(&alignment);

        assert_eq!(summaries[0], ColumnSummary { num_gaps: 0, top_residue: Some(b'A'), top_residue_count: 3 });
        assert_eq!(summaries[1], ColumnSummary { num_gaps: 2, top_residue: Some(b'K'), top_residue_count: 2 });
        // tie between C and D is broken by the lowest character code
        assert_eq!(summaries[2], ColumnSummary { num_gaps: 0, top_residue: Some(b'C'), top_residue_count: 2 });
        assert_eq!(summaries[3], ColumnSummary { num_gaps: 3, top_residue: Some(b'X'), top_residue_count: 1 });
    }

    #[test]
    fn test_summarise_columns_multiple_blocks() {
        use crate::column_counts::COL_BLOCK;

        let ncols = 2 * COL_BLOCK + 5;
        let alignment = Array2::from_shape_fn((3, ncols), |(r, c)| if r == 0 && c % 2 == 0 { b'-' } else { b'L' });
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
        use ColumnClass::*;
        let classes = vec![Retained, Gappy, Retained, Retained, Conserved, Retained, Unconserved];
        let mut filter_info = MsaTrimmer::filter_markers(&[4, 3], &classes);

        assert_eq!(filter_info[0].retained_cols, vec![0, 2, 3]);
        assert_eq!(filter_info[0].gappy, 1);
        assert_eq!(filter_info[1].retained_cols, vec![5]);
        assert_eq!((filter_info[1].conserved, filter_info[1].unconserved), (1, 1));

        let trimmer = MsaTrimmer::new(String::new(), String::new(), None, 3, 50.0, 25.0, 95.0, 1, String::new());
        let sampled = trimmer.subsample_columns(&mut filter_info, 2);
        assert_eq!(sampled.len(), 3);
        assert!(sampled.contains(&5));
        assert!(sampled.windows(2).all(|w| w[0] < w[1]));
        assert_eq!((filter_info[0].sampled_cols, filter_info[1].sampled_cols), (2, 1));

        // sampling is deterministic for a given seed
        assert_eq!(sampled, trimmer.subsample_columns(&mut filter_info, 2));
    }

    #[test]
    fn test_masked_marker_lengths() {
        let markers = vec![
            MarkerInfo { id: "M1".into(), description: "".into(), length: 4 },
            MarkerInfo { id: "M2".into(), description: "".into(), length: 3 },
            MarkerInfo { id: "M3".into(), description: "".into(), length: 2 },
        ];
        // the first column of M2 is removed, which GenomeTreeTk msa_filter attributes to M1
        let mask = [true, false, true, true, false, true, true, false, false];
        assert_eq!(MsaTrimmer::masked_marker_lengths(&markers, &mask), vec![3, 2, 0]);
        assert_eq!(MsaTrimmer::masked_marker_lengths(&markers, &[true; 9]), vec![4, 3, 2]);
    }
}
