use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use log::info;
use needletail::parse_fastx_file;

use crate::alignment_writer::{check_compressed_input, create_writer};
use super::aai::{genome_aai, LocusComparison, NeighbourComparer, OverlapFilter};
use super::locus_set::{genome_id_from_header, LocusSet};
use super::outliers::{expected_values, tally_votes, OutlierMethod, VoteCounts};
use super::plots::{write_plots, AxisScale, LocusPlotData, PlotSettings, PointClass};
use super::taxonomy::Taxonomy;

/// Test used to identify outlier loci.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutlierTest {
    /// Neighbours selected by genome-wide AAI. Outlier if locus identity is unexpectedly low.
    LowIdentity,
    /// Neighbours selected by locus identity. Outlier if genome-wide AAI is unexpectedly low, i.e. the
    /// locus is unexpectedly similar to genomes from another lineage (e.g. contamination or lateral transfer).
    HighIdentity,
}

impl OutlierTest {
    /// Value used to order comparisons when fitting the trend.
    fn x(&self, c: &LocusComparison) -> f64 {
        match self {
            OutlierTest::LowIdentity => c.similarity,
            OutlierTest::HighIdentity => c.locus_identity,
        }
    }

    /// Value tested against its expected value from the trend.
    fn y(&self, c: &LocusComparison) -> f64 {
        match self {
            OutlierTest::LowIdentity => c.locus_identity,
            OutlierTest::HighIdentity => c.similarity,
        }
    }

    fn description(&self) -> &'static str {
        match self {
            OutlierTest::LowIdentity => "low identity",
            OutlierTest::HighIdentity => "high identity",
        }
    }
}

/// Locus comparisons with their expected values, outlier status, and outlier votes for an outlier test.
struct TestResults {
    test: OutlierTest,
    comparisons: Vec<LocusComparison>,
    expected: Vec<f64>,
    outliers: Vec<bool>,
    votes: BTreeMap<(usize, usize), VoteCounts>,
}

impl TestResults {
    fn is_masked(&self, genome_idx: usize, locus_idx: usize) -> bool {
        self.votes.get(&(genome_idx, locus_idx)).is_some_and(|v| v.is_masked())
    }

    fn num_masked(&self) -> usize {
        self.votes.values().filter(|v| v.is_masked()).count()
    }
}

/// Minimum genome-wide AAI (%) between all neighbours responsible for a high-identity mask for them to be
/// reported as a single source lineage.
pub const SOURCE_COHERENT_AAI: f64 = 90.0;

/// Neighbours of a genome locus ranked by locus identity, with their high-identity outlier status.
type IdentityMatches = HashMap<(usize, usize), Vec<(LocusComparison, bool)>>;

/// Minimum number of masked loci for a genome to be flagged as suspect.
pub const MIN_SUSPECT_MASKED: usize = 3;

/// Neighbours by locus identity of each masked genome locus, in rank order, with their high-identity outlier
/// status (false if the high-identity test wasn't performed).
fn identity_matches(comparisons: &[LocusComparison], outliers: Option<&[bool]>, masked: &BTreeSet<(usize, usize)>) -> IdentityMatches {
    let mut matches: IdentityMatches = HashMap::new();
    for (i, c) in comparisons.iter().enumerate() {
        let key = (c.genome_idx, c.locus_idx);
        if masked.contains(&key) {
            matches.entry(key).or_default().push((c.clone(), outliers.is_some_and(|o| o[i])));
        }
    }

    for ranked in matches.values_mut() {
        ranked.sort_by_key(|(c, _)| c.neighbour_rank);
    }

    matches
}

/// Main source of each genome's high-identity masks: the most common best match among them (ties: first genome
/// in ID order), and the number of those masks whose best match is in the same lineage, i.e. is the main source or
/// has at least SOURCE_COHERENT_AAI genome-wide AAI to it. Genomes without high-identity masks are omitted.
fn main_sources(locus_set: &LocusSet, upper: &TestResults, matches: &IdentityMatches) -> HashMap<usize, (usize, usize)> {
    let mut best_matches: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (key, v) in &upper.votes {
        if let Some((c, _)) = matches.get(key).and_then(|m| m.first()).filter(|_| v.is_masked()) {
            best_matches.entry(key.0).or_default().push(c.neighbour_idx);
        }
    }

    best_matches
        .into_iter()
        .map(|(genome_idx, sources)| {
            let same_lineage = |a: usize, b: usize| genome_aai(locus_set, a, b).is_some_and(|aai| aai >= SOURCE_COHERENT_AAI);
            (genome_idx, main_source(&sources, same_lineage))
        })
        .collect()
}

/// Most common source genome (ties: lowest index) and the number of sources that are it or in its lineage.
fn main_source(sources: &[usize], same_lineage: impl Fn(usize, usize) -> bool) -> (usize, usize) {
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    for &s in sources {
        *counts.entry(s).or_default() += 1;
    }
    let (&main, _) = counts
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
        .expect("at least one source");
    let num_in_lineage = sources.iter().filter(|&&s| s == main || same_lineage(main, s)).count();

    (main, num_in_lineage)
}

/// Check if a genome has enough masked loci to be suspect (e.g. contaminated or chimeric): at least
/// `suspect_masked_perc` percent of its loci, and at least MIN_SUSPECT_MASKED loci.
fn is_suspect(num_present: usize, num_masked: usize, suspect_masked_perc: f64) -> bool {
    num_masked >= MIN_SUSPECT_MASKED && 100.0 * num_masked as f64 >= suspect_masked_perc * num_present as f64
}

/// Parameters of an outlier test.
#[derive(Debug, Clone, Copy)]
pub struct TestParams {
    /// Number of closest neighbours compared to each genome locus.
    pub neighbours: usize,
    /// Number of nearest observations used to calculate the running median.
    pub window_size: usize,
    /// Percentage points below the expected value for a comparison to be an outlier.
    pub outlier_threshold: f64,
}

pub struct DivergentLociMasker {
    msa_dir: String,
    method: OutlierMethod,
    lower: TestParams,
    upper: Option<TestParams>,
    overlap: OverlapFilter,
    suspect_masked_perc: f64,
    taxonomy_file: Option<String>,
    plot: bool,
    out_dir: String,
}

impl DivergentLociMasker {
    /// Create masker. The high-identity test is performed if `upper` parameters are given. Genomes with at least
    /// `suspect_masked_perc` percent of their loci masked are flagged as suspect in the genome summary. If a
    /// taxonomy file is given, the lowest common taxon of the genomes most similar to each masked locus is reported.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        msa_dir: String,
        method: OutlierMethod,
        lower: TestParams,
        upper: Option<TestParams>,
        overlap: OverlapFilter,
        suspect_masked_perc: f64,
        taxonomy_file: Option<String>,
        plot: bool,
        out_dir: String,
    ) -> Self {
        Self {
            msa_dir,
            method,
            lower,
            upper,
            overlap,
            suspect_masked_perc,
            taxonomy_file,
            plot,
            out_dir,
        }
    }

    pub fn run(&mut self) -> Result<()> {
        info!("Starting divergent loci masking workflow.");

        info!("Reading MSAs from directory: {}", self.msa_dir);
        let locus_set = LocusSet::from_dir(&self.msa_dir)?;
        info!("Sequence type detected: {}", if locus_set.is_dna { "DNA" } else { "amino acid" });

        let taxonomy = match &self.taxonomy_file {
            Some(file) => {
                let taxonomy = Taxonomy::from_file(file)?;
                let num_classified = locus_set.genome_ids.iter().filter(|id| taxonomy.contains(id)).count();
                info!("Read taxonomy for {} of {} genomes from: {}", num_classified, locus_set.num_genomes(), file);
                Some(taxonomy)
            }
            None => None,
        };

        info!("Comparing loci to {} closest neighbour(s) by AAI:", self.lower.neighbours);
        if let Some(upper) = &self.upper {
            info!("Also comparing loci to {} closest neighbour(s) by locus identity.", upper.neighbours);
        }
        // neighbours are always ranked by locus identity so the most similar copy of each masked locus can be
        // reported; without the high-identity test only the most similar one is needed
        let comparer = NeighbourComparer::new(
            &locus_set,
            self.lower.neighbours,
            self.upper.map_or(1, |p| p.neighbours),
            self.overlap,
        );
        let (lower_comparisons, upper_comparisons) = comparer.run();

        let lower = self.evaluate(&locus_set, lower_comparisons, OutlierTest::LowIdentity, &self.lower);
        let (upper, unevaluated_upper) = match self.upper {
            Some(params) => (Some(self.evaluate(&locus_set, upper_comparisons, OutlierTest::HighIdentity, &params)), Vec::new()),
            None => (None, upper_comparisons),
        };

        let num_present = locus_set.present.iter().filter(|&&p| p).count();
        if lower.votes.len() < num_present {
            info!(
                "{} of {} genome loci had no neighbours meeting the selection criteria and were not evaluated.",
                num_present - lower.votes.len(),
                num_present
            );
        }

        let mut masked: BTreeSet<(usize, usize)> = lower.votes.iter().filter(|(_, v)| v.is_masked()).map(|(k, _)| *k).collect();
        if let Some(upper) = &upper {
            masked.extend(upper.votes.iter().filter(|(_, v)| v.is_masked()).map(|(k, _)| *k));
        }
        info!(
            "Masked {} of {} genome loci ({} with low identity, {} with high identity).",
            masked.len(),
            num_present,
            lower.num_masked(),
            upper.as_ref().map_or_else(|| "not tested".to_string(), |u| u.num_masked().to_string())
        );

        // write results to output files
        let out_dir = Path::new(&self.out_dir);
        self.write_comparisons(&locus_set, &lower, &out_dir.join("divergent_loci_comparisons.tsv"))?;
        if let Some(upper) = &upper {
            self.write_comparisons(&locus_set, upper, &out_dir.join("divergent_loci_upper_comparisons.tsv"))?;
        }
        let matches = match &upper {
            Some(u) => identity_matches(&u.comparisons, Some(&u.outliers), &masked),
            None => identity_matches(&unevaluated_upper, None, &masked),
        };
        self.write_masked_loci(
            &locus_set,
            &lower,
            upper.as_ref(),
            &masked,
            &matches,
            taxonomy.as_ref(),
            &out_dir.join("divergent_loci_masked.tsv"),
        )?;
        self.write_locus_summary(&locus_set, &lower, upper.as_ref(), &masked, &out_dir.join("divergent_loci_summary.tsv"))?;
        let main_sources = upper.as_ref().map(|u| main_sources(&locus_set, u, &matches));
        self.write_genome_summary(
            &locus_set,
            &lower,
            upper.as_ref(),
            &masked,
            main_sources.as_ref(),
            &out_dir.join("divergent_loci_genome_summary.tsv"),
        )?;

        info!("Writing masked MSAs to file.");
        let msa_out_dir = out_dir.join("masked_msas");
        std::fs::create_dir_all(&msa_out_dir)
            .with_context(|| format!("Failed to create directory: {}", msa_out_dir.display()))?;
        self.write_masked_msas(&locus_set, &masked, &msa_out_dir)?;

        if self.plot {
            info!("Creating plots.");
            let plot_dir = out_dir.join("plots");
            std::fs::create_dir_all(&plot_dir)
                .with_context(|| format!("Failed to create directory: {}", plot_dir.display()))?;

            let settings = PlotSettings {
                title: "Divergent loci (low identity)".to_string(),
                x_label: "genome-wide AAI to neighbour excluding locus (%)".to_string(),
                y_label: "locus identity (%)".to_string(),
                window_size: self.lower.window_size,
                outlier_threshold: self.lower.outlier_threshold,
                total_masked: lower.num_masked(),
                total_genome_loci: lower.votes.len(),
            };
            let plot_data = self.plot_data(&locus_set, &lower);
            write_plots(&plot_data, AxisScale::Linear, &settings, &plot_dir.join("divergent_loci_plots"))?;
            write_plots(&plot_data, AxisScale::Log, &settings, &plot_dir.join("divergent_loci_plots_log"))?;

            if let (Some(upper), Some(params)) = (&upper, &self.upper) {
                let settings = PlotSettings {
                    title: "Divergent loci (high identity)".to_string(),
                    x_label: "locus identity to neighbour (%)".to_string(),
                    y_label: "genome-wide AAI excluding locus (%)".to_string(),
                    window_size: params.window_size,
                    outlier_threshold: params.outlier_threshold,
                    total_masked: upper.num_masked(),
                    total_genome_loci: upper.votes.len(),
                };
                let plot_data = self.plot_data(&locus_set, upper);
                write_plots(&plot_data, AxisScale::Linear, &settings, &plot_dir.join("divergent_loci_upper_plots"))?;
                write_plots(&plot_data, AxisScale::Log, &settings, &plot_dir.join("divergent_loci_upper_plots_log"))?;
            }
        }

        Ok(())
    }

    /// Calculate expected values, identify outlier comparisons, and tally outlier votes for an outlier test.
    fn evaluate(
        &self,
        locus_set: &LocusSet,
        mut comparisons: Vec<LocusComparison>,
        test: OutlierTest,
        params: &TestParams,
    ) -> TestResults {
        comparisons.sort_by(|a, b| {
            (a.locus_idx, a.genome_idx, a.neighbour_rank).cmp(&(b.locus_idx, b.genome_idx, b.neighbour_rank))
        });

        let expected = self.calc_expected(locus_set, &comparisons, test, params.window_size);
        let outliers: Vec<bool> = comparisons
            .iter()
            .zip(&expected)
            .map(|(c, &e)| test.y(c) < e - params.outlier_threshold)
            .collect();

        let votes = tally_votes(
            comparisons.iter().map(|c| (c.genome_idx, c.locus_idx)),
            outliers.iter().copied(),
        );

        let results = TestResults {
            test,
            comparisons,
            expected,
            outliers,
            votes,
        };
        info!(
            "Loci with {}: performed {} comparisons, identified {} outlier comparisons, and masked {} of {} evaluated genome loci.",
            test.description(),
            results.comparisons.len(),
            results.outliers.iter().filter(|&&o| o).count(),
            results.num_masked(),
            results.votes.len()
        );

        results
    }

    /// Calculate the expected value of each comparison from the trend within each locus.
    fn calc_expected(&self, locus_set: &LocusSet, comparisons: &[LocusComparison], test: OutlierTest, window_size: usize) -> Vec<f64> {
        let mut expected = vec![f64::NAN; comparisons.len()];

        let mut start = 0;
        for locus_idx in 0..locus_set.loci.len() {
            // comparisons are sorted by locus so each locus is a contiguous block
            let end = start + comparisons[start..].iter().take_while(|c| c.locus_idx == locus_idx).count();

            let mut order: Vec<usize> = (start..end).collect();
            order.sort_by(|&a, &b| {
                test.x(&comparisons[a])
                    .partial_cmp(&test.x(&comparisons[b]))
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| a.cmp(&b))
            });

            let xs: Vec<f64> = order.iter().map(|&i| test.x(&comparisons[i])).collect();
            let ys: Vec<f64> = order.iter().map(|&i| test.y(&comparisons[i])).collect();
            for (&i, e) in order.iter().zip(expected_values(self.method, &xs, &ys, window_size)) {
                expected[i] = e;
            }

            start = end;
        }

        expected
    }

    /// Collect data for plotting each locus, ordered by decreasing number of masked genomes.
    fn plot_data(&self, locus_set: &LocusSet, results: &TestResults) -> Vec<LocusPlotData> {
        let mut plot_data: Vec<LocusPlotData> = locus_set
            .loci
            .iter()
            .enumerate()
            .map(|(locus_idx, locus)| LocusPlotData {
                name: locus.name.clone(),
                points: Vec::new(),
                trend: Vec::new(),
                num_genomes: locus_set.present.column(locus_idx).iter().filter(|&&p| p).count(),
                num_masked: 0,
                num_outliers: 0,
            })
            .collect();

        let test = results.test;
        for ((c, &e), &is_outlier) in results.comparisons.iter().zip(&results.expected).zip(&results.outliers) {
            let class = if !is_outlier {
                PointClass::Normal
            } else if results.is_masked(c.genome_idx, c.locus_idx) {
                PointClass::Masked
            } else {
                PointClass::Outlier
            };

            let data = &mut plot_data[c.locus_idx];
            data.points.push((test.x(c), test.y(c), class));
            data.trend.push((test.x(c), e));
            data.num_outliers += is_outlier as usize;
        }

        for ((_, locus_idx), v) in &results.votes {
            plot_data[*locus_idx].num_masked += v.is_masked() as usize;
        }

        for data in &mut plot_data {
            data.trend.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
        }

        plot_data.sort_by(|a, b| {
            b.num_masked
                .cmp(&a.num_masked)
                .then_with(|| b.num_outliers.cmp(&a.num_outliers))
                .then_with(|| a.name.cmp(&b.name))
        });

        plot_data
    }

    /// Write all locus comparisons for an outlier test.
    fn write_comparisons(&self, locus_set: &LocusSet, results: &TestResults, out_file: &Path) -> Result<()> {
        let file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file.display()))?;
        let mut writer = BufWriter::new(file);

        match results.test {
            OutlierTest::LowIdentity => writeln!(
                writer,
                "locus\tgenome_id\tneighbour_id\tneighbour_rank\tgenome_aai_loo\tlocus_identity\tlocus_aligned\texpected_identity\tidentity_deficit\toutlier"
            )?,
            OutlierTest::HighIdentity => writeln!(
                writer,
                "locus\tgenome_id\tneighbour_id\tneighbour_rank\tlocus_identity\tlocus_aligned\tgenome_aai_loo\texpected_aai\taai_deficit\toutlier"
            )?,
        }

        for ((c, e), is_outlier) in results.comparisons.iter().zip(&results.expected).zip(&results.outliers) {
            let (first, second) = match results.test {
                OutlierTest::LowIdentity => (format!("{:.4}\t{:.4}\t{}", c.similarity, c.locus_identity, c.locus_aligned), e - c.locus_identity),
                OutlierTest::HighIdentity => (format!("{:.4}\t{}\t{:.4}", c.locus_identity, c.locus_aligned, c.similarity), e - c.similarity),
            };
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}",
                locus_set.loci[c.locus_idx].name,
                locus_set.genome_ids[c.genome_idx],
                locus_set.genome_ids[c.neighbour_idx],
                c.neighbour_rank,
                first,
                e,
                second,
                is_outlier
            )?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write genome loci that were masked, along with the outlier votes from each test and the genomes whose
    /// copies of the locus are most similar.
    #[allow(clippy::too_many_arguments)]
    fn write_masked_loci(
        &self,
        locus_set: &LocusSet,
        lower: &TestResults,
        upper: Option<&TestResults>,
        masked: &BTreeSet<(usize, usize)>,
        matches: &IdentityMatches,
        taxonomy: Option<&Taxonomy>,
        out_file: &Path,
    ) -> Result<()> {
        let file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file.display()))?;
        let mut writer = BufWriter::new(file);

        let mut keys: Vec<&(usize, usize)> = masked.iter().collect();
        keys.sort_by_key(|(genome_idx, locus_idx)| (*locus_idx, *genome_idx));

        let votes = |results: Option<&TestResults>, key: &(usize, usize)| -> (String, String) {
            match results {
                Some(r) => {
                    let v = r.votes.get(key).copied().unwrap_or_default();
                    (v.comparisons.to_string(), v.outliers.to_string())
                }
                None => ("NA".to_string(), "NA".to_string()),
            }
        };

        writeln!(
            writer,
            "locus\tgenome_id\tnum_comparisons\tnum_outliers\tnum_upper_comparisons\tnum_upper_outliers\tmasked_by\t\
             best_match_id\tbest_match_identity\tbest_match_aai\tsource_coherent\tmatch_taxon\tshared_taxon"
        )?;
        for key in keys {
            let (genome_idx, locus_idx) = *key;
            let by_lower = lower.is_masked(genome_idx, locus_idx);
            let by_upper = upper.is_some_and(|u| u.is_masked(genome_idx, locus_idx));
            let masked_by = match (by_lower, by_upper) {
                (true, true) => "both",
                (true, false) => "low_identity",
                _ => "high_identity",
            };

            let (lower_cmps, lower_outliers) = votes(Some(lower), key);
            let (upper_cmps, upper_outliers) = votes(upper, key);

            // genomes whose copies of the locus are most similar: the neighbours responsible for a high-identity
            // mask, otherwise all neighbours by locus identity
            let ranked = matches.get(key).map(Vec::as_slice).unwrap_or_default();
            let sources: Vec<usize> = if by_upper {
                ranked.iter().filter(|(_, is_outlier)| *is_outlier).map(|(c, _)| c.neighbour_idx).collect()
            } else {
                ranked.iter().map(|(c, _)| c.neighbour_idx).collect()
            };

            let best_match = match ranked.first() {
                Some((c, _)) => format!("{}\t{:.4}\t{:.4}", locus_set.genome_ids[c.neighbour_idx], c.locus_identity, c.similarity),
                None => "NA\tNA\tNA".to_string(),
            };
            let source_coherent = if by_upper && sources.len() > 1 {
                let coherent = sources.iter().enumerate().all(|(i, &a)| {
                    sources[i + 1..]
                        .iter()
                        .all(|&b| genome_aai(locus_set, a, b).is_some_and(|aai| aai >= SOURCE_COHERENT_AAI))
                });
                coherent.to_string()
            } else {
                "NA".to_string()
            };
            let ids = |genomes: &[usize]| genomes.iter().map(|&g| locus_set.genome_ids[g].as_str()).collect::<Vec<_>>();
            let lowest_common_taxon = |genomes: &[usize]| {
                taxonomy
                    .filter(|_| !sources.is_empty())
                    .and_then(|t| t.lowest_common_taxon(ids(genomes)))
                    .unwrap_or("NA")
                    .to_string()
            };
            let match_taxon = lowest_common_taxon(&sources);
            let shared_taxon = lowest_common_taxon(&[&[genome_idx][..], &sources].concat());

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                locus_set.loci[locus_idx].name,
                locus_set.genome_ids[genome_idx],
                lower_cmps,
                lower_outliers,
                upper_cmps,
                upper_outliers,
                masked_by,
                best_match,
                source_coherent,
                match_taxon,
                shared_taxon
            )?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write summary statistics for each locus.
    fn write_locus_summary(
        &self,
        locus_set: &LocusSet,
        lower: &TestResults,
        upper: Option<&TestResults>,
        masked: &BTreeSet<(usize, usize)>,
        out_file: &Path,
    ) -> Result<()> {
        let file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file.display()))?;
        let mut writer = BufWriter::new(file);

        // (genomes with comparisons, comparisons, outlier comparisons, masked genomes) for a test
        let stats = |results: &TestResults| {
            let mut stats = vec![(0usize, 0usize, 0usize, 0usize); locus_set.loci.len()];
            for ((_, locus_idx), v) in &results.votes {
                let s = &mut stats[*locus_idx];
                s.0 += 1;
                s.1 += v.comparisons;
                s.2 += v.outliers;
                s.3 += v.is_masked() as usize;
            }
            stats
        };
        let lower_stats = stats(lower);
        let upper_stats = upper.map(stats);

        let mut num_masked = vec![0usize; locus_set.loci.len()];
        for (_, locus_idx) in masked {
            num_masked[*locus_idx] += 1;
        }

        writeln!(
            writer,
            "locus\tnum_genomes\tnum_genomes_compared\tnum_comparisons\tnum_outlier_comparisons\tnum_masked_low_identity\t\
             num_upper_comparisons\tnum_upper_outlier_comparisons\tnum_masked_high_identity\tnum_masked\tmasked_perc"
        )?;
        for (locus_idx, locus) in locus_set.loci.iter().enumerate() {
            let num_genomes = locus_set.present.column(locus_idx).iter().filter(|&&p| p).count();
            let (num_compared, num_comparisons, num_outliers, num_masked_lower) = lower_stats[locus_idx];
            let upper_cols = match &upper_stats {
                Some(s) => format!("{}\t{}\t{}", s[locus_idx].1, s[locus_idx].2, s[locus_idx].3),
                None => "NA\tNA\tNA".to_string(),
            };
            let masked_perc = if num_genomes > 0 {
                100.0 * num_masked[locus_idx] as f64 / num_genomes as f64
            } else {
                0.0
            };

            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}",
                locus.name,
                num_genomes,
                num_compared,
                num_comparisons,
                num_outliers,
                num_masked_lower,
                upper_cols,
                num_masked[locus_idx],
                masked_perc
            )?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write the number and percentage of loci present in each genome before and after masking, if the genome
    /// is suspect, and the main source of its high-identity masks.
    fn write_genome_summary(
        &self,
        locus_set: &LocusSet,
        lower: &TestResults,
        upper: Option<&TestResults>,
        masked: &BTreeSet<(usize, usize)>,
        main_sources: Option<&HashMap<usize, (usize, usize)>>,
        out_file: &Path,
    ) -> Result<()> {
        let file = File::create(out_file).with_context(|| format!("Failed to create file: {}", out_file.display()))?;
        let mut writer = BufWriter::new(file);

        let num_genomes = locus_set.num_genomes();
        let count = |keys: &mut dyn Iterator<Item = usize>| {
            let mut counts = vec![0usize; num_genomes];
            for genome_idx in keys {
                counts[genome_idx] += 1;
            }
            counts
        };
        let masked_lower = count(&mut lower.votes.iter().filter(|(_, v)| v.is_masked()).map(|((g, _), _)| *g));
        let masked_upper = upper.map(|u| count(&mut u.votes.iter().filter(|(_, v)| v.is_masked()).map(|((g, _), _)| *g)));
        let masked_total = count(&mut masked.iter().map(|(g, _)| *g));

        let num_loci = locus_set.loci.len();
        let perc = |n: usize| 100.0 * n as f64 / num_loci as f64;

        writeln!(
            writer,
            "genome_id\tnum_loci\tloci_perc\tnum_masked_low_identity\tnum_masked_high_identity\tnum_masked\t\
             num_loci_after_masking\tloci_perc_after_masking\tsuspect\tmain_source_id\tmain_source_masks"
        )?;
        let mut suspect = Vec::new();
        for (genome_idx, genome_id) in locus_set.genome_ids.iter().enumerate() {
            let num_present = locus_set.present.row(genome_idx).iter().filter(|&&p| p).count();
            let num_after = num_present - masked_total[genome_idx];
            let is_suspect = is_suspect(num_present, masked_total[genome_idx], self.suspect_masked_perc);
            if is_suspect {
                suspect.push(format!("{} ({} of {} loci)", genome_id, masked_total[genome_idx], num_present));
            }
            let main_source = match main_sources {
                Some(sources) => match sources.get(&genome_idx) {
                    Some(&(source_idx, num_masks)) => format!("{}\t{}", locus_set.genome_ids[source_idx], num_masks),
                    None => "NA\t0".to_string(),
                },
                None => "NA\tNA".to_string(),
            };
            writeln!(
                writer,
                "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
                genome_id,
                num_present,
                perc(num_present),
                masked_lower[genome_idx],
                masked_upper.as_ref().map_or_else(|| "NA".to_string(), |m| m[genome_idx].to_string()),
                masked_total[genome_idx],
                num_after,
                perc(num_after),
                is_suspect,
                main_source
            )?;
        }

        writer.flush()?;

        if suspect.is_empty() {
            info!("No genomes have {}% or more of their loci masked.", self.suspect_masked_perc);
        } else {
            info!(
                "{} genome(s) have {}% or more of their loci masked and are flagged as suspect (e.g. contaminated or chimeric): {}",
                suspect.len(),
                self.suspect_masked_perc,
                suspect.join(", ")
            );
        }

        Ok(())
    }

    /// Write MSAs with masked genome loci replaced by gaps.
    fn write_masked_msas(&self, locus_set: &LocusSet, masked: &BTreeSet<(usize, usize)>, msa_out_dir: &Path) -> Result<()> {
        for (locus_idx, locus) in locus_set.loci.iter().enumerate() {
            let input_file = locus.file.to_string_lossy();
            let out_file = msa_out_dir.join(format!("{}_masked.fasta", locus.name));
            let (mut writer, _) = create_writer(&out_file.to_string_lossy(), check_compressed_input(&input_file))?;

            let mut reader = parse_fastx_file(&locus.file)
                .with_context(|| format!("Failed to open alignment file: {}", input_file))?;
            while let Some(record) = reader.next() {
                let rec = record.with_context(|| "Failed to parse FASTA record")?;
                let header = String::from_utf8_lossy(rec.id());
                let genome_id = genome_id_from_header(&header);

                let genome_idx = locus_set
                    .genome_ids
                    .binary_search_by(|id| id.as_str().cmp(genome_id))
                    .map_err(|_| anyhow::anyhow!("Genome {} not found in: {}", genome_id, input_file))?;

                writeln!(writer, ">{}", header)?;
                if masked.contains(&(genome_idx, locus_idx)) {
                    writer.write_all(&vec![b'-'; rec.seq().len()])?;
                } else {
                    writer.write_all(&rec.seq())?;
                }
                writeln!(writer)?;
            }

            writer.flush()?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_main_source() {
        // genomes 1 and 2 are one lineage, 7 another
        let same_lineage = |a: usize, b: usize| (a <= 2) == (b <= 2);
        assert_eq!(main_source(&[2, 1, 2, 7], same_lineage), (2, 3));
        // ties go to the lowest genome index
        assert_eq!(main_source(&[7, 1], same_lineage), (1, 1));
        assert_eq!(main_source(&[7], same_lineage), (7, 1));
    }

    #[test]
    fn test_is_suspect() {
        assert!(is_suspect(30, 3, 10.0));
        assert!(!is_suspect(31, 3, 10.0));
        // too few masked loci, regardless of percentage
        assert!(!is_suspect(4, 2, 10.0));
        assert!(is_suspect(44, 29, 10.0));
    }
}
