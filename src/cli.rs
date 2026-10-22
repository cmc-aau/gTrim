use clap::{Parser, Subcommand};
use std::path::Path;

use crate::prune::divergent::outliers::OutlierMethod;

const DEFAULT_GAP_THRESHOLD: f64 = 100.0;

const DEFAULT_MAX_RESIDUE_PRUNED: usize = 100000;
const DEFAULT_NUM_PERMUTATIONS: usize = 100;
const DEFAULT_RESIDUES_TO_PRUNE_PERC: f64 = 1.0;
const DEFAULT_NUM_THREADS: usize = 1;

const DEFAULT_NUM_NEIGHBOURS: usize = 3;
const DEFAULT_WINDOW_SIZE: usize = 120;
const DEFAULT_OUTLIER_THRESHOLD: f64 = 25.0;
const DEFAULT_MIN_ALIGNED_FRAC: f64 = 0.5;
const DEFAULT_MIN_ALIGNED_COLUMNS: u32 = 50;
const DEFAULT_SUSPECT_MASKED_PERC: f64 = 10.0;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(disable_help_subcommand = true)]
#[command(arg_required_else_help = true)]
#[command(styles=get_styles())]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Prune columns based on gap threshold
    PruneGap(PruneGapArgs),

    /// Prune columns based on chi-square test.
    PruneChi2(PruneChi2Args),

    /// Apply mask to multiple sequence alignment.
    Mask(MaskArgs),

    /// Calculate column statistics.
    ColumnStats(ColumnStatsArgs),

    /// Mask loci that are highly divergent relative to closest neighbouring genomes.
    MaskDivergentLoci(MaskDivergentLociArgs),
}

#[derive(Parser)]
pub struct PruneGapArgs {
    /// File containing the multiple sequence alignment.
    #[arg(short = 'i', long, help_heading = "Required", value_parser = validate_input_file, required = true)]
    pub alignment_file: String,

    /// Output directory.
    #[arg(short = 'o', long, help_heading = "Required", required = true)]
    pub out_dir: String,

    /// Maximum percentage of gaps to retain column.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_GAP_THRESHOLD, value_parser = validate_gap_threshold)]
    pub gap_threshold: f64,

    /// Maximum number of columns to prune.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MAX_RESIDUE_PRUNED, value_parser = validate_max_residue_pruned)]
    pub max_residue_pruned: usize,

    /// Number of threads to use for pruning columns.
    #[arg(short = 't', long, help_heading = "Optional", default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct PruneChi2Args {
    /// File containing the multiple sequence alignment.
    #[arg(short = 'i', long, help_heading = "Required", value_parser = validate_input_file, required = true)]
    pub alignment_file: String,

    /// Output directory.
    #[arg(short = 'o', long, help_heading = "Required", required = true)]
    pub out_dir: String,

    /// Maximum number of columns to prune.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MAX_RESIDUE_PRUNED, value_parser = validate_max_residue_pruned)]
    pub max_residue_pruned: usize,

    /// Percentage of top delta-chi2 columns to remove per iteration [0, 100]. Set to 0 to prune 1 column per iteration.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_RESIDUES_TO_PRUNE_PERC, value_parser = validate_residues_to_prune_perc)]
    pub residues_to_prune_perc: f64,

    /// Number of alignment permutations.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_NUM_PERMUTATIONS, value_parser = validate_permutations)]
    pub permutations: usize,

    /// Prunes columns until no significant reduction in compositional bias is observed.
    /// Recommend setting residues_to_prune_perc to 1-2%. The number of columns to prune 
    /// will be reduced to 0.1% of the original alignment length per iteration after a 
    /// threshold of 99.0% is reached.
    #[arg(long, help_heading = "Optional", conflicts_with = "rapid_convergence")]
    pub touchdown: bool,

    /// Prunes columns until no significant reduction in compositional bias is observed.
    /// Recommend setting residues_to_prune_perc to 1-2%. The number of columns pruned
    /// is always residues_to_prune_perc multipled by the length of the remaining 
    /// alignment. This will generally converge faster than touchdown mode.
    #[arg(long, help_heading = "Optional", conflicts_with = "touchdown")]
    pub rapid_convergence: bool,

    /// Number of threads to use when calculating permutations. Threads replicate the alignment in memory,
    /// which results in higher peak memory requirements.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads_permutations: usize,

    /// Number of threads to use for pruning columns.
    #[arg(short = 't', long, help_heading = "Optional", default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads_pruning: usize,
}

#[derive(Parser)]
pub struct MaskArgs {
    /// File containing the multiple sequence alignment.
    #[arg(short = 'i', long, help_heading = "Required", value_parser = validate_input_file, required = true)]
    pub alignment_file: String,

    /// File containing the binary mask.
    #[arg(short = 'm', long, help_heading = "Required", value_parser = validate_file, required = true)]
    pub mask_file: String,

    /// Output directory.
    #[arg(short = 'o', long, help_heading = "Required", required = true)]
    pub out_dir: String,

    /// Number of threads to use for pruning columns.
    #[arg(short = 't', long, help_heading = "Optional", default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct ColumnStatsArgs {
    /// File containing the multiple sequence alignment.
    #[arg(short = 'i', long, help_heading = "Required", value_parser = validate_input_file, required = true)]
    pub alignment_file: String,

    /// Output directory.
    #[arg(short = 'o', long, help_heading = "Required", required = true)]
    pub out_dir: String,

    /// Number of threads to use for pruning columns.
    #[arg(short = 't', long, help_heading = "Optional", default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads: usize,
}

#[derive(Parser)]
pub struct MaskDivergentLociArgs {
    /// Directory containing one MSA per locus (.faa, .fna, .fasta, optionally .gz). FASTA headers must start with the genome ID.
    #[arg(short = 'i', long, help_heading = "Required", value_parser = validate_dir, required = true)]
    pub msa_dir: String,

    /// Output directory.
    #[arg(short = 'o', long, help_heading = "Required", required = true)]
    pub out_dir: String,

    /// Low-identity test: number of closest neighbouring genomes (by AAI excluding the evaluated locus) to compare each locus against.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_NUM_NEIGHBOURS, value_parser = validate_neighbours)]
    pub neighbours: usize,

    /// Method used to calculate the expected value of each comparison, for both tests.
    #[arg(long, help_heading = "Optional", value_enum, default_value_t = OutlierMethod::MonotoneMedian)]
    pub method: OutlierMethod,

    /// Low-identity test: number of observations nearest in genome-wide AAI used to calculate the running median of locus identity.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_WINDOW_SIZE, value_parser = validate_window_size)]
    pub window_size: usize,

    /// Low-identity test: percentage points below the expected locus identity for a comparison to be an outlier.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_OUTLIER_THRESHOLD, value_parser = validate_outlier_threshold)]
    pub outlier_threshold: f64,

    /// Minimum fraction of the shorter sequence that must be aligned for a comparison to be considered.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MIN_ALIGNED_FRAC, value_parser = validate_min_aligned_frac)]
    pub min_aligned_frac: f64,

    /// Minimum number of aligned columns at a locus for a comparison to be considered (0 to disable).
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MIN_ALIGNED_COLUMNS)]
    pub min_aligned_columns: u32,

    /// High-identity test: number of closest neighbouring genomes (by identity of the evaluated locus) to compare each locus against.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_NUM_NEIGHBOURS, value_parser = validate_neighbours)]
    pub upper_neighbours: usize,

    /// High-identity test: number of observations nearest in locus identity used to calculate the running median of genome-wide AAI.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_WINDOW_SIZE, value_parser = validate_window_size)]
    pub upper_window_size: usize,

    /// High-identity test: percentage points below the expected genome-wide AAI for a comparison to be an outlier.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_OUTLIER_THRESHOLD, value_parser = validate_outlier_threshold)]
    pub upper_outlier_threshold: f64,

    /// Do not mask loci with unexpectedly high identity, i.e. loci whose closest matches by locus identity
    /// have unexpectedly low genome-wide AAI to the genome (e.g. contamination or lateral transfer).
    #[arg(long = "no-upper-mask", alias = "no_upper_mask", help_heading = "Optional")]
    pub no_upper_mask: bool,

    /// Tab-separated file of genome ID and taxonomy, with ranks separated by ';' (e.g. GTDB taxonomy file). Used to
    /// report the lowest common taxon of the genomes whose copies of each masked locus are most similar.
    #[arg(long, help_heading = "Optional", value_parser = validate_file)]
    pub taxonomy: Option<String>,

    /// Flag a genome as suspect (e.g. contaminated or chimeric) in the genome summary if at least this percentage
    /// of its loci, and at least 3 loci, are masked.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_SUSPECT_MASKED_PERC, value_parser = validate_suspect_masked_perc)]
    pub suspect_masked_perc: f64,

    /// Create SVG and PNG scatter plots of locus identity against genome-wide AAI for each locus (linear and log scale).
    #[arg(long, help_heading = "Optional")]
    pub plot: bool,

    /// Number of threads.
    #[arg(short = 't', long, help_heading = "Optional", default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads: usize,
}

fn validate_neighbours(value: &str) -> Result<usize, String> {
    let neighbours: usize = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if neighbours == 0 {
        return Err("Number of neighbours must be greater than 0".to_string());
    }

    Ok(neighbours)
}

fn validate_window_size(value: &str) -> Result<usize, String> {
    let window: usize = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if window == 0 {
        return Err("Window size must be greater than 0".to_string());
    }

    Ok(window)
}

fn validate_outlier_threshold(value: &str) -> Result<f64, String> {
    let threshold: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if !(0.0..=100.0).contains(&threshold) {
        return Err("Outlier threshold must be in the range [0, 100]".to_string());
    }

    Ok(threshold)
}

fn validate_min_aligned_frac(value: &str) -> Result<f64, String> {
    let frac: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if !(0.0..=1.0).contains(&frac) {
        return Err("Minimum aligned fraction must be in the range [0, 1]".to_string());
    }

    Ok(frac)
}

fn validate_suspect_masked_perc(value: &str) -> Result<f64, String> {
    let perc: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if !(0.0..=100.0).contains(&perc) {
        return Err("Suspect masked percentage must be in the range [0, 100]".to_string());
    }

    Ok(perc)
}

fn validate_dir(dir: &str) -> Result<String, String> {
    let path = Path::new(dir);
    if !path.is_dir() {
        return Err(format!("Directory '{}' does not exist", dir));
    }

    Ok(dir.to_string())
}

fn validate_gap_threshold(value: &str) -> Result<f64, String> {
    let perc: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if !(0.0..=100.0).contains(&perc) {
        return Err("Gap threshold must be in the range [0, 100]".to_string());
    }

    Ok(perc)
}

fn validate_threads(threads: &str) -> Result<usize, String> {
    let threads: usize = threads
        .parse()
        .map_err(|_| format!("`{threads}` isn't a valid value"))?;

    if !(1..=1024).contains(&threads) {
        return Err("Threads must be in the range [1, 1024]".to_string());
    }

    Ok(threads)
}

fn validate_max_residue_pruned(value: &str) -> Result<usize, String> {
    let max_residue: usize = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if max_residue == 0 {
        return Err("Maximum residues to prune must be greater than 0".to_string());
    }

    Ok(max_residue)
}

fn validate_residues_to_prune_perc(value: &str) -> Result<f64, String> {
    let perc: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if !(0.0..=100.0).contains(&perc) {
        return Err("Residues to prune percentage must be in the range [0, 100]".to_string());
    }

    Ok(perc)
}

fn validate_permutations(value: &str) -> Result<usize, String> {
    let perms: usize = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if perms < 5 {
        return Err("Number of permutations should be at least 5".to_string());
    }

    Ok(perms)
}

fn validate_input_file(file: &str) -> Result<String, String> {
    let path = Path::new(file);
    if !path.exists() {
        return Err(format!("File '{}' does not exist", file));
    }
    if !path.is_file() {
        return Err(format!("'{}' is not a file", file));
    }
    if path.extension().and_then(|s| s.to_str()) != Some("faa") 
        && path.extension().and_then(|s| s.to_str()) != Some("fna") 
        && path.extension().and_then(|s| s.to_str()) != Some("fasta")
        && path.extension().and_then(|s| s.to_str()) != Some("gz") {
        return Err("Input file must have .faa, .fna, .fasta, or .gz extension".to_string());
    }
    Ok(file.to_string())
}

fn validate_file(file: &str) -> Result<String, String> {
    let path = Path::new(file);
    if !path.exists() {
        return Err(format!("File '{}' does not exist", file));
    }
    if !path.is_file() {
        return Err(format!("'{}' is not a file", file));
    }

    Ok(file.to_string())
}

fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}

#[test]
fn test_verify_cli() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}
