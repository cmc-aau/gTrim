use clap::{Parser, Subcommand};
use std::path::Path;

const DEFAULT_GAP_THRESHOLD: f64 = 100.0;

const DEFAULT_MAX_RESIDUE_PRUNED: usize = 100000;
const DEFAULT_NUM_THREADS: usize = 1;

const DEFAULT_MSA_LENGTH: usize = 5000;
const DEFAULT_MIN_PERC_AA: f64 = 50.0;
const DEFAULT_MIN_CONSENSUS: f64 = 25.0;
const DEFAULT_MAX_CONSENSUS: f64 = 95.0;
const DEFAULT_SEED: u64 = 12345;

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

    /// Apply mask to multiple sequence alignment.
    Mask(MaskArgs),

    /// Calculate column statistics.
    ColumnStats(ColumnStatsArgs),

    /// Filter gappy, highly conserved, or poorly conserved columns and subsample columns from each marker.
    Trim(TrimArgs),
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
pub struct TrimArgs {
    /// File containing the concatenated marker gene alignment.
    #[arg(short = 'i', long, help_heading = "Required", value_parser = validate_input_file, required = true)]
    pub alignment_file: String,

    /// Tab-separated file listing the markers in the order they are concatenated in the alignment.
    /// Requires the columns 'Marker Id', 'Description' and 'Length (bp)' (e.g. gtdb_r<release>_bac120_markers_info.tsv).
    #[arg(short = 'm', long, help_heading = "Required", value_parser = validate_file, required = true)]
    pub marker_info_file: String,

    /// Output directory.
    #[arg(short = 'o', long, help_heading = "Required", required = true)]
    pub out_dir: String,

    /// Desired length of the alignment after randomly subsampling the same number of columns from each marker.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MSA_LENGTH, value_parser = validate_msa_length)]
    pub msa_length: usize,

    /// Minimum percentage of sequences with a residue (not a gap) required to retain a column.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MIN_PERC_AA, value_parser = validate_perc)]
    pub min_perc_aa: f64,

    /// Minimum percentage of residues that must be the most common residue to retain a column (inclusive).
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MIN_CONSENSUS, value_parser = validate_perc)]
    pub min_consensus: f64,

    /// Maximum percentage of residues that can be the most common residue to retain a column (exclusive).
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_MAX_CONSENSUS, value_parser = validate_perc)]
    pub max_consensus: f64,

    /// Seed for randomly subsampling columns.
    #[arg(long, help_heading = "Optional", default_value_t = DEFAULT_SEED)]
    pub seed: u64,

    /// Number of threads.
    #[arg(short = 't', long, help_heading = "Optional", default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads: usize,
}

fn validate_msa_length(value: &str) -> Result<usize, String> {
    let len: usize = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if len == 0 {
        return Err("MSA length must be greater than 0".to_string());
    }

    Ok(len)
}

fn validate_perc(value: &str) -> Result<f64, String> {
    let perc: f64 = value
        .parse()
        .map_err(|_| format!("`{value}` isn't a valid value"))?;

    if !(0.0..=100.0).contains(&perc) {
        return Err("Percentage must be in the range [0, 100]".to_string());
    }

    Ok(perc)
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
