use clap::{Parser, Subcommand};
use std::path::Path;


const DEFAULT_MAX_RESIDUE_PRUNED: usize = 100000;
const DEFAULT_NUM_PERMUTATIONS: usize = 100;
const DEFAULT_RESIDUES_TO_PRUNE_PERC: f64 = 1.0;
const DEFAULT_NUM_THREADS: usize = 1;

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
    /// Prune alignment based on chi-square test.
    Chi2Prune(Chi2PruneArgs),
}

#[derive(Parser)]
pub struct Chi2PruneArgs {
    /// File containing the multiple sequence alignment.
    #[arg(long, value_parser = validate_input_file, required = true)]
    pub file: String,

    /// Output directory.If not specified, output files will be written 
    /// to the same directory as the input file.
    #[arg(long)]
    pub out_dir: Option<String>,

    /// Maximum number of columns to prune.
    #[arg(long, default_value_t = DEFAULT_MAX_RESIDUE_PRUNED, value_parser = validate_max_residue_pruned)]
    pub max_residue_pruned: usize,

    /// Percentage of top delta-chi2 columns to remove per iteration [0, 100]. Set to 0 to prune 1 column per iteration.
    #[arg(long, default_value_t = DEFAULT_RESIDUES_TO_PRUNE_PERC, value_parser = validate_residues_to_prune_perc)]
    pub residues_to_prune_perc: f64,

    /// Number of alignment permutations.
    #[arg(long, default_value_t = DEFAULT_NUM_PERMUTATIONS, value_parser = validate_permutations)]
    pub permutations: usize,

    /// Prunes columns until no significant reduction in compositional bias is observed.
    /// Recommend setting residues_to_prune_perc to 1-2%. The number of columns to prune 
    /// will be reduced to 0.1% of the original alignment length per iteration after a 
    /// threshold of 99.0% is reached.
    #[arg(long, conflicts_with = "rapid_convergence")]
    pub touchdown: bool,

    /// Prunes columns until no significant reduction in compositional bias is observed.
    /// Recommend setting residues_to_prune_perc to 1-2%. The number of columns pruned
    /// is always residues_to_prune_perc multipled by the length of the remaining 
    /// alignment. This will generally converge faster than touchdown mode.
    #[arg(long, conflicts_with = "touchdown")]
    pub rapid_convergence: bool,

    /// Number of threads to use when calculating permutations. Threads replicate the alignment in memory,
    /// which results in higher peak memory requirements.
    #[arg(long, default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads_permutations: usize,

    /// Number of threads to use for pruning columns.
    #[arg(long, default_value_t = DEFAULT_NUM_THREADS, value_parser = validate_threads)]
    pub threads_pruning: usize,
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
