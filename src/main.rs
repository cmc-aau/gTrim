use std::env;
use std::path::Path;
use std::time::Instant;

use log::{error, info};
use clap::Parser;
use anyhow::Result;

use gtrim::logging;
use gtrim::cli;
use gtrim::prune::gap::gap_pruner;
use gtrim::prune::mask::mask_pruner;
use gtrim::stats::column_stats;

/// Common initialization required by all commands.
fn init(out_dir: &String, threads: usize) -> Result<()> {
    // create output directory if it doesn't already exist
    std::fs::create_dir_all(out_dir)
        .map_err(|e| anyhow::anyhow!("Failed to create output directory '{}': {}", out_dir, e))?;

    // setup logger
    if let Err(e) = logging::setup_logger(Path::new(out_dir)) {
        eprintln!("Failed to setup logger: {}", e);
        std::process::exit(1);
    }

    // report version and echo executed command
    const VERSION: &str = env!("CARGO_PKG_VERSION");
    info!("{} v{}", env!("CARGO_PKG_NAME"), VERSION);
    info!("{}", env::args().collect::<Vec<String>>().join(" "));

    // setup threadpool for rayon
    info!("Using {} threads.", threads);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;

    Ok(())
}

/// Prune columns based on gap threshold
fn prune_gap(args: &cli::PruneGapArgs) -> Result<()> {
    init(&args.out_dir, args.threads)?;

    let mut pruner = gap_pruner::GapPruner::new(
        args.alignment_file.clone(),
        args.gap_threshold,
        args.max_residue_pruned,
        args.out_dir.clone(),
    );

    if let Err(e) = pruner.run() {
        error!("Error: {}", e);
        std::process::exit(1);
    }

    Ok(())
}

/// Apply mask to multiple sequence alignment.
fn mask(args: &cli::MaskArgs) -> Result<()> {
    init(&args.out_dir, 1)?;

    let mut pruner = mask_pruner::MaskPruner::new(
        args.alignment_file.clone(),
        args.mask_file.clone(),
        args.out_dir.clone(),
    );

    if let Err(e) = pruner.run() {
        error!("Error: {}", e);
        std::process::exit(1);
    }

    Ok(())
} 

/// Calculate column statistics.
fn column_stats(args: &cli::ColumnStatsArgs) -> Result<()> {
    init(&args.out_dir, 1)?;

    let mut stats = column_stats::ColumnStats::new(
        args.alignment_file.clone(),
        args.out_dir.clone(),
    );

    if let Err(e) = stats.run() {
        error!("Error: {}", e);
        std::process::exit(1);
    }

    Ok(())
}

fn main() -> Result<()> {
    let start_time = Instant::now();

    let cli = cli::Cli::parse();

    match &cli.command {
        cli::Commands::PruneGap(args) => {
            prune_gap(args)?;
        }
        cli::Commands::Mask(args) => {
            mask(args)?;
        }
        cli::Commands::ColumnStats(args) => {
            column_stats(args)?;
        }
    }
    
    let elapsed = start_time.elapsed();
    info!("Execution time: {:.2} seconds", elapsed.as_secs_f64());

    Ok(())
}
