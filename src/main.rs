use std::env;
use std::path::Path;
use std::time::Instant;

use log::{error, info};
use clap::Parser;
use anyhow::Result;

use gtrim::logging;
use gtrim::cli::{Cli, Commands, PruneGapArgs, PruneChi2Args};
use gtrim::prune::gap::gap_pruner;
use gtrim::prune::chi2::chi_square_pruner;

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
fn prune_gap(args: &PruneGapArgs) -> Result<()> {
    init(&args.out_dir, args.threads)?;

    let mut pruner = gap_pruner::GapPruner::new(
        args.file.clone(),
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

/// Prune columns based on chi-square test.
fn prune_chi2(args: &PruneChi2Args) -> Result<()> {
    init(&args.out_dir, args.threads_pruning)?;

    info!("Maximum residues to prune: {}", args.max_residue_pruned);
    info!("Permutations: {}", args.permutations);
    info!("Percentage of residues to prune per iteration: {}", args.residues_to_prune_perc);
    info!("Permutation threads: {}", args.threads_permutations);
    info!("Pruning threads: {}", args.threads_pruning);

    let mut pruner = chi_square_pruner::ChiSquarePruner::new(
        args.file.clone(),
        args.max_residue_pruned,
        args.residues_to_prune_perc,
        args.permutations,
        args.touchdown,
        args.rapid_convergence,
        args.threads_permutations,
        args.out_dir.clone(),
    );

    if let Err(e) = pruner.run() {
        error!("Error: {}", e);
        std::process::exit(1);
    }

    Ok(())
}

fn main() -> Result<()> {
    let start_time = Instant::now();

    let cli = Cli::parse();

    match &cli.command {
        Commands::PruneGap(args) => {
            prune_gap(args)?;
        }
        Commands::PruneChi2(args) => {
            prune_chi2(args)?;
        } 
    }
    
    let elapsed = start_time.elapsed();
    info!("Execution time: {:.2} seconds", elapsed.as_secs_f64());

    Ok(())
}
