use std::env;
use std::path::Path;
use std::time::Instant;

use log::{error, info};
use clap::Parser;
use anyhow::Result;

use gtrim::logging;
use gtrim::cli;
use gtrim::prune::gap::gap_pruner;
use gtrim::prune::chi2::chi_square_pruner;
use gtrim::prune::mask::mask_pruner;
use gtrim::prune::divergent::aai::OverlapFilter;
use gtrim::prune::divergent::divergent_loci_masker::{self, TestParams};
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

/// Prune columns based on chi-square test.
fn prune_chi2(args: &cli::PruneChi2Args) -> Result<()> {
    init(&args.out_dir, args.threads_pruning)?;

    info!("Maximum residues to prune: {}", args.max_residue_pruned);
    info!("Permutations: {}", args.permutations);
    info!("Percentage of residues to prune per iteration: {}", args.residues_to_prune_perc);
    info!("Permutation threads: {}", args.threads_permutations);
    info!("Pruning threads: {}", args.threads_pruning);

    let mut pruner = chi_square_pruner::ChiSquarePruner::new(
        args.alignment_file.clone(),
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

/// Mask loci that are highly divergent relative to closest neighbouring genomes.
fn mask_divergent_loci(args: &cli::MaskDivergentLociArgs) -> Result<()> {
    init(&args.out_dir, args.threads)?;

    let lower = TestParams {
        neighbours: args.neighbours,
        window_size: args.window_size,
        outlier_threshold: args.outlier_threshold,
    };
    let upper = (!args.no_upper_mask).then_some(TestParams {
        neighbours: args.upper_neighbours,
        window_size: args.upper_window_size,
        outlier_threshold: args.upper_outlier_threshold,
    });

    info!("Outlier method: {:?}", args.method);
    info!("Minimum aligned fraction: {}", args.min_aligned_frac);
    info!("Minimum aligned columns: {}", args.min_aligned_columns);
    info!("Suspect genome threshold: {}% of loci masked (minimum {} loci)", args.suspect_masked_perc, divergent_loci_masker::MIN_SUSPECT_MASKED);
    info!(
        "Low-identity test: {} closest neighbour(s), running median window size {}, outlier threshold {}",
        lower.neighbours, lower.window_size, lower.outlier_threshold
    );
    match &upper {
        Some(p) => info!(
            "High-identity test: {} closest neighbour(s), running median window size {}, outlier threshold {}",
            p.neighbours, p.window_size, p.outlier_threshold
        ),
        None => info!("High-identity test: disabled"),
    }

    let mut masker = divergent_loci_masker::DivergentLociMasker::new(
        args.msa_dir.clone(),
        args.method,
        lower,
        upper,
        OverlapFilter {
            min_aligned_frac: args.min_aligned_frac,
            min_aligned_columns: args.min_aligned_columns,
        },
        args.suspect_masked_perc,
        args.taxonomy.clone(),
        args.plot,
        args.out_dir.clone(),
    );

    if let Err(e) = masker.run() {
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
        cli::Commands::PruneChi2(args) => {
            prune_chi2(args)?;
        } 
        cli::Commands::Mask(args) => {
            mask(args)?;
        }
        cli::Commands::ColumnStats(args) => {
            column_stats(args)?;
        }
        cli::Commands::MaskDivergentLoci(args) => {
            mask_divergent_loci(args)?;
        }
    }
    
    let elapsed = start_time.elapsed();
    info!("Execution time: {:.2} seconds", elapsed.as_secs_f64());

    Ok(())
}
