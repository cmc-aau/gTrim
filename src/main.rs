use std::env;
use std::path::Path;
use std::time::Instant;

use log::{error, info};
use clap::Parser;
use anyhow::Result;

use gtrim::{alignment_pruner, logging};
use gtrim::cli::{Cli, Commands};


fn main() -> Result<()> {
    let start_time = Instant::now();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Chi2Prune(args) => {
            // Extract directory from file path for logging
            let file_path = Path::new(&args.file);
            let out_dir = file_path.parent().unwrap_or_else(|| Path::new("."));

            // Setup logger
            if let Err(e) = logging::setup_logger(out_dir) {
                eprintln!("Failed to setup logger: {}", e);
                std::process::exit(1);
            }

            info!("WitChi version {}", env!("CARGO_PKG_VERSION"));
            info!("{}", env::args().collect::<Vec<String>>().join(" "));

            info!("Maximum residues to prune: {}", args.max_residue_pruned);
            info!("Permutations: {}", args.permutations);
            info!("Percentage of residues to prune per iteration: {}", args.residues_to_prune_perc);
            info!("Permutation threads: {}", args.threads_permutations);
            info!("Pruning threads: {}", args.threads_pruning);

            // build global thread pool that will be used for all
            // multi-threaded operations, except the permutations
            // required for permutation testing as this parallelization
            // heavily impacts peak memory requirements
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.threads_pruning)
                .build_global()?;

            let mut pruner = alignment_pruner::AlignmentPruner::new(
                args.file.clone(),
                args.max_residue_pruned as usize,
                args.permutations,
                args.residues_to_prune_perc,
                args.touchdown,
                args.threads_permutations,
            );

            if let Err(e) = pruner.run() {
                error!("Error: {}", e);
                std::process::exit(1);
            }
        } 
    }
    
    let elapsed = start_time.elapsed();
    info!("Execution time: {:.2} seconds", elapsed.as_secs_f64());

    Ok(())
}
