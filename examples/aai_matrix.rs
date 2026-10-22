//! Calculate the pairwise AAI matrix between all genomes in a directory of per-locus MSAs.
//!
//! AAI is pooled over all loci: identical / aligned positions, where a position is aligned if both
//! genomes have a valid residue (the same definition used by `gtrim mask-divergent-loci`).
//!
//! Usage:
//!   cargo run --release --example aai_matrix -- <msa_dir> <out_dir> <threads>
//!
//! Writes:
//!   <out_dir>/genome_ids.txt    genome IDs, one per line, in matrix order
//!   <out_dir>/aai_matrix.f32    n x n row-major little-endian f32 matrix of AAI (%), NaN if no aligned positions

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use rayon::prelude::*;

use gtrim::prune::divergent::aai::identity_counts;
use gtrim::prune::divergent::locus_set::LocusSet;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        anyhow::bail!("Usage: aai_matrix <msa_dir> <out_dir> <threads>");
    }
    let (msa_dir, out_dir) = (&args[1], Path::new(&args[2]));
    let threads: usize = args[3].parse().context("Invalid number of threads")?;
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global()?;
    std::fs::create_dir_all(out_dir)?;

    eprintln!("Reading MSAs from: {}", msa_dir);
    let ls = LocusSet::from_dir(msa_dir)?;
    let n = ls.num_genomes();
    eprintln!("Calculating AAI between {} genomes across {} loci.", n, ls.loci.len());

    let row = |i: usize| ls.matrix.row(i).to_slice().expect("matrix is in standard layout");
    let rows: Vec<Vec<f32>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..n)
                .map(|j| {
                    let (identical, aligned) = identity_counts(row(i), row(j));
                    if aligned == 0 {
                        f32::NAN
                    } else {
                        100.0 * identical as f32 / aligned as f32
                    }
                })
                .collect()
        })
        .collect();

    let mut ids = BufWriter::new(File::create(out_dir.join("genome_ids.txt"))?);
    for id in &ls.genome_ids {
        writeln!(ids, "{}", id)?;
    }
    ids.flush()?;

    let mut matrix = BufWriter::new(File::create(out_dir.join("aai_matrix.f32"))?);
    for r in &rows {
        for v in r {
            matrix.write_all(&v.to_le_bytes())?;
        }
    }
    matrix.flush()?;

    eprintln!("Wrote AAI matrix to: {}", out_dir.display());
    Ok(())
}
