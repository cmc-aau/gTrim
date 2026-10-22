# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

gTrim is a Rust CLI (crate `gtrim`, edition 2024) for pruning columns from multiple sequence alignments. Its main feature is an optimized reimplementation of [WitChi](https://github.com/stephkoest/witchi), which prunes compositionally biased columns using an empirical chi-squared test. It targets very large alignments: tens of thousands of columns and up to hundreds of thousands of taxa. Memory use and parallelism are therefore primary design concerns.

## Commands

```bash
cargo build --release                 # binary: target/release/gtrim
cargo test                            # unit tests + integration tests in tests/
cargo test --test prune_chi2_test     # one integration test file
cargo test test_verify_cli            # one test by name (clap CLI consistency check in src/cli.rs)
cargo clippy
```

Subcommands: `prune-gap`, `prune-chi2`, `mask`, `column-stats`, `mask-divergent-loci`. All take `-o/--out-dir`. All except `mask-divergent-loci` take a single MSA via `-i/--alignment-file`. `mask-divergent-loci` takes a directory of per-locus MSAs via `-i/--msa-dir`. Example:

```bash
cargo run --release -- prune-chi2 -i tests/input_file/gtdb_r226-n100.faa.gz -o out/ \
  --residues-to-prune-perc 1 --max-residue-pruned 100 --permutations 100 --touchdown -t 4 --threads-permutations 1
```

The README's usage section is outdated: it shows WitChi-style flags such as `--file` and `--num_workers_chisq`. `src/cli.rs` is the source of truth.

## Workflow rules

- **Pull requests:** every PR must bump the version in `Cargo.toml` (which also changes the `--version` output and the version logged at startup) and add an entry for that version to `CHANGELOG.md`. Create `CHANGELOG.md` if it doesn't exist yet.
- **CURRENT_WORK.md:** keep a `CURRENT_WORK.md` file in the repo root that describes the task currently in progress. Update it as the work moves forward. It must contain enough to resume the work if the connection drops or the server crashes: the goal, what is done, the next step, the files involved, and any running jobs with their commands, output locations and how to check on them.
- **Example programs:** `examples/aai_matrix.rs` writes the genome-by-genome AAI matrix used by the bucketed evaluation (`cargo build --release --examples`).
- **MASK_DIVERGENT_LOCI.md:** this file records the methodology, evaluation design and results for `mask-divergent-loci`. Update it whenever the method, its defaults or its evaluation change, or new tests are run: results, conclusions, the reproduction commands in section 6, and a change-log entry. Analysis scripts live in `analysis/mask_divergent_loci/` and write to a `WORK` directory outside the repo.
- **Long-running tasks:** before starting anything expected to take more than 10 minutes (for example, runs on the large alignments in `test_sets/`), state the expected total run time and the expected finish time in Vancouver time (`TZ=America/Vancouver date`).

## Architecture

- `main.rs` dispatches each subcommand. A shared `init()` creates the output directory, sets up log4rs logging to stderr and `<out_dir>/gtrim.log`, and builds the **global** rayon thread pool.
- Every command module (`prune/gap`, `prune/chi2`, `prune/mask`, `stats`) follows the same pattern:
  1. `Struct::new(...)` is followed by `run()`.
  2. `run()` reads the alignment through `AlignmentReader` into `(Vec<String> ids, Array2<char>)` (rows = sequences, columns = alignment positions, residues uppercased).
  3. It calls `SequenceTypeDetector::detect` to choose DNA or amino acids. Unknown characters default to amino acids.
  4. It computes the result and writes it out.
- Output naming: `alignment_writer::get_output_prefix` strips `.gz/.faa/.fna/.fasta` from the input name and joins the result with `out_dir`. Each command adds a parameter-encoded suffix, e.g. `_p{perms}_s{cols_per_iter}_m{max}[_touchdown|_rapid]_pruned.fasta` for chi2 and `_gt{thresh}_m{max}_pruned.fasta` for gap. Companion files (`.tsv`, `.mask`, `_scores.tsv`, `_score_dict.json`) are derived by string-replacing `.fasta` in that suffix. If the input ends in `.gz`, the output FASTA is gzip-compressed (multithreaded via `gzp`) and gets a `.gz` appended.
- A mask file is a single line of `0`/`1` characters, one per original column, where `1` = kept. `mask` applies such a file to an alignment.

### Chi-square pruning (`src/prune/chi2/`)

- `chi_square_pruner.rs` runs the workflow:
  1. A permutation test on the original alignment builds the null distribution of per-taxon chi² values.
  2. An iterative loop then repeats:
     - compute per-row stats;
     - compute empirical p-values;
     - check the stopping criteria (alignment p ≥ 0.95, or no significant taxa);
     - rank columns by quartic chi² difference;
     - remove the top N columns.
  3. The loop tracks `original_indices` so that pruned columns are reported in original coordinates.
- N per iteration is `residues_to_prune_perc` × the original length, with a minimum of 1. In `--touchdown` mode, N drops to 0.1% of the original length once the median chi² falls below the 99th percentile of the permuted distribution. In `--rapid-convergence` mode, N is recomputed from the current alignment length every iteration.
- `permutation_test.rs` uses its **own** rayon pool, sized by `--threads-permutations`. The global pool, sized by `-t/--threads-pruning`, is used for everything else. Each permutation thread clones the whole alignment, so peak memory scales with `--threads-permutations`. Keep it small.
- Permutations are seeded deterministically (`ChaCha8Rng`, seed `12345 + i`), so results are reproducible. Results do not depend on the thread count.
- `chi_square_calculator.rs` has both serial and `_parallel` versions of some routines. The serial ones are used inside permutation workers, which are already parallel. The `_parallel` ones use the global pool. If any character count is zero, 1 is added to every count to avoid division by zero.
- Column ranking uses the column index as a tiebreaker so that output is deterministic.
- The intentional difference from WitChi: commit f0badba fixed an off-by-one error, so gTrim prunes one fewer iteration than WitChi. That outcome is correct. Do not "fix" it to match WitChi.

### Divergent loci masking (`src/prune/divergent/`)

`mask-divergent-loci` operates on a directory of per-locus MSAs, not a single alignment. Each FASTA header starts with the genome ID. The command masks (replaces with gaps) a genome's locus when that locus is much less similar to the genome's closest neighbours than genome-wide similarity predicts.

- `locus_set.rs` concatenates all loci into one genomes × columns `Array2<u8>`. Genomes are sorted by ID. Valid residues are stored as their byte; gaps, ambiguous residues and absent loci are stored as `GAP_CODE` (0). It also records non-gap counts per (genome, locus).
- `aai.rs`: for each focal genome, in parallel, it computes per-locus (identical, aligned) counts against every other genome on the fly. Nothing is stored as an n² matrix: the work is O(n²) but memory is O(n × loci) per thread. For each locus, neighbours are ranked by **leave-one-out AAI**: pooled counts over all *other* loci, so the locus being evaluated never influences neighbour choice or its x-value. The k closest genomes that have the locus and pass the `--min-aligned-frac` overlap test (applied both to the locus and to the leave-one-out totals) and have at least `--min-aligned-columns` (default 50) aligned columns at the locus (`OverlapFilter`) are kept. Each comparison records x = leave-one-out AAI (`LocusComparison::similarity`) and y = locus identity, both in percent.
- **High-identity ("upper") test:** on by default; `--no-upper-mask` disables it. It has its own `--upper-neighbours`, `--upper-window-size` and `--upper-outlier-threshold`, passed as a separate `TestParams` (the low-identity test uses `--neighbours`, `--window-size` and `--outlier-threshold`). The `aai.rs` pass also selects the k neighbours with the highest *locus identity*, breaking ties by leave-one-out AAI. The masker fits a second per-locus trend, of leave-one-out AAI against locus identity, and flags AAI < expected − T, catching copies acquired from another lineage present in the data. `OutlierTest` in `divergent_loci_masker.rs` sets which value is x and which is y, so both tests share one evaluation path. A genome locus is masked if either test masks it.
- `outliers.rs`: within each locus, the expected y is the running median of the `--window-size` points nearest in x, made non-decreasing in x by isotonic regression (pool adjacent violators) under the default `--method monotone-median`; `--method median` skips that step. The window is extended to include every point tied with its end x values, so tied points share one expected value whatever their order. A comparison is an outlier if y < expected − `--outlier-threshold`, where the threshold is in absolute percentage points. A (genome, locus) is masked if more than 50% of its comparisons are outliers. `OutlierMethod` is a clap `ValueEnum`; new methods go there and in `expected_values`.
- `taxonomy.rs` and the masked TSV: for each masked genome locus, `divergent_loci_masked.tsv` reports the best match by locus identity, `source_coherent` (outlier neighbours of the high-identity test all ≥ `SOURCE_COHERENT_AAI` = 90% genome-wide AAI to each other) and, with `--taxonomy` (genome ID → `;`-separated ranks, GTDB style), the lowest common taxon of the matching genomes. The genome summary reports each genome's main source (most common best match among its high-identity masks) and how many of those masks come from that lineage (`main_source`). `NeighbourComparer` therefore always ranks neighbours by locus identity (k′ = 1 when `--no-upper-mask` is set).
- The genome summary's `suspect` column flags genomes with at least `--suspect-masked-perc` (default 10) percent of their loci, and at least `MIN_SUSPECT_MASKED` (3) loci, masked. It only reports; nothing is dropped.
- `plots.rs` (`--plot`): writes `plots/divergent_loci_plots{,_log}.{svg,png}`, a 7-column grid with one panel per locus. The SVG is built by hand as a string and rasterised to PNG with `resvg`, which is pure Rust and uses system fonts (DejaVu) for text. Grey points that overlap at half-pixel resolution are drawn once, and line vertices are reduced to one per pixel column, to keep the SVGs manageable (13–18 MB on the full example set). Log axes plot −log10(100 − identity), floored at 0.1 so that 100% is drawn at 99.9. Example reference figures are in `example_files/`.
- On the 77-locus, 5,813-genome example set (`/srv/projects/donovan/gtrim/data/gtdb_r232_fun77_markers`), a run takes about 25 s and 1.0 GB with 96 threads; `--plot` adds about 7 s and raises peak memory to 1.4 GB.

## Tests

The integration tests in `tests/` run the compiled binary via `assert_cmd` (`cargo_bin_cmd!("gtrim")`), write to a tempdir, and compare the output against `tests/ground_truth/` using the helpers in `tests/test_utils.rs`. Those helpers compare TSVs with a float tolerance and the first N columns only, compare FASTA files by ID-to-sequence map, and compare masks exactly. Inputs are in `tests/input_file/`. Score files that depend on permutations are not compared. The `mask-divergent-loci` test data (`tests/input_file/divergent_test/`) is a subset of the example set with one planted divergent locus. It was generated by `tests/scripts/make_divergent_test_set.py`. The exact command is in the comments in `tests/mask_divergent_loci_test.rs`. The larger alignments in `test_sets/` are for manual benchmarking and are not used by `cargo test`.
