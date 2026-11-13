# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

gTrim is a Rust CLI (crate `gtrim`, edition 2024) for pruning columns from multiple sequence alignments. It targets very large alignments: tens of thousands of columns and up to hundreds of thousands of taxa. Memory use and parallelism are therefore primary design concerns.

## Commands

```bash
cargo build --release                 # binary: target/release/gtrim
cargo test                            # unit tests + integration tests in tests/
cargo test --test prune_gap_test      # one integration test file
cargo test test_verify_cli            # one test by name (clap CLI consistency check in src/cli.rs)
cargo clippy
```

Subcommands: `prune-gap`, `mask`, `column-stats`, `trim`. All take a single MSA via `-i/--alignment-file` and `-o/--out-dir`; `trim` also requires `-m/--marker-info-file`. Examples:

```bash
cargo run --release -- prune-gap -i tests/input_file/gap_test.fna.gz -o out/ --gap-threshold 50
cargo run --release -- trim -i tests/input_file/trim_test.faa.gz -m tests/input_file/trim_test_markers_info.tsv -o out/ --msa-length 500
```

`src/cli.rs` is the source of truth for options.

## Workflow rules

- **Pull requests:** every PR must bump the version in `Cargo.toml` (which also changes the `--version` output and the version logged at startup) and add an entry for that version to `CHANGELOG.md`. Create `CHANGELOG.md` if it doesn't exist yet.
- **CURRENT_WORK.md:** keep a `CURRENT_WORK.md` file in the repo root that describes the task currently in progress. Update it as the work moves forward. It must contain enough to resume the work if the connection drops or the server crashes: the goal, what is done, the next step, the files involved, and any running jobs with their commands, output locations and how to check on them.
- **Long-running tasks:** before starting anything expected to take more than 10 minutes (for example, runs on the large alignments in `test_sets/`), state the expected total run time and the expected finish time in Vancouver time (`TZ=America/Vancouver date`).

## Architecture

- `main.rs` dispatches each subcommand. A shared `init()` creates the output directory, sets up log4rs logging to stderr and `<out_dir>/gtrim.log`, and builds the **global** rayon thread pool.
- Every command module (`prune/gap`, `prune/mask`, `prune/trim`, `stats`) follows the same pattern:
  1. `Struct::new(...)` is followed by `run()`.
  2. `run()` reads the alignment through `AlignmentReader` into `(Vec<String> ids, Array2<char>)` (rows = sequences, columns = alignment positions, residues uppercased).
  3. It calls `SequenceTypeDetector::detect` to choose DNA or amino acids. Unknown characters default to amino acids.
  4. It computes the result and writes it out.
- Output naming: `alignment_writer::get_output_prefix` strips `.gz/.faa/.fna/.fasta` from the input name and joins the result with `out_dir`. Each command adds a parameter-encoded suffix, e.g. `_gt{thresh}_m{max}_pruned.fasta` for gap and `_l{len}_aa{perc}_c{min}-{max}_trimmed.fasta` for trim. Companion files (`.tsv`, `.mask`, `_filtered.mask`, `_markers.tsv`, `_genomes.tsv`) are derived by string-replacing `.fasta` in that suffix. If the input ends in `.gz`, the output FASTA is gzip-compressed (multithreaded via `gzp`) and gets a `.gz` appended.
- A mask file is a single line of `0`/`1` characters, one per original column, where `1` = kept. `mask` applies such a file to an alignment.

### Trimming (`src/prune/trim/`)

`trim` reimplements `msa_filter` from GenomeTreeTk (`../genome_tree_tk/genometreetk/msa_filter.py`) without its optional WitChi input.

- `marker_info.rs` reads the marker info file. Header names are matched case-insensitively, so both `Marker Id` (GTDB) and `Marker ID` (GTDB-Tk) work. The marker lengths must sum to the alignment length.
- `msa_trimmer.rs` summarises columns in parallel over blocks of 64 adjacent columns, so each sequence is read contiguously. Gaps are `-`, `.`, `_`, `*` (`TRIM_GAP_CHARS`, which differs from `GAP_CHARS` in `sequence_type_detector.rs`). A column is gappy if more than (100 − `--min-perc-aa`)% of sequences have a gap, or all do. Otherwise consensus is the count of the most common residue divided by the number of non-gap characters: below `--min-consensus` is unconserved, at or above `--max-consensus` is conserved. The fractions are computed exactly as in the Python code so that the filtered mask is identical.
- Subsampling takes the smallest k for which the sum over markers of min(retained columns, k) reaches `--msa-length`, then samples min(k, retained) columns from each marker with a `ChaCha8Rng` seeded by `--seed`, one RNG used across markers in order. Output does not depend on the thread count.
- Intentional differences from `msa_filter`: sampling is seeded (Python is unseeded), so the sampled columns differ but the number per marker is the same; if `--msa-length` cannot be reached, every retained column is used (Python caps k at 999); a run with no retained columns is an error; residues in `_genomes.tsv` are ASCII letters, as in Python's `isalpha()`. `--min-perc-aa` filters columns, not genomes, despite the wording of the Python help text.
- Verified against `msa_filter` on `test_sets/gtdb_n5000_test` (bac120 r226 markers info built from the GTDB-Tk r226 HMMs) with the defaults and with `--msa-length 12000 --min-perc-aa 70 --min-consensus 10 --max-consensus 90`: identical filtered masks, per-marker counts and sampled columns per marker. 2 s versus 26 s for Python. On `gtdb_n10000_test`, 4 s and 3.2 GB with 16 threads; most of the time and memory is the alignment reader.

## Tests

The integration tests in `tests/` run the compiled binary via `assert_cmd` (`cargo_bin_cmd!("gtrim")`), write to a tempdir, and compare the output against `tests/ground_truth/` using the helpers in `tests/test_utils.rs`. Those helpers compare TSVs with a float tolerance and the first N columns only, compare FASTA files by ID-to-sequence map, and compare masks exactly. Inputs are in `tests/input_file/`. The `trim` test data (`trim_test.faa.gz`, `trim_test_markers_info.tsv`) is the first 200 genomes and first 10 bac120 markers of `test_sets/gtdb_n5000_test`; its ground truth was checked against GenomeTreeTk `msa_filter` with `--msa_length 500`. The larger alignments in `test_sets/` are for manual benchmarking and are not used by `cargo test`.
