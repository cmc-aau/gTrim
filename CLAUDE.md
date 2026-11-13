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

Subcommands: `prune-gap`, `mask`, `column-stats`. All take a single MSA via `-i/--alignment-file` and `-o/--out-dir`. Example:

```bash
cargo run --release -- prune-gap -i tests/input_file/gap_test.fna.gz -o out/ --gap-threshold 50
```

`src/cli.rs` is the source of truth for options.

## Workflow rules

- **Pull requests:** every PR must bump the version in `Cargo.toml` (which also changes the `--version` output and the version logged at startup) and add an entry for that version to `CHANGELOG.md`. Create `CHANGELOG.md` if it doesn't exist yet.
- **CURRENT_WORK.md:** keep a `CURRENT_WORK.md` file in the repo root that describes the task currently in progress. Update it as the work moves forward. It must contain enough to resume the work if the connection drops or the server crashes: the goal, what is done, the next step, the files involved, and any running jobs with their commands, output locations and how to check on them.
- **Long-running tasks:** before starting anything expected to take more than 10 minutes (for example, runs on the large alignments in `test_sets/`), state the expected total run time and the expected finish time in Vancouver time (`TZ=America/Vancouver date`).

## Architecture

- `main.rs` dispatches each subcommand. A shared `init()` creates the output directory, sets up log4rs logging to stderr and `<out_dir>/gtrim.log`, and builds the **global** rayon thread pool.
- Every command module (`prune/gap`, `prune/mask`, `stats`) follows the same pattern:
  1. `Struct::new(...)` is followed by `run()`.
  2. `run()` reads the alignment through `AlignmentReader` into `(Vec<String> ids, Array2<char>)` (rows = sequences, columns = alignment positions, residues uppercased).
  3. It calls `SequenceTypeDetector::detect` to choose DNA or amino acids. Unknown characters default to amino acids.
  4. It computes the result and writes it out.
- Output naming: `alignment_writer::get_output_prefix` strips `.gz/.faa/.fna/.fasta` from the input name and joins the result with `out_dir`. Each command adds a parameter-encoded suffix, e.g. `_gt{thresh}_m{max}_pruned.fasta` for gap. Companion files (`.tsv`, `.mask`) are derived by string-replacing `.fasta` in that suffix. If the input ends in `.gz`, the output FASTA is gzip-compressed (multithreaded via `gzp`) and gets a `.gz` appended.
- A mask file is a single line of `0`/`1` characters, one per original column, where `1` = kept. `mask` applies such a file to an alignment.

## Tests

The integration tests in `tests/` run the compiled binary via `assert_cmd` (`cargo_bin_cmd!("gtrim")`), write to a tempdir, and compare the output against `tests/ground_truth/` using the helpers in `tests/test_utils.rs`. Those helpers compare TSVs with a float tolerance and the first N columns only, compare FASTA files by ID-to-sequence map, and compare masks exactly. Inputs are in `tests/input_file/`. The larger alignments in `test_sets/` are for manual benchmarking and are not used by `cargo test`.
