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

Subcommands: `prune-gap`, `trim`, `mask`, `column-stats`. All take a single MSA via `-i/--alignment-file` and `-o/--out-dir`; `trim` also requires `-m/--marker-info-file`. Examples:

```bash
cargo run --release -- prune-gap -i tests/input_file/gap_test.fna.gz -o out/ --gap-threshold 50
cargo run --release -- trim -i tests/input_file/trim_test.faa.gz -m tests/input_file/trim_test_markers_info.tsv -o out/ --msa-length 500
```

`src/cli.rs` is the source of truth for options. Keep each clap help string to one sentence of at most 120 characters; file formats and other detail belong in the docs, not in `--help` (the README also stays brief and links to the docs).

Known quirks: `prune-gap --max-residue-pruned` is accepted but has no effect, and `mask` and `column-stats` accept `-t/--threads` but always run with one thread (`main.rs` passes 1 to `init`).

Documentation (Sphinx, `sphinx_rtd_theme`, styled like gTranslate's docs) is in `docs/src/` and is published to https://cmc-aau.github.io/gTrim/ by `.github/workflows/docs.yml` (pushes the build to the `gh_pages` branch). Command options are written by hand in `docs/src/commands/*.rst` with `.. option::` directives, so any change to `src/cli.rs`, output files or log messages must also be made there and in `docs/src/files/`. The version is read from `Cargo.toml`; `{{ gtrim_version }}` in a page is replaced with it. The changelog page includes `CHANGELOG.md`. Build locally with `make html SPHINXOPTS="-W --keep-going"` in `docs/` (`sphinx-build` is not installed system-wide; use a venv with `docs/requirements.txt`). The logo files (`gtrim_mark.svg`, `gtrim_icon.svg`, `gtrim_lockup.svg`, `gtrim_stacked.svg`, `gtrim_wordmark.svg`, `gtrim_favicon.png`) are in `docs/src/_static/`; they follow the gTranslate logo (same wordmark letters, palette and red #A33A2C), with "Trim" in burnt orange #C1561B.

## Workflow rules

- **Pull requests:** every PR must bump the version in `Cargo.toml` (which also changes the `--version` output and the version logged at startup) and add an entry for that version to `CHANGELOG.md`. Create `CHANGELOG.md` if it doesn't exist yet.
- **CURRENT_WORK.md:** keep a `CURRENT_WORK.md` file in the repo root that describes the task currently in progress. Update it as the work moves forward. It must contain enough to resume the work if the connection drops or the server crashes: the goal, what is done, the next step, the files involved, and any running jobs with their commands, output locations and how to check on them.
- **Long-running tasks:** before starting anything expected to take more than 10 minutes (for example, runs on the large alignments in `test_sets/`), state the expected total run time and the expected finish time in Vancouver time (`TZ=America/Vancouver date`).

## Architecture

- `main.rs` dispatches each subcommand. A shared `init()` creates the output directory, sets up log4rs logging to stderr and `<out_dir>/gtrim.log`, and builds the **global** rayon thread pool.
- Every command module (`prune/gap`, `prune/trim`, `prune/mask`, `stats`) follows the same pattern:
  1. `Struct::new(...)` is followed by `run()`.
  2. `run()` reads the alignment through `AlignmentReader` into `(Vec<String> ids, Array2<u8>)` (rows = sequences, columns = alignment positions, residues as uppercased ASCII bytes). The reader appends every sequence to one row-major buffer that becomes the array's storage, so the alignment costs one byte per cell and is never copied. Keep it that way: don't introduce per-sequence `Vec`s or `char` arrays, and use byte literals (`b'-'`, `*b"ACGT"`) for residue constants.
  3. It calls `SequenceTypeDetector::detect` to choose DNA or amino acids. Unknown characters default to amino acids.
  4. It computes the result and writes it out.
- Output naming: `alignment_writer::get_output_prefix` strips `.gz/.faa/.fna/.fasta` from the input name and joins the result with `out_dir`. Each command then adds `.<word>` for the command (`.pruned`, `.trimmed`, `.masked`, `.column_stats`) and the file extension, e.g. `x.faa.gz` → `x.trimmed.fasta.gz`, `x.trimmed.mask`, `x.trimmed.filtered.mask`, `x.trimmed.markers.tsv`, `x.trimmed.genomes.tsv`. Parameter values are deliberately not encoded in file names. If the input ends in `.gz`, the output FASTA is gzip-compressed (multithreaded via `gzp`) and gets a `.gz` appended.
- Per-column statistics go through `column_counts::map_column_counts`, which counts byte values over blocks of 64 adjacent columns in parallel so each sequence is read contiguously. Walking the row-major array one column at a time (`alignment.column(i)`) is several times slower on large alignments.
- A mask file is a single line of `0`/`1` characters, one per original column, where `1` = kept. `mask` applies such a file to an alignment, and `trim --mask-file` accepts one.

### Trimming (`src/prune/trim/`)

`trim` reimplements `msa_filter` from GenomeTreeTk (`../genome_tree_tk/genometreetk/msa_filter.py`), including its optional WitChi input (`--mask-file`).

- `marker_info.rs` reads the marker info file. Header names are matched case-insensitively, so both `Marker Id` (GTDB) and `Marker ID` (GTDB-Tk) work. The marker lengths must sum to the alignment length, or, with `--mask-file`, to the length of the full alignment before masking.
- `column_mask.rs` reads `--mask-file`: a WitChi `*_pruned.tsv` (0-based `Original Index` of each removed column, header matched case-insensitively) or a one-line binary mask, detected from the content. Both describe the full concatenated alignment, i.e. the summed marker length. The input MSA must contain only the retained columns; `run()` checks this and reduces each marker's length by its removed columns. Output masks are written in full-alignment coordinates (as `msa_filter` does), and `.markers.tsv` has a `masked_columns` column.
- `msa_trimmer.rs` summarises columns with `map_column_counts`. Gaps are `-`, `.`, `_`, `*` (`TRIM_GAP_CHARS`, which differs from `GAP_CHARS` in `sequence_type_detector.rs`). A column is gappy if more than (100 − `--min-perc-aa`)% of sequences have a gap, or all do. Otherwise consensus is the count of the most common residue divided by the number of non-gap characters: below `--min-consensus` is unconserved, at or above `--max-consensus` is conserved. The fractions are computed exactly as in the Python code so that the filtered mask is identical.
- Subsampling takes the smallest k for which the sum over markers of min(retained columns, k) reaches `--msa-length`, then samples min(k, retained) columns from each marker with a `ChaCha8Rng` seeded by `--seed` (default 42), one RNG used across markers in order. Output does not depend on the thread count. The default `--msa-length` is 1,000,000, so by default every retained column is kept unless more than that pass the filters; `msa_filter`'s default is 5000.
- Intentional differences from `msa_filter`: sampling is seeded (Python is unseeded), so the sampled columns differ but the number per marker is the same; if `--msa-length` cannot be reached, every retained column is used (Python caps k at 999); `msa_filter`'s `find_marker` uses `<=`, so a WitChi-removed column that is the first column of a marker is attributed to the preceding marker (gTrim attributes it correctly); a run with no retained columns is an error; residues in `.genomes.tsv` are ASCII letters, as in Python's `isalpha()`. `--min-perc-aa` filters columns, not genomes, despite the wording of the Python help text.
- Verified against `msa_filter` on `test_sets/gtdb_n5000_test` (bac120 r226 markers info built from the GTDB-Tk r226 HMMs) with `--msa-length 5000` and with `--msa-length 12000 --min-perc-aa 70 --min-consensus 10 --max-consensus 90`: identical filtered masks, per-marker counts and sampled columns per marker. 2 s versus 26 s for Python. With a WitChi file (2,500 random interior columns removed) the outputs were also identical; with marker-start columns removed, only the affected markers' statistics differ. On `gtdb_n10000_test` (10,000 × 41,084), 2.1 s and 0.42 GB peak RSS with 16 threads (4.1 s and 3.2 GB before 0.5.0 stored residues as bytes). Nearly all of the remaining time is single-threaded gzip decompression of the input.

## Tests

The integration tests in `tests/` run the compiled binary via `assert_cmd` (`cargo_bin_cmd!("gtrim")`), write to a tempdir, and compare the output against `tests/ground_truth/` using the helpers in `tests/test_utils.rs`. Those helpers compare TSVs with a float tolerance and the first N columns only, compare FASTA files by ID-to-sequence map, and compare masks exactly. Inputs are in `tests/input_file/`. The `trim` test data (`trim_test.faa.gz`, `trim_test_markers_info.tsv`) is the first 200 genomes and first 10 bac120 markers of `test_sets/gtdb_n5000_test`; its ground truth was checked against GenomeTreeTk `msa_filter` with `--msa_length 500`. `trim_test_witchi.faa.gz` is the same alignment with the 150 columns in `trim_test_witchi.tsv` (WitChi format; `trim_test_witchi.mask` is the binary equivalent) removed, none at a marker start, so its ground truth also matches `msa_filter`. A change to an output format means updating the ground truth; check that only the intended part changed (e.g. with `cmp` against a fresh run). The larger alignments in `test_sets/` are for manual benchmarking and are not used by `cargo test`.

To compare against `msa_filter`, use `/opt/miniforge/bin/python3` (PyPy; the system `python3` is 3.6) with `../genome_tree_tk` on `PYTHONPATH`. `msa_filter` imports `biolib.seq_io`, so a small shim providing `read_fasta` may be needed.
