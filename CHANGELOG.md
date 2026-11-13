# Changelog

## 0.6.0

### Changed

- Subcommands are listed in the order `prune-gap`, `trim`, `mask`, `column-stats`.
- `trim`: the default `--msa-length` is now 1000000, so by default every column passing the filters is retained and no subsampling is done (unless more than 1,000,000 columns pass). Pass `--msa-length 5000` for the previous behaviour, which matches GenomeTreeTk `msa_filter`.
- `trim`: the default `--seed` is now 42 (was 12345). The `trim` test ground truth was regenerated; the filtered columns are unchanged.

## 0.5.0

### Changed

- Alignments are stored as `Array2<u8>` (one byte per residue) instead of `Array2<char>`. The reader appends each sequence to a single buffer that becomes the array's storage, rather than collecting a `Vec<Vec<char>>` and copying it into the array, which kept two 4-byte copies of every residue in memory at once. Peak memory drops by about 8× for every command. On `test_sets/gtdb_n10000_test` (10,000 × 41,084, 16 threads), measured with `/usr/bin/time -v`:

  | Command | Peak RSS before → after | Wall time before → after |
  | --- | --- | --- |
  | `trim` | 3.25 GB → 0.42 GB | 4.1 s → 2.1 s |
  | `prune-gap --gap-threshold 50` | 3.25 GB → 0.43 GB | 14.9 s → 3.1 s |
  | `mask` | 3.25 GB → 0.41 GB | 5.1 s → 2.8 s |
  | `column-stats` | 3.25 GB → 0.43 GB | 19.4 s → 3.0 s |

  The output files of `trim`, `prune-gap` and `mask` are byte-identical to 0.4.0.
- `prune-gap` and `column-stats` count residues over blocks of adjacent columns in parallel (new `column_counts` module, shared with `trim`) instead of walking the alignment one column at a time.
- `column-stats` lists residues with equal counts in character order. Previously their order depended on hash map iteration and could change between runs. Its TSV is now written through a buffered writer.
- Non-ASCII bytes in an alignment are written back out unchanged. Previously each was treated as a Latin-1 character and written as UTF-8.

## 0.4.0

### Added

- `trim` subcommand: a reimplementation of `msa_filter` from GenomeTreeTk, without its optional WitChi input. It takes a concatenated marker gene alignment and a marker info file, removes gappy (`--min-perc-aa`), highly conserved (`--max-consensus`) and poorly conserved (`--min-consensus`) columns, and randomly subsamples the same number of columns from each marker to reach about `--msa-length` columns. Outputs the trimmed alignment, masks before and after subsampling, per-marker statistics and per-genome residue counts.
  - Filtered columns and the number of columns sampled per marker are identical to `msa_filter`. Sampling is seeded (`--seed`), so results are reproducible.
  - When `--msa-length` cannot be reached, every retained column is used, rather than at most 999 per marker.
- `rand` and `rand_chacha` dependencies.

## 0.3.0

### Removed

- `prune-chi2` subcommand and its module (`src/prune/chi2/`), integration test and test data.
- `mask-divergent-loci` subcommand and its module (`src/prune/divergent/`), integration test, test data and test-set script, `examples/aai_matrix.rs`, the example plots in `example_files/`, the methodology document `MASK_DIVERGENT_LOCI.md` and the evaluation scripts in `analysis/mask_divergent_loci/`.
- `src/progress.rs` and the dependencies used only by the removed commands: `csv`, `indicatif`, `rand`, `rand_chacha`, `resvg`, `serde`, `serde_json` and `flate2` (dev).

### Changed

- `README.md` and `CLAUDE.md` describe the remaining commands (`prune-gap`, `mask`, `column-stats`).

## 0.2.0

### Added

- `mask-divergent-loci` subcommand. It takes a directory of per-locus MSAs and masks (replaces with gaps) the loci of a genome that are inconsistent with vertical descent, such as contaminants, laterally transferred copies, paralogs and misaligned sequences. The method and its evaluation are documented in `MASK_DIVERGENT_LOCI.md`.
  - **Low-identity test:** a locus is masked when its identity to the genome's nearest relatives (by AAI with the locus left out) falls more than `--outlier-threshold` points below the trend for that locus, in a majority of `--neighbours` comparisons.
  - **High-identity test** (on by default; `--no-upper-mask` disables it): a locus is masked when the genomes with the most similar copies are much less related to the genome than that identity predicts (`--upper-neighbours`, `--upper-window-size`, `--upper-outlier-threshold`).
  - **Trend line:** `--method` sets how the expected value is calculated. The default, `monotone-median`, is a running median over the `--window-size` nearest points, made non-decreasing. Points with tied x values share one expected value.
  - **Overlap filters:** `--min-aligned-frac` (default 0.5) and `--min-aligned-columns` (default 50).
  - **Outputs:** masked MSAs; every comparison of both tests; masked loci with vote counts and the source of each masked copy (best match, whether the matching genomes form one lineage, and, with `--taxonomy`, their lowest common taxon); per-locus and per-genome summaries, the latter flagging suspect genomes (`--suspect-masked-perc`) and the main source of their high-identity masks; and, with `--plot`, SVG and PNG plots of each locus.
- `examples/aai_matrix.rs`: writes the genome-by-genome AAI matrix of a directory of MSAs.
- `analysis/mask_divergent_loci/`: scripts for the evaluation in `MASK_DIVERGENT_LOCI.md`.
- `CLAUDE.md`: guidance for working on the code base.

### Fixed

- The `column-stats` module (`src/stats/`) was missing from the repository, so 0.1.0 did not build from a clean checkout.
- `prune-chi2` integration test: the ground truth had been stale since the off-by-one fix in f0badba.
