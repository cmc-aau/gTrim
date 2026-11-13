# Changelog

## 0.8.0

### Added

- `trim --mask-file`: the columns of the full concatenated alignment that were removed before trimming, e.g. by WitChi. This is the `--witchi_info_file` input of GenomeTreeTk `msa_filter`. The file can be a WitChi `*_pruned.tsv` (the 0-based index of each removed column in an `Original Index` column) or a binary mask as used by `mask`; the format is detected from the content. The alignment then contains only the retained columns, and each marker's length is reduced by the number of its columns that were removed so that the alignment's columns are assigned to the correct markers.
  - Checks that the inputs line up: a binary mask must have one character per column described by the marker info file; WitChi indices must be within the summed marker length and listed only once; and the number of retained columns must equal the alignment length.
  - Output masks (`.trimmed.mask`, `.trimmed.filtered.mask`) cover the full alignment, as in `msa_filter`, so removed columns are `0`.
  - The output matches `msa_filter` (checked on 5,000 genomes with 2,500 columns removed), except when a removed column is the first column of a marker: `msa_filter` attributes such a column to the preceding marker (its `find_marker` compares with `<=`), so the per-marker statistics of the affected markers differ from gTrim's.
- Test data for `trim` with a mask (`tests/input_file/trim_test_witchi.*`) and integration tests for WitChi and binary masks and for mismatched inputs.
- Documentation of the WitChi pruned TSV format (`docs/src/files/witchi.tsv.rst`).

### Changed

- `trim --marker-info-file` help is shortened to one line; the file format is described in the documentation.
- `.trimmed.markers.tsv` has a new `masked_columns` column after `length`: the number of columns of the marker removed by `--mask-file` (0 without one). `masked_columns`, `retained_columns`, `gappy_columns`, `conserved_columns` and `unconserved_columns` sum to `length`.

## 0.7.2

### Changed

- The README shows the stacked gTrim logo, with "gTrim" centred below the mark, as in the gPlace and gTranslate READMEs. A dark-mode version (`docs/src/_static/gtrim_stacked_dark.svg`, with a light "g") is shown when GitHub uses a dark theme.

## 0.7.1

### Added

- Sphinx documentation in `docs/`, in the style of the gTranslate documentation, covering installation, each command and its options, the input and output files, how `trim` works and how it differs from GenomeTreeTk `msa_filter`, an FAQ, references, and this changelog. The version shown is read from `Cargo.toml`. Build it with `make html` in `docs/` (requirements in `docs/requirements.txt`).
- gTrim logo in the style of the gTranslate logo: nine alignment rows whose columns share one colour per nucleotide, with red bars marking the filtered columns, and a "gTrim" wordmark. The files (mark, icon, horizontal and stacked lockups, wordmark and favicon) are in `docs/src/_static/`; the stacked lockup and favicon are used by the documentation and the mark heads the README.
- GitHub Actions workflow (`.github/workflows/docs.yml`) that builds the documentation on every push to `main` and publishes it to the `gh_pages` branch for GitHub Pages at https://cmc-aau.github.io/gTrim/.

### Changed

- The README states more prominently that gTrim is a work in progress and links to the documentation.

## 0.7.0

### Changed

- Output file names no longer encode parameter values. Each command adds a word for the command to the input name (with `.gz`, `.faa`, `.fna` or `.fasta` removed): `.pruned` (`prune-gap`), `.trimmed` (`trim`), `.masked` (`mask`) and `.column_stats` (`column-stats`). For example, `trim` on `gtdb_r237_bac120_concatenated.masked.faa.gz` now writes `gtdb_r237_bac120_concatenated.masked.trimmed.fasta.gz` instead of `gtdb_r237_bac120_concatenated.masked_l5000_aa50_c25-95_trimmed.fasta.gz`. Companion files use `.` instead of `_`: `.trimmed.filtered.mask`, `.trimmed.markers.tsv`, `.trimmed.genomes.tsv`, `.pruned.tsv`, `.pruned.mask`.
- Test ground-truth files renamed to match; their contents are unchanged.

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
