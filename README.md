# gTrim

[Work in Progress. In Active Development. Recommend contacting Donovan before use.]

Toolkit for pruning multiple sequence alignments.

## Functionality

- `prune-gap`: prune columns whose percentage of gaps exceeds a threshold.
- `trim`: filter gappy, highly conserved and poorly conserved columns from a concatenated marker gene alignment, and optionally randomly subsample the same number of columns from each marker (a reimplementation of `msa_filter` in [GenomeTreeTk](https://github.com/dparks1134/GenomeTreeTk)).
- `mask`: apply a binary mask (a single line of `0`/`1` characters, one per column, `1` = kept) to an alignment.
- `column-stats`: report gap and residue statistics for each column of an alignment.

Input alignments can be FASTA files (`.faa`, `.fna`, `.fasta`), optionally gzip-compressed (`.gz`). If the input is compressed, the output alignment is also compressed.

## Usage

Prune columns with more than 50% gaps:

```bash
gtrim prune-gap -i <msa_file> -o <out_dir> --gap-threshold 50
```

Options:

- `-i/--alignment-file`: file containing the multiple sequence alignment.
- `-o/--out-dir`: output directory.
- `--gap-threshold`: maximum percentage of gaps to retain a column (default: 100).
- `--max-residue-pruned`: maximum number of columns to prune (default: 100000).
- `-t/--threads`: number of threads (default: 1).

Outputs `<prefix>_gt{threshold}_m{max}_pruned.fasta` (the pruned alignment), a `.tsv` with per-column gap statistics and a `.mask` file.

Trim a concatenated marker gene alignment:

```bash
gtrim trim -i <msa_file> -m <marker_info_file> -o <out_dir>
```

Options:

- `-i/--alignment-file`: file containing the concatenated marker gene alignment.
- `-m/--marker-info-file`: tab-separated file listing the markers in the order they are concatenated, with the columns `Marker Id`, `Description` and `Length (bp)` (e.g. `gtdb_r<release>_bac120_markers_info.tsv`). The marker lengths must sum to the alignment length.
- `-o/--out-dir`: output directory.
- `--min-perc-aa`: minimum percentage of sequences with a residue (not a gap) required to retain a column (default: 50).
- `--min-consensus`: minimum percentage of residues that must be the most common residue to retain a column, inclusive (default: 25).
- `--max-consensus`: maximum percentage of residues that can be the most common residue to retain a column, exclusive (default: 95).
- `--msa-length`: desired alignment length. The smallest number of columns per marker that gives at least this many columns is sampled from each marker. The default (1000000) retains every column that passes the filters, unless more than 1,000,000 columns pass. Use 5000 to match the default of `msa_filter`.
- `--seed`: seed for the random subsampling (default: 42).
- `-t/--threads`: number of threads (default: 1).

Gaps are `-`, `.`, `_` and `*`. Outputs, with the prefix `<input>_l{msa_length}_aa{min_perc_aa}_c{min}-{max}_trimmed`:

- `.fasta`: the trimmed alignment.
- `.mask`: columns in the trimmed alignment.
- `_filtered.mask`: columns that passed the filters, before subsampling.
- `_markers.tsv`: per-marker counts of retained, sampled, gappy, conserved and unconserved columns.
- `_genomes.tsv`: number and percentage of residues in each sequence of the trimmed alignment.

Apply a mask to an alignment:

```bash
gtrim mask -i <msa_file> -m <mask_file> -o <out_dir>
```

Outputs `<prefix>_masked.fasta`.

Calculate column statistics:

```bash
gtrim column-stats -i <msa_file> -o <out_dir>
```

Outputs `<prefix>_column_stats.tsv`.
