# gTrim

[Work in Progress. In Active Development. Recommend contacting Donovan before use.]

Toolkit for pruning multiple sequence alignments.

## Functionality

- `prune-gap`: prune columns whose percentage of gaps exceeds a threshold.
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
