# gTrim

[Work in Progress. In Active Development. Recommend contacting Donovan before use.]

Toolkit for pruning multiple sequence alignments.

## Functionality

The toolkit currently only supports an optimized Rust implementation of the [WitChi python](https://github.com/stephkoest/witchi) program for pruning columns based on compositional bias using an empirical chi-squared test ([Köstlbacher et al., 2025](https://www.biorxiv.org/content/10.1101/2025.07.14.663642v1.full)). It is optimized for large alignments (tens of thousands) spanning large numbers of taxa (hundreds of thousands).

## Resource Requirements

- Peak memory is determined by the size of the multiple sequence alignment. This is replicate for each thread used to calculate permutations (--threads-permutations). Calculating permutations is not the most time consuming step so in general the number of threads should be kept small to reduce peak memory requirements.

## Usage

Compositionally biased columns can be pruned with:

```bash
gTrim chi2-prune --file <msa_file> --residues-to-prune-perc <#> --max_residue_pruned <#> --touchdown --threads-pruning <#> --threads-permutations <#> \
  --num_workers_chisq 2 --num_workers_permute 1 --residues_to_prune 2 --pruning_algorithm quartic
```

#### Options:

- `--file`: Path to the alignment file.
- `--max-residue-pruned`: Maximum columns to prune (default: 100).
- `--residues-to-prune-perc`: Percentage of biased columns to prune per iteration [0, 100] (default = 1.0). Set to 0 to prune 1 column per iteration.
- `--permutations`: Number of permutations for empirical distribution (default: 100).
- `--touchdown`: Prune columns until no significant reduction in compositional bias is observed. Recommend setting --residues_to_prune to 1-2% of the original alignment length. The number of columns to prune will be reduced to 0.1% of the original alignment length per iteration after a threshold of 99.0% is reached.
- `--threads-pruning`: Number of CPU threads for chi-square calculations (default: 1).
- `--threads-permutations`: Number of CPU threads for permutation parallelization (default: 1).
