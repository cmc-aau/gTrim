#!/usr/bin/env python3
"""Create a subset of per-locus MSAs for testing the mask-divergent-loci command.

Randomly samples genomes and loci from a directory of per-locus MSAs (headers are
genome IDs) and optionally plants divergent loci by replacing the sequence of a
genome with the sequence of a donor genome for that locus.

Example:
  make_divergent_test_set.py <src_dir> <out_dir> --num-loci 6 --num-genomes 150 \
      --plant <locus>:<genome>:<donor>
"""

import argparse
import gzip
import os
import random


def read_fasta(path):
    seqs = {}
    order = []
    with gzip.open(path, 'rt') if path.endswith('.gz') else open(path) as f:
        cur_id = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                cur_id = line[1:].split()[0]
                order.append(cur_id)
                seqs[cur_id] = []
            elif cur_id is not None:
                seqs[cur_id].append(line)
    return order, {k: ''.join(v) for k, v in seqs.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('src_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--num-loci', type=int, default=None, help='number of loci to sample (default: all)')
    parser.add_argument('--shortest', action='store_true', help='select the shortest loci instead of random loci')
    parser.add_argument('--num-genomes', type=int, required=True)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--plant', action='append', default=[],
                        help='locus:genome:donor - replace genome sequence with donor sequence for locus')
    args = parser.parse_args()

    rng = random.Random(args.seed)

    files = sorted(f for f in os.listdir(args.src_dir) if f.endswith(('.faa', '.faa.gz', '.fna', '.fna.gz', '.fasta', '.fasta.gz')))
    msas = {f: read_fasta(os.path.join(args.src_dir, f)) for f in files}

    if args.num_loci is not None:
        if args.shortest:
            files = sorted(files, key=lambda f: (len(next(iter(msas[f][1].values()))), f))[:args.num_loci]
        else:
            files = sorted(rng.sample(files, args.num_loci))

    genomes = sorted(set(g for f in files for g in msas[f][0]))
    sampled = set(rng.sample(genomes, min(args.num_genomes, len(genomes))))

    plants = {}
    for p in args.plant:
        locus, genome, donor = p.split(':')
        plants[(locus, genome)] = donor

    os.makedirs(args.out_dir, exist_ok=True)
    for f in files:
        locus = f.split('.')[0]
        order, seqs = msas[f]
        out_name = f if f.endswith('.gz') else f + '.gz'
        with gzip.open(os.path.join(args.out_dir, out_name), 'wt') as out:
            for gid in order:
                if gid not in sampled:
                    continue
                seq = seqs[gid]
                donor = plants.get((locus, gid))
                if donor is not None:
                    seq = seqs[donor]
                out.write(f'>{gid}\n{seq}\n')

    print(f'Wrote {len(files)} loci and {len(sampled)} genomes to {args.out_dir}')


if __name__ == '__main__':
    main()
