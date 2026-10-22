#!/usr/bin/env python3
"""Create planted datasets for the bucketed evaluation of mask-divergent-loci (MASK_DIVERGENT_LOCI.md).

Each plant replaces the sequence of a focal genome f at locus L with the aligned sequence of a donor genome D
from another lineage. Plants are assigned to buckets using values measured before planting, and buckets are
filled to a target size over successive runs (planted datasets).

Experiment "low" (low-identity test):
  A = median leave-one-out AAI of f to its lower-test neighbours at L (baseline run)
  G = identity between the planted (donor) copy and the native copy of f at L
  Donors must be foreign: AAI(f, D) <= A - 5 and G <= A - 5. Cells are used if the upper bound of the
  G bin is <= the upper bound of the A bin - 5. The donor's own sequence at L is removed, as the source of a
  contaminant is rarely in the dataset (a donor that keeps its copy sees its own sequence in a foreign genome
  and can be masked by the high-identity test).

Experiment "high" (high-identity test):
  B = AAI between f and the best remaining match of the planted copy at L
  I = identity of the planted copy to that best remaining match
  The donor's sequence at L, and every sequence at L more identical to the planted copy than the upper bound
  of the I bin, are removed (held out) so the best remaining match falls in the I bin (at most --max-holdouts per
  plant, and at most --max-locus-holdout-frac of a locus's sequences per run). Donors must be foreign:
  AAI(f, D) <= A - 10.

For both experiments, a control dataset with the same holdouts but no plants is written per run, and is
used to measure knock-on masks (masks of other genome loci caused by the plants).

Across runs, unused (locus, focal genome) pairs are preferred, so genomes are reused at different loci rather than
the same locus with a different donor; each bucket takes at most --max-locus-per-bucket plants from any one locus.

Constraints within a run: at most --max-locus-frac of the genomes with a locus are planted at that locus,
at most --max-per-genome plants per focal genome, focal and donor genomes are disjoint, each donor copy is
planted at most once per run, donors planted at the same locus in a run are from different lineages (copies
below --max-donor-identity and genomes below --max-donor-aai to each other; foreign copies from one lineage
in several genomes make sparsely sampled members of that lineage look contaminated to the high-identity
test), genomes that are
low-identity test neighbours of each other at a locus (baseline run) are not both planted there, and (high) holdouts
never remove a planted locus or another plant's best match. Focal loci and donor copies masked in the
baseline run are not used.

Outputs ($WORK/buckets_<experiment>/): runNN_msas/ and runNN_control_msas/, plants.tsv, job lists for
run_job.sh (jobs.txt: planted datasets, run with KEEP=1 as comparison tables are needed; jobs_control.txt:
control datasets, run without KEEP), and a bucket fill summary on stdout.
"""

import argparse
import collections
import gzip
import os
import random
import statistics
import sys

import numpy as np

AA = b'ACDEFGHIKLMNPQRSTVWY'
DNA = b'ACGT'

A_BINS = [(0, 70), (70, 80), (80, 85), (85, 90), (90, 95), (95, 98), (98, 100.01)]
G_BINS = [(0, 40), (40, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 95)]
B_BINS = [(0, 50), (50, 60), (60, 70), (70, 80), (80, 90)]
I_BINS = [(0, 80), (80, 90), (90, 95), (95, 98), (98, 100.01)]


def bin_label(lo, hi):
    if lo == 0:
        return f'<{hi:g}'
    return f'{lo:g}-{min(hi, 100):g}'


def find_bin(bins, v):
    for i, (lo, hi) in enumerate(bins):
        if lo <= v < hi:
            return i
    return None


class Loci:
    """Per-locus MSAs: raw sequences for writing and encoded arrays for identity calculations."""

    def __init__(self, msa_dir, min_aligned_frac):
        self.min_aligned_frac = min_aligned_frac
        self.files = sorted(f for f in os.listdir(msa_dir)
                            if f.endswith(('.faa', '.fna', '.fasta', '.faa.gz', '.fna.gz', '.fasta.gz')))
        self.names, self.order, self.seqs, self.arr, self.row, self.nongap = [], [], [], [], [], []
        for f in self.files:
            order, seqs, cur = [], {}, None
            fh = gzip.open(os.path.join(msa_dir, f), 'rt') if f.endswith('.gz') else open(os.path.join(msa_dir, f))
            for line in fh:
                line = line.strip()
                if line.startswith('>'):
                    cur = line[1:].split()[0]
                    order.append(cur)
                    seqs[cur] = []
                elif cur is not None:
                    seqs[cur].append(line)
            seqs = {k: ''.join(v) for k, v in seqs.items()}
            self.names.append(f.split('.')[0])
            self.order.append(order)
            self.seqs.append(seqs)

        # residue alphabet (amino acid unless every locus looks like DNA)
        sample = set(''.join(next(iter(s.values())) for s in self.seqs).upper().encode())
        alphabet = DNA if sample <= set(b'ACGTN-.') else AA
        valid = np.zeros(256, dtype=bool)
        valid[list(alphabet)] = True

        for order, seqs in zip(self.order, self.seqs):
            raw = np.frombuffer(''.join(seqs[g] for g in order).upper().encode(), dtype=np.uint8)
            arr = np.where(valid[raw], raw, 0).astype(np.uint8).reshape(len(order), -1)
            self.arr.append(arr)
            self.row.append({g: i for i, g in enumerate(order)})
            self.nongap.append((arr != 0).sum(axis=1))

    def identity(self, li, query_row):
        """Identity (%) of a sequence (row index at locus li) to all sequences at the locus; NaN if insufficient overlap."""
        arr, q = self.arr[li], self.arr[li][query_row]
        both = (arr != 0) & (q != 0)
        aligned = both.sum(axis=1)
        identical = (both & (arr == q)).sum(axis=1)
        min_len = np.minimum(self.nongap[li], self.nongap[li][query_row])
        ok = (aligned > 0) & (aligned >= self.min_aligned_frac * min_len)
        out = np.full(arr.shape[0], np.nan)
        out[ok] = 100.0 * identical[ok] / aligned[ok]
        return out


def read_baseline(baseline_dir):
    """Median LOO AAI of each (locus, genome) to its lower-test neighbours, neighbour relationships
    (in either direction) between genomes at each locus, and masked (locus, genome) pairs."""
    aai = collections.defaultdict(list)
    related = collections.defaultdict(set)
    with open(os.path.join(baseline_dir, 'divergent_loci_comparisons.tsv')) as f:
        header = f.readline().rstrip('\n').split('\t')
        li, gi, ni, ai = (header.index(c) for c in ('locus', 'genome_id', 'neighbour_id', 'genome_aai_loo'))
        for line in f:
            t = line.rstrip('\n').split('\t')
            aai[(t[li], t[gi])].append(float(t[ai]))
            related[(t[li], t[gi])].add(t[ni])
            related[(t[li], t[ni])].add(t[gi])
    a_value = {k: statistics.median(v) for k, v in aai.items()}

    masked = set()
    with open(os.path.join(baseline_dir, 'divergent_loci_masked.tsv')) as f:
        f.readline()
        for line in f:
            t = line.split('\t')
            masked.add((t[0], t[1]))
    return a_value, related, masked


class AaiMatrix:
    def __init__(self, aai_dir):
        self.ids = open(os.path.join(aai_dir, 'genome_ids.txt')).read().split()
        self.idx = {g: i for i, g in enumerate(self.ids)}
        n = len(self.ids)
        self.m = np.fromfile(os.path.join(aai_dir, 'aai_matrix.f32'), dtype='<f4').reshape(n, n)

    def index_array(self, genomes):
        """Matrix indices of a list of genomes (-1 if unknown)."""
        return np.array([self.idx.get(g, -1) for g in genomes], dtype=np.int64)

    def row(self, genome, index_array):
        """AAI of genome to each genome given by an index array (NaN if unknown)."""
        r = self.m[self.idx[genome]][index_array].astype(float)
        r[index_array < 0] = np.nan
        return r


def write_dataset(loci, out_dir, replace, holdouts):
    """Write MSAs with planted sequences (replace: (locus_idx, genome) -> donor) and held-out sequences removed."""
    os.makedirs(out_dir, exist_ok=True)
    for li, f in enumerate(loci.files):
        name = f if f.endswith('.gz') else f + '.gz'
        with gzip.open(os.path.join(out_dir, name), 'wt', compresslevel=1) as out:
            seqs = loci.seqs[li]
            for g in loci.order[li]:
                if (li, g) in holdouts:
                    continue
                donor = replace.get((li, g))
                out.write(f'>{g}\n{seqs[donor] if donor else seqs[g]}\n')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--experiment', choices=['low', 'high'], required=True)
    p.add_argument('--msa-dir', required=True, help='original per-locus MSAs')
    p.add_argument('--work', required=True, help='working directory (WORK)')
    p.add_argument('--baseline', help='baseline run with comparison tables (default: $WORK/details/orig_k3_n120)')
    p.add_argument('--aai-dir', help='AAI matrix from examples/aai_matrix.rs (default: $WORK/aai)')
    p.add_argument('--target', type=int, default=500, help='plants per bucket')
    p.add_argument('--target-overrides', default='',
                   help='comma-separated row bin targets, e.g. "<70=200" (rows: A bins for low, B bins for high)')
    p.add_argument('--max-runs', type=int, default=30)
    p.add_argument('--max-locus-frac', type=float, default=0.0025,
                   help='maximum fraction of the genomes with a locus planted at that locus in a run')
    p.add_argument('--max-per-genome', type=int, default=1,
                   help='plants per focal genome per run (a planted locus lowers the genome-wide AAI used at its other loci)')
    p.add_argument('--max-holdouts', type=int, default=100, help='maximum sequences removed per plant (high)')
    p.add_argument('--max-locus-holdout-frac', type=float, default=0.02,
                   help='maximum fraction of a locus\'s sequences removed per run (high)')
    p.add_argument('--max-locus-per-bucket', type=int, default=15, help='maximum plants from one locus in a bucket')
    p.add_argument('--donor-relatives-aai', type=float, default=90.0,
                   help='AAI above which another genome counts as a relative of the donor (reported as donor_relatives)')
    p.add_argument('--donor-tries', type=int, default=10, help='donors tried per focal locus')
    p.add_argument('--max-donor-identity', type=float, default=90.0,
                   help='maximum identity between donor copies planted at the same locus within a run')
    p.add_argument('--max-donor-aai', type=float, default=80.0,
                   help='maximum AAI between donor genomes planted at the same locus within a run')
    p.add_argument('--min-aligned-frac', type=float, default=0.5)
    p.add_argument('--seed', type=int, default=1)
    args = p.parse_args()
    np.seterr(invalid='ignore')  # comparisons with NaN identities are expected

    baseline = args.baseline or os.path.join(args.work, 'details', 'orig_k3_n120')
    aai_dir = args.aai_dir or os.path.join(args.work, 'aai')
    out_base = os.path.join(args.work, f'buckets_{args.experiment}')
    if os.path.exists(os.path.join(out_base, 'plants.tsv')):
        sys.exit(f'{out_base}/plants.tsv already exists; remove {out_base} to recreate the plants.')
    os.makedirs(out_base, exist_ok=True)

    print('Reading MSAs, baseline and AAI matrix.', flush=True)
    loci = Loci(args.msa_dir, args.min_aligned_frac)
    a_value, related, masked = read_baseline(baseline)
    aai = AaiMatrix(aai_dir)
    locus_idx = {n: i for i, n in enumerate(loci.names)}
    aai_index = [aai.index_array(o) for o in loci.order]

    eligible = [(locus_idx[l], g) for (l, g) in a_value if l in locus_idx and (l, g) not in masked and g in aai.idx]
    eligible.sort()
    print(f'{len(eligible)} eligible focal genome loci.', flush=True)

    if args.experiment == 'low':
        row_bins, col_bins = A_BINS, G_BINS
        cells = {(a, c) for a, (alo, ahi) in enumerate(A_BINS) for c, (clo, chi) in enumerate(G_BINS) if chi <= ahi - 5}
    else:
        row_bins, col_bins = B_BINS, I_BINS
        cells = {(b, c) for b in range(len(B_BINS)) for c in range(len(I_BINS))}
    row_target = {}
    for item in filter(None, args.target_overrides.split(',')):
        label, value = item.split('=')
        labels = [bin_label(*b) for b in row_bins]
        if label not in labels:
            sys.exit(f'Unknown row bin "{label}" in --target-overrides; expected one of: {", ".join(labels)}')
        row_target[labels.index(label)] = int(value)
    target = {c: row_target.get(c[0], args.target) for c in cells}
    remaining = dict(target)

    # number of eligible focal genome loci in each A bin (used to fill the scarcest buckets first)
    supply = collections.Counter(find_bin(A_BINS, a_value[(loci.names[li], g)]) for li, g in eligible)

    plants_out = open(os.path.join(out_base, 'plants.tsv'), 'w')
    if args.experiment == 'low':
        plants_out.write('run\tlocus\tgenome_id\tdonor_id\tA\tA_bin\tG\tG_bin\tdonor_aai\tdonor_relatives\n')
    else:
        plants_out.write('run\tlocus\tgenome_id\tdonor_id\tA\tdonor_aai\tdonor_relatives\tbest_match_id\tB_design\t'
                         'B_bin_design\tI_design\tI_bin_design\tnum_holdouts\tholdouts\n')
    jobs = open(os.path.join(out_base, 'jobs.txt'), 'w')
    jobs_control = open(os.path.join(out_base, 'jobs_control.txt'), 'w')
    used = set()
    used_focal = set()                      # (locus, genome) pairs planted in any run
    bucket_locus = collections.Counter()    # plants per (cell, locus) over all runs

    def donor_relatives(d):
        """Number of other genomes with AAI >= --donor-relatives-aai to the donor."""
        r = aai.m[aai.idx[d]]
        return int(np.sum(r >= args.donor_relatives_aai)) - 1

    for run in range(1, args.max_runs + 1):
        if all(v <= 0 for v in remaining.values()):
            break
        rng = random.Random(args.seed * 1000 + run)
        order = eligible[:]
        rng.shuffle(order)
        if args.experiment == 'low':
            # process A bins with the most remaining plants relative to their supply first, and within a bin
            # focal genome loci that have not been planted in an earlier run
            need = collections.Counter()
            for (r, c), v in remaining.items():
                need[r] += max(v, 0)
            priority = {r: need[r] / supply[r] if supply[r] else 0 for r in range(len(A_BINS))}
            order.sort(key=lambda lg: (-priority[find_bin(A_BINS, a_value[(loci.names[lg[0]], lg[1])])], lg in used_focal))
        else:
            order.sort(key=lambda lg: lg in used_focal)

        per_locus = collections.Counter()
        cap = [max(1, int(args.max_locus_frac * len(o))) for o in loci.order]
        holdout_budget = [max(1, int(args.max_locus_holdout_frac * len(o))) for o in loci.order]
        per_locus_holdouts = collections.Counter()
        per_genome = collections.Counter()
        focal_genomes, donor_genomes = set(), set()
        replace, holdouts, reserved = {}, set(), set()
        donor_loci = set()  # each donor copy is planted at most once per run
        donor_rows = collections.defaultdict(list)  # rows of donor copies planted at each locus in this run
        donor_ids = collections.defaultdict(list)   # matrix indices of donor genomes planted at each locus in this run

        def unrelated_donor(li, ident_to_donor, donor):
            """True if a donor is from a different lineage to all donors already planted at the locus: its copy is
            below --max-donor-identity and its genome below --max-donor-aai to each of them."""
            rows = donor_rows[li]
            if not rows:
                return True
            if np.nanmax(np.r_[ident_to_donor[rows], -1.0]) >= args.max_donor_identity:
                return False
            return np.nanmax(np.r_[aai.m[aai.idx[donor], donor_ids[li]], -1.0]) < args.max_donor_aai
        n_run = 0

        for li, f in order:
            if all(v <= 0 for v in remaining.values()):
                break
            if per_locus[li] >= cap[li] or per_genome[f] >= args.max_per_genome or f in donor_genomes:
                continue
            if (li, f) in holdouts or (li, f) in reserved:
                continue
            # avoid planting genomes that are neighbours of each other at the same locus
            if any((li, g) in replace for g in related.get((loci.names[li], f), ())):
                continue
            a = a_value[(loci.names[li], f)]
            row_f = loci.row[li][f]
            genomes = loci.order[li]

            if args.experiment == 'low':
                a_bin = find_bin(A_BINS, a)
                open_cols = [c for (r, c) in cells if r == a_bin and remaining[(r, c)] > 0
                             and bucket_locus[((r, c), li)] < args.max_locus_per_bucket]
                if not open_cols:
                    continue
                ident = loci.identity(li, row_f)
                donor_aai = aai.row(f, aai_index[li])
                ok = (~np.isnan(ident)) & (~np.isnan(donor_aai)) & (donor_aai <= a - 5) & (ident <= a - 5)
                candidates = collections.defaultdict(list)
                for j in np.flatnonzero(ok):
                    d = genomes[j]
                    if d == f or d in focal_genomes or (loci.names[li], d) in masked or (li, f, d) in used \
                            or (li, d) in donor_loci:
                        continue
                    c = find_bin(G_BINS, ident[j])
                    if c in open_cols:
                        candidates[c].append(j)
                if not candidates:
                    continue
                best = max(remaining[(a_bin, c)] for c in candidates)
                c = rng.choice(sorted(c for c in candidates if remaining[(a_bin, c)] == best))
                j = None
                for k in rng.sample(candidates[c], min(len(candidates[c]), args.donor_tries)):
                    if unrelated_donor(li, loci.identity(li, k), genomes[k]):
                        j = k
                        break
                if j is None:
                    continue
                d = genomes[j]
                if (li, d) in replace or (li, d) in reserved:
                    continue
                holdouts.add((li, d))
                plants_out.write(f'{run}\t{loci.names[li]}\t{f}\t{d}\t{a:.4f}\t{bin_label(*A_BINS[a_bin])}\t'
                                 f'{ident[j]:.4f}\t{bin_label(*G_BINS[c])}\t{donor_aai[j]:.4f}\t{donor_relatives(d)}\n')
                remaining[(a_bin, c)] -= 1
                bucket_locus[((a_bin, c), li)] += 1
            else:
                if per_locus_holdouts[li] >= holdout_budget[li]:
                    continue  # every plant holds out at least the donor copy, so none fits once the budget is spent
                open_cells = [cell for cell in cells if remaining[cell] > 0]
                donor_aai = aai.row(f, aai_index[li])
                donor_ok = (~np.isnan(donor_aai)) & (donor_aai <= a - 10) & (donor_aai < B_BINS[-1][1])
                pool = [j for j in np.flatnonzero(donor_ok)
                        if genomes[j] != f and genomes[j] not in focal_genomes
                        and (loci.names[li], genomes[j]) not in masked and (li, genomes[j]) not in reserved
                        and (li, genomes[j]) not in holdouts and (li, f, genomes[j]) not in used]
                if not pool:
                    continue

                # prefer donors in the B bin with the most remaining plants
                b_need = collections.Counter()
                for (b, c) in open_cells:
                    b_need[b] += remaining[(b, c)]
                by_bin = collections.defaultdict(list)
                for j in pool:
                    by_bin[find_bin(B_BINS, donor_aai[j])].append(j)
                tries = []
                for b, _ in b_need.most_common():
                    js = by_bin.get(b, [])
                    tries.extend(rng.sample(js, min(len(js), args.donor_tries - len(tries))))
                    if len(tries) >= args.donor_tries:
                        break

                best = None
                for j in tries:
                    ident = loci.identity(li, j)
                    if not unrelated_donor(li, ident, genomes[j]):
                        continue
                    ident[row_f] = np.nan  # the focal copy is replaced by the planted copy
                    ident[j] = np.nan      # the donor copy is held out
                    for c, (ilo, ihi) in enumerate(I_BINS):
                        if not any(cc == c and remaining[(b, cc)] > 0 for (b, cc) in cells):
                            continue
                        hold = {j} | set(np.flatnonzero(ident >= ihi).tolist())
                        new_holdouts = sum(1 for k in hold if (li, genomes[k]) not in holdouts)
                        if len(hold) > args.max_holdouts or per_locus_holdouts[li] + new_holdouts > holdout_budget[li]:
                            continue
                        if any((li, genomes[k]) in reserved or (li, genomes[k]) in replace for k in hold):
                            continue
                        in_bin = np.flatnonzero((ident >= ilo) & (ident < ihi))
                        in_bin = [k for k in in_bin if (li, genomes[k]) not in holdouts and (li, genomes[k]) not in replace]
                        if not in_bin:
                            continue
                        r = max(in_bin, key=lambda k: (ident[k], -k))
                        b_val = aai.m[aai.idx[f], aai.idx[genomes[r]]]
                        b = find_bin(B_BINS, b_val) if not np.isnan(b_val) else None
                        if b is None or remaining.get((b, c), 0) <= 0 or bucket_locus[((b, c), li)] >= args.max_locus_per_bucket:
                            continue
                        score = remaining[(b, c)]
                        if best is None or score > best[0]:
                            best = (score, j, c, b, r, b_val, ident[r], sorted(hold))
                if best is None:
                    continue
                _, j, c, b, r, b_val, i_val, hold = best
                d = genomes[j]
                hold_ids = [genomes[k] for k in hold]
                plants_out.write(f'{run}\t{loci.names[li]}\t{f}\t{d}\t{a:.4f}\t{donor_aai[j]:.4f}\t{donor_relatives(d)}\t'
                                 f'{genomes[r]}\t{b_val:.4f}\t{bin_label(*B_BINS[b])}\t{i_val:.4f}\t{bin_label(*I_BINS[c])}\t'
                                 f'{len(hold_ids)}\t{";".join(hold_ids)}\n')
                remaining[(b, c)] -= 1
                bucket_locus[((b, c), li)] += 1
                per_locus_holdouts[li] += sum(1 for g in hold_ids if (li, g) not in holdouts)
                holdouts.update((li, g) for g in hold_ids)
                reserved.add((li, genomes[r]))

            replace[(li, f)] = d
            used_focal.add((li, f))
            donor_loci.add((li, d))
            donor_rows[li].append(loci.row[li][d])
            donor_ids[li].append(aai.idx[d])
            used.add((li, f, d))
            focal_genomes.add(f)
            donor_genomes.add(d)
            per_locus[li] += 1
            per_genome[f] += 1
            n_run += 1

        if n_run == 0:
            print(f'Run {run}: no further plants possible; stopping.', flush=True)
            break

        name = f'run{run:02d}'
        write_dataset(loci, os.path.join(out_base, f'{name}_msas'), replace, holdouts)
        jobs.write(f'buckets_{args.experiment}/{name} 3 120 25\n')
        write_dataset(loci, os.path.join(out_base, f'{name}_control_msas'), {}, holdouts)
        jobs_control.write(f'buckets_{args.experiment}/{name}_control 3 120 25\n')
        plants_out.flush()
        jobs.flush()
        jobs_control.flush()
        filled = sum(target[c] - max(v, 0) for c, v in remaining.items())
        print(f'Run {run}: {n_run} plants; buckets filled {filled}/{sum(target.values())} '
              f'({sum(1 for v in remaining.values() if v <= 0)}/{len(cells)} complete).', flush=True)

    plants_out.close()
    jobs.close()
    jobs_control.close()

    print('\nPlants per bucket (rows: ' + ('A' if args.experiment == 'low' else 'B') + ', columns: '
          + ('G' if args.experiment == 'low' else 'I') + '; - = bucket not used)')
    print('        ' + ''.join(f'{bin_label(*b):>8}' for b in col_bins))
    for r, rb in enumerate(row_bins):
        vals = [f'{target[(r, c)] - remaining[(r, c)]:>8}' if (r, c) in cells else f'{"-":>8}' for c in range(len(col_bins))]
        print(f'{bin_label(*rb):>8}' + ''.join(vals))


if __name__ == '__main__':
    main()
