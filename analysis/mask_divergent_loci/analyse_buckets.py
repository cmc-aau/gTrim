#!/usr/bin/env python3
"""Analyse the bucketed evaluation of mask-divergent-loci (MASK_DIVERGENT_LOCI.md).

Reads plants from $WORK/buckets_<experiment>/plants.tsv and the run outputs (with comparison tables) in
$WORK/details/buckets_<experiment>/runNN_k3_n120, and the control runs in $WORK/runs/buckets_<experiment>/
runNN_control_k3_n120. Knock-on masks are measured against the matching control run (same held-out sequences,
no plants).

Plants are excluded if one of their neighbours (in the comparisons of either test) was itself planted at the
same locus in the same run. For "high", the bucket is assigned from the actual best match in the run
output: B = leave-one-out AAI and I = locus identity of the rank 1 neighbour of the high-identity test.

Writes $WORK/buckets_<experiment>/results.tsv and prints summary tables.
"""

import argparse
import collections
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_bucket_plants import A_BINS, G_BINS, B_BINS, I_BINS, bin_label, find_bin  # noqa: E402

TESTS = ['low', 'high', 'either']
GNB_BINS = G_BINS + [(95, 100.01)]          # identity of the planted copy to the actual neighbours can exceed 95
DEPTH_BINS = [(0, 1), (1, 2), (2, 5), (5, 10**9)]  # donor relatives (genomes with AAI >= 90 to the donor)
DEFICIT_BINS = [(-1000, 0), (0, 10), (10, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 50), (50, 1000)]


def read_masked(run_dir):
    sets = {t: set() for t in TESTS}
    with open(os.path.join(run_dir, 'divergent_loci_masked.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            key = (r['locus'], r['genome_id'])
            sets['either'].add(key)
            if r['masked_by'] in ('low_identity', 'both'):
                sets['low'].add(key)
            if r['masked_by'] in ('high_identity', 'both'):
                sets['high'].add(key)
    return sets


def stream_comparisons(path, keys, columns):
    """Rows of a comparisons table for (locus, genome) keys, as dicts of the requested columns."""
    rows = collections.defaultdict(list)
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        li, gi = header.index('locus'), header.index('genome_id')
        idx = {c: header.index(c) for c in columns}
        for line in f:
            t = line.rstrip('\n').split('\t')
            key = (t[li], t[gi])
            if key in keys:
                rows[key].append({c: t[i] for c, i in idx.items()})
    return rows


def bootstrap_ci(by_genome, reps, rng):
    """95% CI of recall by resampling focal genomes (by_genome: list of (n, detected))."""
    if not by_genome:
        return float('nan'), float('nan')
    n = np.array([x[0] for x in by_genome], dtype=float)
    d = np.array([x[1] for x in by_genome], dtype=float)
    pick = rng.randint(0, len(n), size=(reps, len(n)))
    recall = d[pick].sum(axis=1) / n[pick].sum(axis=1)
    return float(np.percentile(recall, 2.5)), float(np.percentile(recall, 97.5))


def print_grid(title, row_bins, col_bins, values, pooled=None):
    """Print a grid of per-cell strings; `pooled` gives an extra per-row column (e.g. all cells of the row)."""
    print(title)
    print(f'{"":>8}' + ''.join(f'{bin_label(*b):>14}' for b in col_bins) + (f'{"row":>14}' if pooled else ''))
    for r, rb in enumerate(row_bins):
        print(f'{bin_label(*rb):>8}' + ''.join(f'{values.get((r, c), "-"):>14}' for c in range(len(col_bins)))
              + (f'{pooled.get(r, "-"):>14}' if pooled else ''))
    print()


def recall_grid(groups, test, with_n=True):
    """Per-key recall strings for a dict key -> list of result rows."""
    return {k: f'{sum(r["detected"][test] for r in rs) / len(rs):.2f}' + (f' ({len(rs)})' if with_n else '')
            for k, rs in groups.items() if rs}


def locus_strata(work):
    """Tertiles of loci by median locus identity between lower-test neighbours in the baseline run."""
    ident = collections.defaultdict(list)
    with open(os.path.join(work, 'details', 'orig_k3_n120', 'divergent_loci_comparisons.tsv')) as f:
        header = f.readline().rstrip('\n').split('\t')
        li, ii = header.index('locus'), header.index('locus_identity')
        for line in f:
            t = line.rstrip('\n').split('\t')
            ident[t[li]].append(float(t[ii]))
    med = sorted((float(np.median(v)), l) for l, v in ident.items())
    n = len(med)
    strata = {}
    for i, (m, l) in enumerate(med):
        strata[l] = 0 if i < n / 3 else (1 if i < 2 * n / 3 else 2)
    bounds = [med[0][0], med[n // 3][0], med[2 * n // 3][0], med[-1][0]]
    return strata, bounds


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--experiment', choices=['low', 'high'], required=True)
    p.add_argument('--work', required=True)
    p.add_argument('--thresholds', default='15,20,25,30', help='thresholds for recomputed recall of the evaluated test')
    p.add_argument('--bootstrap', type=int, default=1000)
    p.add_argument('--seed', type=int, default=1)
    args = p.parse_args()

    exp = args.experiment
    base = os.path.join(args.work, f'buckets_{exp}')
    details = os.path.join(args.work, 'details', f'buckets_{exp}')
    thresholds = [float(t) for t in args.thresholds.split(',')]
    rng = np.random.RandomState(args.seed)

    if exp == 'low':
        table, deficit_col = 'divergent_loci_comparisons.tsv', 'identity_deficit'
        row_bins, col_bins, row_name, col_name = A_BINS, G_BINS, 'A', 'G'
    else:
        table, deficit_col = 'divergent_loci_upper_comparisons.tsv', 'aai_deficit'
        row_bins, col_bins, row_name, col_name = B_BINS, I_BINS, 'B', 'I'

    plants_by_run = collections.defaultdict(list)
    with open(os.path.join(base, 'plants.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            plants_by_run[int(r['run'])].append(r)

    results = []  # per plant
    knock_on = []  # per run: (run, plants, {test: knock-on})
    ko_by_depth = {'low': collections.Counter(), 'high': collections.Counter()}
    plants_by_depth = collections.Counter()
    strata, strata_bounds = locus_strata(args.work)
    excluded = collections.Counter()
    for run in sorted(plants_by_run):
        name = f'run{run:02d}'
        run_dir = os.path.join(details, f'{name}_k3_n120')
        if not os.path.exists(os.path.join(run_dir, 'divergent_loci_masked.tsv')):
            print(f'Warning: missing output for {name}; skipping.', file=sys.stderr)
            continue
        plants = plants_by_run[run]
        keys = {(pl['locus'], pl['genome_id']) for pl in plants}
        masked = read_masked(run_dir)
        comps = stream_comparisons(os.path.join(run_dir, table), keys,
                                   ['neighbour_id', 'neighbour_rank', 'locus_identity', 'genome_aai_loo', deficit_col])
        other_table = 'divergent_loci_upper_comparisons.tsv' if exp == 'low' else 'divergent_loci_comparisons.tsv'
        other_comps = stream_comparisons(os.path.join(run_dir, other_table), keys, ['neighbour_id'])

        control = read_masked(os.path.join(args.work, 'runs', f'buckets_{exp}', f'{name}_control_k3_n120'))
        knock_on.append((run, len(plants), {t: len(masked[t] - keys) - len(control[t] - keys) for t in TESTS}))

        # attribute knock-on masks to the plant whose genome is an outlier neighbour of the masked genome locus,
        # and bin them by the sampling depth of that plant's donor lineage
        depth_of = {(pl['locus'], pl['genome_id']): int(pl['donor_relatives']) for pl in plants}
        for t, tbl in [('low', 'divergent_loci_comparisons.tsv'), ('high', 'divergent_loci_upper_comparisons.tsv')]:
            ko_keys = (masked[t] - keys) - control[t]
            if not ko_keys:
                continue
            ko_rows = stream_comparisons(os.path.join(run_dir, tbl), ko_keys, ['neighbour_id', 'outlier'])
            for k in ko_keys:
                culprits = {(k[0], r['neighbour_id']) for r in ko_rows.get(k, []) if r['outlier'] == 'true' and (k[0], r['neighbour_id']) in keys}
                if culprits:
                    d = min(depth_of[c] for c in culprits)
                    ko_by_depth[t][next(i for i, (lo, hi) in enumerate(DEPTH_BINS) if lo <= d < hi)] += 1
                else:
                    ko_by_depth[t]['unattributed'] += 1

        for pl in plants:
            key = (pl['locus'], pl['genome_id'])
            rows = comps.get(key, [])
            if not rows:
                excluded['no comparisons'] += 1
                continue
            if any((pl['locus'], r['neighbour_id']) in keys for r in rows + other_comps.get(key, [])):
                excluded['neighbour also planted'] += 1
                continue
            if exp == 'low':
                rb, cb = find_bin(A_BINS, float(pl['A'])), find_bin(G_BINS, float(pl['G']))
            else:
                top = min(rows, key=lambda r: int(r['neighbour_rank']))
                rb = find_bin(B_BINS, float(top['genome_aai_loo']))
                cb = find_bin(I_BINS, float(top['locus_identity']))
                if rb is None or cb is None:
                    excluded['best match outside bucket ranges'] += 1
                    continue
            deficits = [float(r[deficit_col]) for r in rows]
            plants_by_depth[next(i for i, (lo, hi) in enumerate(DEPTH_BINS) if lo <= int(pl['donor_relatives']) < hi)] += 1
            results.append({
                'run': run, 'locus': pl['locus'], 'genome': pl['genome_id'], 'cell': (rb, cb),
                'stratum': strata.get(pl['locus']),
                'nb_identity': float(np.mean([float(r['locus_identity']) for r in rows])),
                'detected': {t: key in masked[t] for t in TESTS},
                'deficit': float(np.mean(deficits)),
                'at_threshold': {t: 2 * sum(d > t for d in deficits) > len(deficits) for t in thresholds},
            })

    print(f'Experiment: {exp}; plants evaluated {len(results)}; excluded {dict(excluded)}\n')
    evaluated_test = exp

    # per-bucket recall with cluster bootstrap CIs
    by_cell = collections.defaultdict(list)
    for r in results:
        by_cell[r['cell']].append(r)
    out_rows = []
    grids = {t: {} for t in TESTS}
    n_grid = {}
    for cell, rs in sorted(by_cell.items()):
        genomes = collections.defaultdict(list)
        for r in rs:
            genomes[r['genome']].append(r)
        row = {row_name: bin_label(*row_bins[cell[0]]), col_name: bin_label(*col_bins[cell[1]]),
               'n': len(rs), 'n_genomes': len(genomes), 'n_loci': len({r['locus'] for r in rs}),
               'mean_deficit': f"{np.mean([r['deficit'] for r in rs]):.2f}"}
        for t in TESTS:
            recall = sum(r['detected'][t] for r in rs) / len(rs)
            lo, hi = bootstrap_ci([(len(v), sum(x['detected'][t] for x in v)) for v in genomes.values()], args.bootstrap, rng)
            row[f'recall_{t}'] = f'{recall:.3f}'
            row[f'ci_{t}'] = f'{lo:.3f}-{hi:.3f}'
            grids[t][cell] = f'{recall:.2f} ({len(rs)})'
        n_grid[cell] = f'{len(rs)}/{len(genomes)}/{len({r["locus"] for r in rs})}'
        for thr in thresholds:
            row[f'recall_{evaluated_test}_T{thr:g}'] = f"{sum(r['at_threshold'][thr] for r in rs) / len(rs):.3f}"
        out_rows.append(row)

    with open(os.path.join(base, 'results.tsv'), 'w', newline='') as f:
        if out_rows:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter='\t')
            w.writeheader()
            w.writerows(out_rows)

    by_row = collections.defaultdict(list)
    for r in results:
        by_row[r['cell'][0]].append(r)
    print_grid(f'== Plants / distinct focal genomes / distinct loci per bucket (rows {row_name}, columns {col_name})', row_bins, col_bins, n_grid,
               pooled={r: f'{len(rs)}/{len({x["genome"] for x in rs})}/{len({x["locus"] for x in rs})}' for r, rs in by_row.items()})
    for t in [evaluated_test] + [t for t in TESTS if t != evaluated_test]:
        print_grid(f'== Recall of the {t}-identity test (n)' if t != 'either' else '== Recall of either test (n)',
                   row_bins, col_bins, grids[t], pooled=recall_grid(by_row, t))
    print(f'95% CIs (resampling focal genomes) are in {os.path.join(base, "results.tsv")}\n')

    if exp == 'low':
        by_nb = collections.defaultdict(list)
        for r in results:
            cb = find_bin(GNB_BINS, r['nb_identity'])
            if cb is not None:
                by_nb[(r['cell'][0], cb)].append(r)
        print_grid('== Recall of the low-identity test by A and mean identity of the planted copy to its actual neighbours (n)',
                   A_BINS, GNB_BINS, recall_grid(by_nb, 'low'))

    print(f'== Recall of the {evaluated_test}-identity test by locus conservation stratum '
          f'(tertiles of median neighbour identity per locus in the baseline: {strata_bounds[0]:.1f}-{strata_bounds[1]:.1f}, '
          f'{strata_bounds[1]:.1f}-{strata_bounds[2]:.1f}, {strata_bounds[2]:.1f}-{strata_bounds[3]:.1f})')
    for s, name in enumerate(['variable loci', 'intermediate loci', 'conserved loci']):
        by_cell_s = collections.defaultdict(list)
        for r in results:
            if r['stratum'] == s:
                by_cell_s[r['cell']].append(r)
        print_grid(f'-- {name}: {sum(len(v) for v in by_cell_s.values())} plants, '
                   f'{len({r["locus"] for v in by_cell_s.values() for r in v})} loci', row_bins, col_bins, recall_grid(by_cell_s, evaluated_test))

    for thr in thresholds:
        g = {cell: f"{sum(r['at_threshold'][thr] for r in rs) / len(rs):.2f}" for cell, rs in by_cell.items()}
        print_grid(f'== Recall of the {evaluated_test}-identity test recomputed at threshold {thr:g}', row_bins, col_bins, g)

    print(f'== Recall of the {evaluated_test}-identity test by realised deficit (mean expected - observed over comparisons)')
    for lo, hi in DEFICIT_BINS:
        rs = [r for r in results if lo <= r['deficit'] < hi]
        if rs:
            label = f'<{hi}' if lo < -100 else (f'>={lo}' if hi > 100 else f'{lo}-{hi}')
            print(f'  {label:>7}: n={len(rs):>6}  recall={sum(r["detected"][evaluated_test] for r in rs) / len(rs):.3f}')
    print()

    print(f'== Per-locus recall (all buckets pooled; lowest {evaluated_test}-identity recall first)')
    by_locus = collections.defaultdict(list)
    for r in results:
        by_locus[r['locus']].append(r)
    loci = sorted(by_locus.items(), key=lambda kv: sum(r['detected'][evaluated_test] for r in kv[1]) / len(kv[1]))
    for locus, rs in loci[:15]:
        print(f'  {locus:<16} n={len(rs):>5}  {evaluated_test}={sum(r["detected"][evaluated_test] for r in rs) / len(rs):.3f}'
              f'  either={sum(r["detected"]["either"] for r in rs) / len(rs):.3f}')
    print()

    print('== Knock-on masks (masks of non-planted genome loci minus control masks)')
    total = sum(n for _, n, _ in knock_on)
    for t in TESTS:
        k = sum(ko[t] for _, _, ko in knock_on)
        print(f'  {t:>6}: {k:+d} over {total} plants ({1000 * k / total if total else 0:+.1f} per 1000 plants)')
    print('  per run: ' + ', '.join(f'run{run:02d} ' + '/'.join(f'{ko[t]:+d}' for t in TESTS) for run, _, ko in knock_on)
          + '  (low/high/either)')
    print('\n== Knock-on masks attributed to a plant, by sampling depth of the donor lineage (other genomes with AAI >= 90 to the donor)')
    print(f'{"donor relatives":>16} {"plants":>7} {"low":>6} {"high":>6}   (knock-on per 1000 plants: low / high)')
    for i, (lo, hi) in enumerate(DEPTH_BINS):
        label = f'{lo}' if hi == lo + 1 else (f'{lo}+' if hi > 10**8 else f'{lo}-{hi - 1}')
        n = plants_by_depth[i]
        kl, kh = ko_by_depth['low'][i], ko_by_depth['high'][i]
        print(f'{label:>16} {n:>7} {kl:>6} {kh:>6}   {1000 * kl / n if n else 0:6.1f} / {1000 * kh / n if n else 0:6.1f}')
    print(f'{"unattributed":>16} {"":>7} {ko_by_depth["low"]["unattributed"]:>6} {ko_by_depth["high"]["unattributed"]:>6}')


if __name__ == '__main__':
    main()
