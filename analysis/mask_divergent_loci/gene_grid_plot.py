#!/usr/bin/env python3
"""Scatter plots of each test's comparisons for all 77 genes (MASK_DIVERGENT_LOCI.md section 4.15).

One point per genome: of its k = 3 comparisons in the baseline run ($WORK/details/orig_k3_n120), the one with the
middle locus identity (the higher of 2 when only 2 exist). Trend lines use all comparisons.

Low-identity test: comparisons are each genome against its 3 nearest relatives by leave-one-out AAI;
X = leave-one-out AAI, Y = locus identity. Comparisons of genome loci
masked by this test (low_identity or both) are orange. The blue lines are the expected identity (running median,
solid) and the outlier threshold (expected - 25, dashed).

High-identity test: comparisons are each genome against its 3 most identical copies at the locus; X = locus identity,
Y = leave-one-out AAI (the axes are swapped, as in the test). Comparisons of genome loci masked by this test
(high_identity or both) are green. The blue lines are the expected AAI (solid) and the threshold (expected - 25,
dashed).

Each plot shows the number of genome loci masked by that test in its bottom right, and genes are ordered by it.
Gene IDs in grids are coloured by the gene review recommendation (red remove, blue review, black keep).
Needs matplotlib and gene_quality.tsv from gene_quality.py. Environment: RUN_DIR (default $WORK/details/orig_k3_n120)
and GENE_QUALITY (default $WORK/gene_quality/gene_quality.tsv).

Usage:
  WORK=... python3 gene_grid_plot.py grid TEST OUT.png [FIRST LAST]   grid of genes FIRST..LAST (1-based, 7 per row)
  WORK=... python3 gene_grid_plot.py panels OUT_DIR                   one PNG per gene and test plus genes.json
TEST is low or high.
"""
import collections
import csv
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

W = os.environ['WORK']
B = os.environ.get('RUN_DIR', os.path.join(W, 'details', 'orig_k3_n120'))      # run with comparison tables
Q = os.environ.get('GENE_QUALITY', os.path.join(W, 'gene_quality', 'gene_quality.tsv'))
REMOVE = ['2934341at4751', '520361at4751', '3472996at4751']
REVIEW = ['3505740at4751', '3494367at4751', '3453803at4751']
TIER_COLOUR = {'remove': '#C62828', 'review': '#1565C0', 'keep': '#1B1B1B'}
KEPT, TREND, THRESHOLD = '#8C9196', '#1E5AA8', 25.0
TESTS = {
    'low': dict(table='divergent_loci_comparisons.tsv', x='genome_aai_loo', y='locus_identity', e='expected_identity',
                by=('low_identity', 'both'), colour='#E69F00', xlim=(40, 100), ylim=(0, 100),
                xticks=[40, 50, 60, 70, 80, 90, 100], yticks=(0, 50, 100), yticks_big=(0, 25, 50, 75, 100)),
    'high': dict(table='divergent_loci_upper_comparisons.tsv', x='locus_identity', y='genome_aai_loo', e='expected_aai',
                 by=('high_identity', 'both'), colour='#009E73', xlim=(30, 100), ylim=(30, 100),
                 xticks=[30, 40, 50, 60, 70, 80, 90, 100], yticks=(40, 70, 100), yticks_big=(30, 40, 50, 60, 70, 80, 90, 100)),
}


def tier(g):
    return 'remove' if g in REMOVE else 'review' if g in REVIEW else 'keep'


def load(test):
    """Per gene: one point per genome (its middle neighbour by locus identity), the masked flag, and the expected values
    of all comparisons (for the trend lines)."""
    cfg = TESTS[test]
    masked = set()
    for r in csv.DictReader(open(os.path.join(B, 'divergent_loci_masked.tsv')), delimiter='\t'):
        if r['masked_by'] in cfg['by']:
            masked.add((r['locus'], r['genome_id']))
    comps = collections.defaultdict(list)
    with open(os.path.join(B, cfg['table'])) as f:
        h = f.readline().rstrip('\n').split('\t')
        li, gi, xi, yi, ei, ii = (h.index(c) for c in ('locus', 'genome_id', cfg['x'], cfg['y'], cfg['e'], 'locus_identity'))
        for line in f:
            t = line.rstrip('\n').split('\t')
            comps[(t[li], t[gi])].append((float(t[ii]), float(t[xi]), float(t[yi]), float(t[ei])))
    pts = collections.defaultdict(lambda: ([], [], [], [], [], []))
    for (locus, genome), c in comps.items():
        c.sort()                                  # by locus identity
        mid = c[1] if len(c) >= 2 else c[0]       # middle of 3; the higher-identity of 2
        x, y, m, e, tx, te = pts[locus]
        x.append(mid[1]); y.append(mid[2]); m.append((locus, genome) in masked); e.append(mid[3])
        tx.extend(v[1] for v in c); te.extend(v[3] for v in c)
    pts = {g: tuple(np.array(v) for v in p) for g, p in pts.items()}
    genes = [r['gene'] for r in csv.DictReader(open(Q), delimiter='\t')]
    counts = collections.Counter(l for l, _ in masked)
    order = sorted(genes, key=lambda g: (-counts[g], g))
    return pts, order, counts


def trend(x, e, width=0.5):
    """Expected value (running median) along x, reduced to one vertex per `width` points of x."""
    b = np.floor(x / width)
    keys = np.unique(b)
    return np.array([x[b == k].mean() for k in keys]), np.array([np.median(e[b == k]) for k in keys])


def draw(ax, p, n_masked, test, small):
    cfg = TESTS[test]
    x, y, m, e, ax_x, ax_e = p
    lo, hi = cfg['xlim']
    ax.plot([lo, hi], [lo, hi], ls=':', lw=0.5 if small else 1.0, color='#AAAAAA', zorder=1)
    ax.scatter(x[~m], y[~m], s=0.3 if small else 2.0, c=KEPT, alpha=0.35, linewidths=0, rasterized=True, zorder=2)
    ax.scatter(x[m], y[m], s=1.2 if small else 6, c=cfg['colour'], alpha=0.9, linewidths=0, rasterized=True, zorder=3)
    tx, ty = trend(ax_x, ax_e)                 # trend of all comparisons, as used by the test
    lw = 0.6 if small else 1.4
    ax.plot(tx, ty, ls='-', lw=lw, color=TREND, zorder=4)
    ax.plot(tx, ty - THRESHOLD, ls=(0, (3, 2)), lw=lw, color=TREND, zorder=4)
    ax.set_xlim(*cfg['xlim']); ax.set_ylim(*cfg['ylim'])
    ax.set_xticks(cfg['xticks']); ax.set_yticks(cfg['yticks'] if small else cfg['yticks_big'])
    ax.text(0.985, 0.04, f'{n_masked:,} masked', transform=ax.transAxes, ha='right', va='bottom',
            fontsize=4.8 if small else 9, color='#333333',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.8))
    for sp in ax.spines.values():
        sp.set_linewidth(0.4 if small else 0.8)


def grid(test, out, first, last):
    pts, order, counts = load(test)
    genes = order[first - 1:last]
    rows = (len(genes) + 6) // 7
    fig, axes = plt.subplots(rows, 7, figsize=(17.28, 1.41 * rows), dpi=300, sharex=True, sharey=True, squeeze=False)
    plt.subplots_adjust(left=0.025, right=0.997, top=1 - 0.035 * 6 / rows, bottom=0.045 * 6 / rows, wspace=0.06, hspace=0.32)
    for ax in axes.flat[len(genes):]:
        ax.set_visible(False)
    for i, (ax, g) in enumerate(zip(axes.flat, genes)):
        draw(ax, pts[g], counts[g], test, True)
        ax.set_title(f'{first + i}. {g}', fontsize=6.4, color=TIER_COLOUR[tier(g)], fontweight='bold', pad=1.5)
        ax.tick_params(labelsize=4.6, length=1.5, width=0.4, pad=1, labelbottom=True)
    fig.savefig(out, dpi=300, facecolor='#FBFAF7')


def panels(out_dir):
    q = {r['gene']: r for r in csv.DictReader(open(Q), delimiter='\t')}
    meta = {g: {'gene': g, 'tier': tier(g), 'genomes': int(r['genomes']), 'columns': int(r['columns']),
                'disc_after': round(float(r['disc_low_after_pct']) + float(r['disc_high_after_pct']), 3),
                'problem_genomes': int(r['problem_genomes'])} for g, r in q.items()}
    for test in TESTS:
        pts, order, counts = load(test)
        os.makedirs(os.path.join(out_dir, test), exist_ok=True)
        for rank, g in enumerate(order, 1):
            fig, ax = plt.subplots(figsize=(4.8, 3.0), dpi=150)
            plt.subplots_adjust(left=0.09, right=0.965, top=0.965, bottom=0.11)
            draw(ax, pts[g], counts[g], test, False)
            ax.tick_params(labelsize=8, length=3)
            fig.savefig(os.path.join(out_dir, test, f'{g}.png'), dpi=150, facecolor='#FBFAF7')
            plt.close(fig)
            meta[g][f'{test}_masked'] = counts[g]
            meta[g][f'{test}_rank'] = rank
    json.dump(sorted(meta.values(), key=lambda m: m['low_rank']), open(os.path.join(out_dir, 'genes.json'), 'w'))
    print(f'wrote {2 * len(meta)} panels to {out_dir}')


if __name__ == '__main__':
    if sys.argv[1] == 'grid':
        first, last = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (1, 77)
        grid(sys.argv[2], sys.argv[3], first, last)
    else:
        panels(sys.argv[2])
