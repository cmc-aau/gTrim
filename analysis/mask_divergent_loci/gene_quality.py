#!/usr/bin/env python3
"""Per-gene evidence of paralogy or lateral transfer in the marker set (MASK_DIVERGENT_LOCI.md section 4.15).

Two kinds of evidence are combined for each gene:

1. Nearest-neighbour masks from the baseline run on the unmodified data ($WORK/details/orig_k3_n120):
   genome loci masked by the low-identity test (divergent copies: paralogs, misalignment, distant contamination)
   and by the high-identity test (copies near-identical to a distant lineage: lateral transfer, contamination).

2. Pairwise discordance over the whole AAI range. Masks only compare a genome with its 3 nearest relatives, so a
   paralog or transferred copy shared by a whole clade is invisible to them. Here genome pairs are sampled in
   strata of genome-wide AAI (all strata the same size), and for each gene the locus identity of each pair is
   compared with the identity expected for that gene at the pair's leave-one-out AAI (median of the gene's pairs in
   1-point AAI bins, interpolated). A pair is discordant if |residual| > --threshold. Discordance is measured
   among pairs with leave-one-out AAI >= --min-aai, both with and without the genome loci masked in the baseline.
   Genomes whose sampled pairs are mostly discordant after masking are "problem genomes"; their clustering
   (single linkage at AAI >= --cluster-aai) shows whether a gene's residual problems sit in one lineage or many.

Usage: WORK=... python3 gene_quality.py --msa-dir DATA
Writes $WORK/gene_quality/gene_quality.tsv (one row per gene) and $WORK/gene_quality/residuals.npz.
"""

import argparse
import collections
import csv
import os
import statistics
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_bucket_plants import AaiMatrix, Loci  # noqa: E402

PAIR_BINS = [(0, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 95), (95, 98), (98, 100.01)]


def sample_pairs(m, per_bin, rng):
    n = m.shape[0]
    iu = np.triu_indices(n, k=1)
    vals = m[iu]
    pi, pj = [], []
    for lo, hi in PAIR_BINS:
        idx = np.flatnonzero((vals >= lo) & (vals < hi))
        if len(idx) > per_bin:
            idx = rng.choice(idx, per_bin, replace=False)
        pi.append(iu[0][idx]); pj.append(iu[1][idx])
        print(f'  AAI {lo}-{min(hi, 100):g}: {len(idx)} pairs', flush=True)
    return np.concatenate(pi), np.concatenate(pj)


def expected_by_bins(x, y):
    """Median y in 1-point x bins (bins with >= 20 points), linearly interpolated."""
    b = np.floor(x).astype(int)
    centers, meds = [], []
    for v in np.unique(b):
        sel = b == v
        if sel.sum() >= 20:
            centers.append(v + 0.5); meds.append(np.median(y[sel]))
    return np.interp(x, centers, meds)


def spearman(x, y):
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def single_linkage_groups(idx, m, cut):
    """Group matrix indices by single linkage at AAI >= cut; returns group sizes."""
    idx = list(idx)
    parent = list(range(len(idx)))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    sub = m[np.ix_(idx, idx)]
    for a in range(len(idx)):
        for b in np.flatnonzero(sub[a, a + 1:] >= cut) + a + 1:
            ra, rb = find(a), find(int(b))
            if ra != rb:
                parent[ra] = rb
    return sorted(collections.Counter(find(a) for a in range(len(idx))).values(), reverse=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--msa-dir', required=True)
    p.add_argument('--work', default=os.environ.get('WORK'))
    p.add_argument('--pairs-per-bin', type=int, default=15000)
    p.add_argument('--threshold', type=float, default=25.0)
    p.add_argument('--min-aai', type=float, default=60.0)
    p.add_argument('--cluster-aai', type=float, default=80.0)
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--out-dir', help='output directory (default: $WORK/gene_quality)')
    p.add_argument('--run-dir', help='mask-divergent-loci run with comparison tables (default: $WORK/details/orig_k3_n120)')
    args = p.parse_args()
    rng = np.random.RandomState(args.seed)
    out_dir = args.out_dir or os.path.join(args.work, 'gene_quality')
    os.makedirs(out_dir, exist_ok=True)
    base = args.run_dir or os.path.join(args.work, 'details', 'orig_k3_n120')

    print('Reading MSAs and AAI matrix.', flush=True)
    loci = Loci(args.msa_dir, 0.5)
    aai = AaiMatrix(os.path.join(args.work, 'aai'))
    m = aai.m
    print('Sampling genome pairs by genome-wide AAI:', flush=True)
    pi, pj = sample_pairs(m, args.pairs_per_bin, rng)
    npairs = len(pi)

    # per-gene identical and aligned counts for every sampled pair
    G = len(loci.names)
    ident = np.zeros((G, npairs), dtype=np.int32)
    align = np.zeros((G, npairs), dtype=np.int32)
    ok = np.zeros((G, npairs), dtype=bool)
    for li in range(G):
        rowmap = np.array([loci.row[li].get(g, -1) for g in aai.ids])
        ri, rj = rowmap[pi], rowmap[pj]
        has = (ri >= 0) & (rj >= 0)
        arr, ng = loci.arr[li], loci.nongap[li]
        for s in range(0, npairs, 20000):
            sl = slice(s, s + 20000)
            h = has[sl]
            a, b = arr[ri[sl][h]], arr[rj[sl][h]]
            both = (a != 0) & (b != 0)
            al = both.sum(axis=1); idn = (both & (a == b)).sum(axis=1)
            tmp_al = np.zeros(h.shape[0], dtype=np.int32); tmp_id = np.zeros(h.shape[0], dtype=np.int32)
            tmp_al[h] = al; tmp_id[h] = idn
            align[li, sl] = tmp_al; ident[li, sl] = tmp_id
            minlen = np.minimum(ng[np.maximum(ri[sl], 0)], ng[np.maximum(rj[sl], 0)])
            ok[li, sl] = h & (tmp_al > 0) & (tmp_al >= 0.5 * minlen)
    tot_id, tot_al = ident.sum(axis=0), align.sum(axis=0)

    # baseline masks and comparisons
    summary = {r['locus']: r for r in csv.DictReader(open(os.path.join(base, 'divergent_loci_summary.tsv')), delimiter='\t')}
    masked = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(base, 'divergent_loci_masked.tsv')), delimiter='\t'):
        masked[r['locus']][r['genome_id']] = r['masked_by']
    low_ident = collections.defaultdict(list)   # per gene: mean neighbour identity of low-identity masked loci
    lowkeys = {(l, g) for l, d in masked.items() for g, t in d.items() if t in ('low_identity', 'both')}
    acc = collections.defaultdict(list)
    with open(os.path.join(base, 'divergent_loci_comparisons.tsv')) as f:
        h = f.readline().rstrip('\n').split('\t')
        li_, gi_, ii_ = h.index('locus'), h.index('genome_id'), h.index('locus_identity')
        for line in f:
            t = line.split('\t', 6)
            if (t[li_], t[gi_]) in lowkeys:
                acc[(t[li_], t[gi_])].append(float(t[ii_]))
    for (l, g), v in acc.items():
        low_ident[l].append(sum(v) / len(v))

    rows = []
    resid_store = {}
    for li, name in enumerate(loci.names):
        sel = ok[li]
        x = 100.0 * (tot_id[sel] - ident[li, sel]) / np.maximum(tot_al[sel] - align[li, sel], 1)
        y = 100.0 * ident[li, sel] / align[li, sel]
        r = y - expected_by_bins(x, y)
        gi, gj = pi[sel], pj[sel]
        mk = masked.get(name, {})
        mset = np.array([aai.ids[a] in mk for a in range(len(aai.ids))])
        either_masked = mset[gi] | mset[gj]
        hi_aai = x >= args.min_aai
        def frac(cond, base_sel):
            n = base_sel.sum()
            return float((cond & base_sel).sum() / n) if n else float('nan')
        T = args.threshold
        keep = hi_aai & ~either_masked
        # problem genomes after masking: >= 4 sampled pairs, >= 50% of them discordant
        disc = keep & (np.abs(r) > T)
        npair = np.bincount(np.concatenate([gi[keep], gj[keep]]), minlength=len(aai.ids))
        ndisc = np.bincount(np.concatenate([gi[disc], gj[disc]]), minlength=len(aai.ids))
        prob = np.flatnonzero((npair >= 4) & (ndisc >= 0.5 * npair))
        groups = single_linkage_groups(prob, m, args.cluster_aai) if len(prob) else []
        s = summary[name]
        li_vals = low_ident.get(name, [])
        n_gen = int(s['num_genomes'])
        rows.append({
            'gene': name, 'genomes': n_gen, 'columns': loci.arr[li].shape[1],
            'median_length': int(np.median(loci.nongap[li])),
            'masked_low': int(s['num_masked_low_identity']), 'masked_high': int(s['num_masked_high_identity']),
            'masked_pct': 100.0 * int(s['num_masked']) / n_gen,
            'masked_low_severe': sum(v < 50 for v in li_vals),
            'pairs': int(hi_aai.sum()),
            'robust_sd': float(1.4826 * np.median(np.abs(r[hi_aai] - np.median(r[hi_aai])))),
            'disc_low_pct': 100 * frac(r < -T, hi_aai), 'disc_high_pct': 100 * frac(r > T, hi_aai),
            'disc_low_after_pct': 100 * frac(r < -T, keep), 'disc_high_after_pct': 100 * frac(r > T, keep),
            'spearman': spearman(x, y),
            'problem_genomes': len(prob), 'problem_groups': len(groups),
            'largest_group': groups[0] if groups else 0,
        })
        resid_store[name] = np.stack([x, y, r, gi, gj, either_masked]).astype(np.float32)

    fields = list(rows[0].keys())
    with open(os.path.join(out_dir, 'gene_quality.tsv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in r.items()})
    np.savez_compressed(os.path.join(out_dir, 'residuals.npz'), **resid_store)

    rows.sort(key=lambda r: -(r['masked_pct'] + r['disc_low_after_pct'] + r['disc_high_after_pct']))
    print(f'\n{"gene":<15}{"masked%":>8}{"low":>5}{"high":>5}{"sev":>5}{"sd":>6}{"dL%":>7}{"dH%":>7}{"dLa%":>7}{"dHa%":>7}{"rho":>6}{"prob":>6}{"grp":>5}{"big":>5}{"cols":>6}')
    for r in rows:
        print(f'{r["gene"]:<15}{r["masked_pct"]:>8.2f}{r["masked_low"]:>5}{r["masked_high"]:>5}{r["masked_low_severe"]:>5}'
              f'{r["robust_sd"]:>6.2f}{r["disc_low_pct"]:>7.3f}{r["disc_high_pct"]:>7.3f}{r["disc_low_after_pct"]:>7.3f}'
              f'{r["disc_high_after_pct"]:>7.3f}{r["spearman"]:>6.3f}{r["problem_genomes"]:>6}{r["problem_groups"]:>5}'
              f'{r["largest_group"]:>5}{r["columns"]:>6}')


if __name__ == '__main__':
    main()
