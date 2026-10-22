#!/usr/bin/env python3
"""Distribution of G for real genome loci, by A bucket (MASK_DIVERGENT_LOCI.md section 4.10).

For each genome locus in the baseline run on the unmodified data ($WORK/details/orig_k3_n120), A is the median
leave-one-out AAI to its low-identity test neighbours (as in make_bucket_plants.py) and G is its mean locus
identity to them (the real-data analogue of the planted G). Prints, per A bucket, the number of loci, median and
low percentiles of G, and loci (masked by the low-identity test) per G bin of the Experiment 1 grid plus 95-100.
Writes 2.5-point histograms of G per A bucket (all loci and masked loci) as JSON to the path given, if any.

Usage: WORK=... python3 real_g_distribution.py [out.json]
"""
import collections, csv, json, os, statistics, sys
D = os.path.join(os.environ['WORK'], 'details', 'orig_k3_n120')
A_BINS = [(0, 70), (70, 80), (80, 85), (85, 90), (90, 95), (95, 98), (98, 100.01)]
G_BINS = [(0, 40), (40, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 95), (95, 100.01)]
def fb(bins, v): return next(i for i, (lo, hi) in enumerate(bins) if lo <= v < hi)
aai, ident = collections.defaultdict(list), collections.defaultdict(list)
with open(f'{D}/divergent_loci_comparisons.tsv') as f:
    h = f.readline().rstrip('\n').split('\t'); li, gi, ai, ii = (h.index(c) for c in ('locus', 'genome_id', 'genome_aai_loo', 'locus_identity'))
    for line in f:
        t = line.rstrip('\n').split('\t'); k = (t[li], t[gi])
        aai[k].append(float(t[ai])); ident[k].append(float(t[ii]))
masked = set()
for r in csv.DictReader(open(f'{D}/divergent_loci_masked.tsv'), delimiter='\t'):
    if r['masked_by'] in ('low_identity', 'both'): masked.add((r['locus'], r['genome_id']))
W = 2.5; NB = 40
hist = [[0] * NB for _ in A_BINS]; hist_m = [[0] * NB for _ in A_BINS]
cells = collections.Counter(); cells_m = collections.Counter(); rowG = collections.defaultdict(list)
for k in aai:
    a = statistics.median(aai[k]); g = sum(ident[k]) / len(ident[k])
    r = fb(A_BINS, a); b = min(int(g // W), NB - 1)
    hist[r][b] += 1; cells[(r, fb(G_BINS, g))] += 1; rowG[r].append(g)
    if k in masked: hist_m[r][b] += 1; cells_m[(r, fb(G_BINS, g))] += 1
out = {'hist': hist, 'hist_m': hist_m, 'rows': []}
for r, (lo, hi) in enumerate(A_BINS):
    v = sorted(rowG[r]); n = len(v)
    row = {'A': f'{lo}-{min(hi,100):g}', 'n': n, 'masked': sum(hist_m[r]), 'median_G': v[n // 2], 'p1': v[int(0.01 * n)], 'p01': v[int(0.001 * n)],
           'cells': {f'{G_BINS[c][0]}-{min(G_BINS[c][1],100):g}': [cells[(r, c)], cells_m[(r, c)]] for c in range(len(G_BINS))}}
    out['rows'].append(row)
    print(f"A {row['A']:>7}: n={n:>7} masked={row['masked']:>5} median G={row['median_G']:.1f} 1st pct={row['p1']:.1f} 0.1st pct={row['p01']:.1f}")
    print('    ' + '  '.join(f'{g}:{c[0]}({c[1]})' for g, c in row['cells'].items()))
if len(sys.argv) > 1:
    json.dump(out, open(sys.argv[1], 'w'))
