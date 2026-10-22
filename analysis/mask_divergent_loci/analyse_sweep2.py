"""Summarise sensitivity sweep 2 of mask-divergent-loci (AAI mode): k in {1,2,3,4,5,7} x window size, and outlier threshold sweep at N=120. Requires WORK."""
import collections, csv, os, re, statistics

base = os.environ["WORK"]
KS = [1, 2, 3, 4, 5, 7]
NS = [20, 40, 80, 120, 200, 400]
TS = [10, 15, 20, 25, 30, 35, 40]
REF = (3, 120, 25)


def run_dir(ds, k, n, t):
    d = os.path.join(base, 'runs', f'{ds}_k{k}_n{n}')
    return d if t == 25 else f'{d}_t{t}'


def masked_set(ds, k, n, t):
    with open(os.path.join(run_dir(ds, k, n, t), 'divergent_loci_masked.tsv')) as f:
        return {(r['locus'], r['genome_id']) for r in csv.DictReader(f, delimiter='\t')}


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0


def table(header, rows):
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
    line = lambda r: '  '.join(str(x).rjust(w) for x, w in zip(r, widths))
    print(line(header))
    print('  '.join('-' * w for w in widths))
    for r in rows:
        print(line(r))
    print()


plants = list(csv.DictReader(open(os.path.join(base, 'planted_loci.tsv')), delimiter='\t'))
plant_keys = {(p['locus'], p['genome_id']) for p in plants}

# gap between the planted genome's neighbours (k=3, LOO AAI) and the donor's AAI to the planted genome
nb_aai = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(base, 'details', 'planted_k3_n120', 'divergent_loci_comparisons.tsv')), delimiter='\t'):
    key = (r['locus'], r['genome_id'])
    if key in plant_keys:
        nb_aai[key].append(float(r['genome_aai_loo']))
GAP_BINS = [(-100, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 101)]
bin_name = lambda lo, hi: f'<{hi}' if lo < 0 else (f'{lo}+' if hi > 100 else f'{lo}-{hi}')
gap_bin = {}
for p in plants:
    key = (p['locus'], p['genome_id'])
    gap = statistics.mean(nb_aai[key]) - float(p['donor_aai'])
    gap_bin[key] = next(bin_name(lo, hi) for lo, hi in GAP_BINS if lo <= gap < hi)
bin_names = [bin_name(lo, hi) for lo, hi in GAP_BINS]
bin_keys = {b: {k for k, v in gap_bin.items() if v == b} for b in bin_names}

grid = [(k, n, 25) for k in KS for n in NS]
tgrid = [(k, 120, t) for k in KS for t in TS]
orig = {c: masked_set('orig', *c) for c in set(grid) | set(tgrid)}
planted = {c: masked_set('planted', *c) for c in set(grid) | set(tgrid)}

# ---------------------------------------------------------------------------
print('== A. k x window size (threshold 25): masked genome loci on the original data')
table(['k'] + [f'N={n}' for n in NS], [[k] + [len(orig[(k, n, 25)]) for n in NS] for k in KS])

freq = collections.Counter(x for c in grid for x in orig[c])
consensus = {x for x, c in freq.items() if c >= len(grid) / 2}
print(f'Consensus over the {len(grid)} k x N configurations (masked in >= half): {len(consensus)} genome loci\n')

print('== B. k comparison at N=120, threshold 25')
rows = []
for k in KS:
    c = (k, 120, 25)
    s, ps = orig[c], planted[c]
    others = [orig[(k, n, 25)] for n in NS if n != 120]
    rows.append([
        k, len(s), len({l for l, _ in s}),
        f'{jaccard(s, orig[REF]):.2f}',
        f'{statistics.mean(jaccard(s, o) for o in others):.2f}',
        f'{len(s & consensus) / len(consensus):.2f}', f'{len(s & consensus) / len(s):.2f}',
        f'{len(ps & plant_keys) / len(plant_keys):.2f}',
        len(ps - plant_keys) - len(s),
    ])
table(['k', 'masked', 'genes', 'J_vs_k3', 'J_across_N', 'cons_recall', 'cons_prec', 'planted_recall', 'excess_bg'], rows)
print('excess_bg = non-planted masks in the planted data minus masks in the original data (knock-on masks caused by planting)\n')

print('== C. Jaccard between k values (N=120, threshold 25)')
table(['k'] + [str(k) for k in KS], [[a] + [f'{jaccard(orig[(a, 120, 25)], orig[(b, 120, 25)]):.2f}' for b in KS] for a in KS])

# ---------------------------------------------------------------------------
print('== D. Threshold sweep (N=120): masked genome loci on the original data')
table(['t'] + [f'k={k}' for k in KS], [[t] + [len(orig[(k, 120, t)]) for k in KS] for t in TS])

print('== E. Threshold sweep: overall recall of planted loci')
table(['t'] + [f'k={k}' for k in KS], [[t] + [f'{len(planted[(k, 120, t)] & plant_keys) / len(plant_keys):.2f}' for k in KS] for t in TS])

print('== F. Threshold sweep: excess (knock-on) masks caused by the 300 plants')
table(['t'] + [f'k={k}' for k in KS], [[t] + [len(planted[(k, 120, t)] - plant_keys) - len(orig[(k, 120, t)]) for k in KS] for t in TS])

for k in [2, 3]:
    print(f'== G. Recall by gap between neighbour AAI and donor AAI (k={k}, N=120); n per bin: '
          + ', '.join(f'{b}: {len(bin_keys[b])}' for b in bin_names))
    table(['t'] + [f'gap {b}' for b in bin_names],
          [[t] + [f'{len(planted[(k, 120, t)] & bin_keys[b]) / len(bin_keys[b]):.2f}' if bin_keys[b] else 'na' for b in bin_names] for t in TS])

print('== H. Overlap of masked sets with the threshold-25 result (k=3, N=120, original data)')
ref = orig[(3, 120, 25)]
table(['t', 'masked', 'shared_with_t25', 'J'], [[t, len(orig[(3, 120, t)]), len(orig[(3, 120, t)] & ref), f'{jaccard(orig[(3, 120, t)], ref):.2f}'] for t in TS])

# ---------------------------------------------------------------------------
# Mirror estimate of noise: count comparisons more than t points ABOVE the trend (which divergence
# cannot produce) as an estimate of noise among those more than t points below. Only comparisons
# where both directions are possible (expected - t >= 0 and expected + t <= 100) are counted.
for k in [1, 3]:
    rows_by_key = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(base, 'details', f'orig_k{k}_n120', 'divergent_loci_comparisons.tsv')), delimiter='\t'):
        rows_by_key[(r['locus'], r['genome_id'])].append((float(r['expected_identity']), float(r['identity_deficit'])))
    devs = [d for v in rows_by_key.values() for _, d in v]
    centre = statistics.median(devs)
    mad = statistics.median(abs(d - centre) for d in devs)
    print(f'== I. Mirror noise estimate (original data, k={k}, N=120); robust SD of deviation from trend = {1.4826 * mad:.2f} points')
    rows = []
    for t in TS:
        below = above = 0
        vote_below = vote_above = 0
        for v in rows_by_key.values():
            ok = [d for e, d in v if e - t >= 0 and e + t <= 100]
            b = sum(1 for d in ok if d > t)
            a = sum(1 for d in ok if d < -t)
            below += b
            above += a
            if ok:
                vote_below += 2 * b > len(v)
                vote_above += 2 * a > len(v)
        rows.append([t, below, above, f'{above / below:.3f}' if below else 'na', vote_below, vote_above,
                     f'{vote_above / vote_below:.3f}' if vote_below else 'na'])
    table(['t', 'cmp_below', 'cmp_above', 'ratio', 'loci_below', 'loci_above', 'est_FDR'], rows)
