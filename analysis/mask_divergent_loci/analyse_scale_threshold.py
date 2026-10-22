"""Explore data-driven outlier thresholds: T = c x local robust scale of residuals, estimated per locus over
the same N nearest points (in genome-wide AAI) used for the running median. Uses the k=3, N=120 comparison
tables (run_job.sh with KEEP=1): $WORK/details/{orig,planted}_k3_n120. Requires WORK."""
import collections, csv, os
import numpy as np

base = os.environ['WORK']
N = 120
SCALE_FLOOR = 1.0  # minimum robust SD (percentage points), avoids zero scale where identities are all ~100%
CS = [4, 5, 6, 8, 10, 12]


def load(name):
    d = os.path.join(base, 'details', name)
    by_locus = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(d, 'divergent_loci_comparisons.tsv')), delimiter='\t'):
        by_locus[r['locus']].append((float(r['genome_aai_loo']), float(r['identity_deficit']), r['genome_id']))
    return by_locus


def nearest_windows(xs, w):
    """Start index of the window of the w points nearest in x to each point (xs sorted)."""
    n = len(xs)
    w = min(w, n)
    lo = np.zeros(n, dtype=np.int64)
    s = 0
    for p in range(n):
        while s + w < n and (xs[s + w] - xs[p]) < (xs[p] - xs[s]):
            s += 1
        lo[p] = s
    return lo, w


def scaled_deficits(by_locus):
    """Return list of (locus, genome, deficit, robust scale, expected-region x) for all comparisons."""
    out = []
    for locus, pts in by_locus.items():
        pts.sort(key=lambda t: t[0])
        xs = np.array([p[0] for p in pts])
        res = np.array([p[1] for p in pts])
        lo, w = nearest_windows(xs, N)
        idx = lo[:, None] + np.arange(w)[None, :]
        scale = 1.4826 * np.median(np.abs(res[idx]), axis=1)
        scale = np.maximum(scale, SCALE_FLOOR)
        for (x, d, g), s in zip(pts, scale):
            out.append((locus, g, d, s, x))
    return out


def masked(cmps, rule):
    votes = collections.defaultdict(lambda: [0, 0])
    for locus, g, d, s, x in cmps:
        v = votes[(locus, g)]
        v[0] += 1
        v[1] += rule(d, s)
    return {k for k, (n, o) in votes.items() if 2 * o > n}


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0


q = lambda xs, p: sorted(xs)[int(p * (len(xs) - 1))]

orig = scaled_deficits(load('orig_k3_n120'))
planted = scaled_deficits(load('planted_k3_n120'))

plants = list(csv.DictReader(open(os.path.join(base, 'planted_loci.tsv')), delimiter='\t'))
plant_keys = {(p['locus'], p['genome_id']) for p in plants}

# gap of each plant (neighbour LOO AAI - donor AAI), from the planted comparisons
nb = collections.defaultdict(list)
for locus, g, d, s, x in planted:
    if (locus, g) in plant_keys:
        nb[(locus, g)].append(x)
gap = {(p['locus'], p['genome_id']): np.mean(nb[(p['locus'], p['genome_id'])]) - float(p['donor_aai']) for p in plants}
GAP_BINS = [(-100, 10), (10, 15), (15, 20), (20, 25), (25, 30), (30, 40), (40, 101)]
bin_name = lambda lo, hi: f'<{hi}' if lo < 0 else (f'{lo}+' if hi > 100 else f'{lo}-{hi}')
bins = {bin_name(lo, hi): {k for k, v in gap.items() if lo <= v < hi} for lo, hi in GAP_BINS}

print('== 1. Local robust SD of residuals (running MAD x 1.4826, N=120, floor 1.0), by genome-wide AAI of the comparison')
for lo, hi in [(0, 75), (75, 85), (85, 90), (90, 95), (95, 101)]:
    s = [c[3] for c in orig if lo <= c[4] < hi]
    print(f'  AAI {lo}-{min(hi, 100)}: n={len(s):>7}  scale median {q(s, .5):.2f}  10% {q(s, .1):.2f}  90% {q(s, .9):.2f}')
print()

print('== 2. Histogram of deficits (expected - observed) beyond 5 points, original data')
hist = collections.Counter(int(min(c[2], 99.9) // 5) * 5 for c in orig if c[2] >= 5)
for b in sorted(hist):
    print(f'  {b:>3}-{b + 5:<3} {hist[b]:>6}')
print()

print('== 3. Histogram of robust z = deficit / local scale beyond z = 3, original data')
zs = [c[2] / c[3] for c in orig]
edges = [3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50, 1e9]
for a, b in zip(edges, edges[1:]):
    print(f'  z {a:>3}-{b if b < 1e8 else "inf":<4} {sum(1 for z in zs if a <= z < b):>6}')
print(f'  (z < -3: {sum(1 for z in zs if z < -3)}; for a normal distribution, expected |z| > 4: {len(zs) * 6.3e-5:.0f} per tail)')
print()

ref = masked(orig, lambda d, s: d > 25)
ref_p = masked(planted, lambda d, s: d > 25)
rules = [('T=25 fixed', lambda d, s: d > 25)]
for c in CS:
    rules.append((f'z>{c}', lambda d, s, c=c: d > c * s))
for c in [5, 6, 8]:
    rules.append((f'z>{c} & d>10', lambda d, s, c=c: d > max(c * s, 10)))

print('== 4. Masked genome loci (k=3 majority vote, below only), planted recall, and knock-on masks')
print(f'{"rule":<14} {"masked":>7} {"J_vs_T25":>8} {"recall":>6} ' + ' '.join(f'{b:>6}' for b in bins) + f' {"knock-on":>8}')
for name, rule in rules:
    m = masked(orig, rule)
    mp = masked(planted, rule)
    rec = len(mp & plant_keys) / len(plant_keys)
    by_bin = ' '.join(f'{len(mp & k) / len(k):>6.2f}' if k else f'{"na":>6}' for k in bins.values())
    print(f'{name:<14} {len(m):>7} {jaccard(m, ref):>8.2f} {rec:>6.2f} {by_bin} {len(mp - plant_keys) - len(m):>8}')
print('  gap bins (n): ' + ', '.join(f'{b}: {len(k)}' for b, k in bins.items()))
print()

print('== 5. Implied per-comparison threshold (c x local scale) for z > 6, by genome-wide AAI')
for lo, hi in [(0, 75), (75, 85), (85, 90), (90, 95), (95, 101)]:
    t = [6 * c[3] for c in orig if lo <= c[4] < hi]
    print(f'  AAI {lo}-{min(hi, 100)}: implied T median {q(t, .5):.1f}  10% {q(t, .1):.1f}  90% {q(t, .9):.1f}')
print()

print('== 6. Per-locus implied threshold (z > 6): median over the locus of 6 x local scale')
per_locus = collections.defaultdict(list)
for locus, g, d, s, x in orig:
    per_locus[locus].append(6 * s)
meds = sorted(float(np.median(v)) for v in per_locus.values())
print(f'  loci: {len(meds)}; median implied T across loci: min {meds[0]:.1f}, 25% {q(meds, .25):.1f}, median {q(meds, .5):.1f}, 75% {q(meds, .75):.1f}, max {meds[-1]:.1f}')
