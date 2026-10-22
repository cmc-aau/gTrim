"""What-if: make each locus's expected-value trend non-decreasing (isotonic regression of the running medians)."""
import collections, csv, numpy as np
import os
B=os.path.join(os.environ['WORK'], 'details', 'orig_k3_n120')
def pava(v):
    """Non-decreasing least-squares fit (pool adjacent violators), equal weights."""
    vals, wts, lens = [], [], []
    for x in v:
        vals.append(x); wts.append(1.0); lens.append(1)
        while len(vals) > 1 and vals[-2] > vals[-1]:
            w = wts[-2] + wts[-1]; m = (vals[-2]*wts[-2] + vals[-1]*wts[-1]) / w
            n = lens[-2] + lens[-1]; vals[-2:] = [m]; wts[-2:] = [w]; lens[-2:] = [n]
    return np.repeat(vals, lens)
def run(table, xcol, ycol, ecol):
    rows = collections.defaultdict(list)
    with open(f'{B}/{table}') as f:
        h = f.readline().rstrip('\n').split('\t'); li, gi, xi, yi, ei = (h.index(c) for c in ('locus','genome_id',xcol,ycol,ecol))
        for line in f:
            t = line.rstrip('\n').split('\t'); rows[t[li]].append((float(t[xi]), float(t[yi]), float(t[ei]), t[gi]))
    cur, mono, drops = set(), set(), {}
    for locus, r in rows.items():
        r.sort(key=lambda v: v[0])
        x = np.array([v[0] for v in r]); y = np.array([v[1] for v in r]); e = np.array([v[2] for v in r]); g = [v[3] for v in r]
        em = pava(e)
        drops[locus] = (float(np.max(np.maximum.accumulate(e) - e)), float(np.abs(em - e).mean()))
        for exp, out in ((e, cur), (em, mono)):
            votes = collections.defaultdict(lambda: [0, 0])
            for gg, o in zip(g, y < exp - 25):
                votes[gg][0] += o; votes[gg][1] += 1
            out.update((locus, gg) for gg, (o, n) in votes.items() if 2 * o > n)
    return cur, mono, drops
for name, args in [('low-identity', ('divergent_loci_comparisons.tsv','genome_aai_loo','locus_identity','expected_identity')),
                   ('high-identity', ('divergent_loci_upper_comparisons.tsv','locus_identity','genome_aai_loo','expected_aai'))]:
    cur, mono, drops = run(*args)
    md = sorted(drops.items(), key=lambda kv: -kv[1][0])
    print(f'{name}: current masks {len(cur)}, non-decreasing {len(mono)}, added {len(mono-cur)}, removed {len(cur-mono)}')
    print('   largest dip below an earlier trend value (points):', [(l, round(d, 1)) for l, (d, _) in md[:5]])
    print('   median dip over genes %.1f; median mean change in expected %.2f points' % (np.median([d for d, _ in drops.values()]), np.median([c for _, c in drops.values()])))
    print('   changes by gene:', collections.Counter(l for l, _ in (mono ^ cur)).most_common(6))
