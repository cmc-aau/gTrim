"""Support for choosing --min-aligned-columns: (1) outlier rate and residual scatter by number of aligned columns in
the default run, and (2) how often sampling noise alone would produce outliers over n aligned columns.
Requires WORK with details/orig_k3_n120 (UPPER=1 KEEP=1 run_job.sh orig 3 120 25)."""
import collections, csv, math, os, statistics

base = os.path.join(os.environ['WORK'], 'details', 'orig_k3_n120')
BINS = [0, 25, 50, 75, 100, 150, 200, 300]
T = 25
CANDIDATES = 5000


def bin_label(n):
    i = max(i for i, b in enumerate(BINS) if n >= b)
    return f'{BINS[i]}-{BINS[i + 1] - 1}' if i + 1 < len(BINS) else f'{BINS[i]}+'


def binom_cdf(k, n, p):
    return sum(math.factorial(n) // (math.factorial(i) * math.factorial(n - i)) * p ** i * (1 - p) ** (n - i) for i in range(k + 1))


masked = collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(base, 'divergent_loci_masked.tsv')), delimiter='\t'):
    masked[r['masked_by']].add((r['genome_id'], r['locus']))

for test, filename, y_col, e_col, by in (('low-identity', 'divergent_loci_comparisons.tsv', 'locus_identity', 'expected_identity', 'low_identity'),
                                         ('high-identity', 'divergent_loci_upper_comparisons.tsv', 'genome_aai_loo', 'expected_aai', 'high_identity')):
    resid, outliers, min_aligned = collections.defaultdict(list), collections.Counter(), {}
    for r in csv.DictReader(open(os.path.join(base, filename)), delimiter='\t'):
        a = int(r['locus_aligned'])
        resid[bin_label(a)].append(float(r[y_col]) - float(r[e_col]))
        outliers[bin_label(a)] += r['outlier'] == 'true'
        key = (r['genome_id'], r['locus'])
        min_aligned[key] = min(min_aligned.get(key, a), a)
    test_masked = masked[by] | masked['both']
    loci = collections.Counter(bin_label(a) for a in min_aligned.values())
    loci_masked = collections.Counter(bin_label(a) for k, a in min_aligned.items() if k in test_masked)
    print(f'== {test} test: by aligned columns (genome loci binned by their fewest aligned columns)')
    print('| Aligned | Comparisons | Outlier % | Robust SD of residual | Genome loci | Masked % |')
    print('|---|---|---|---|---|---|')
    for b in [bin_label(x) for x in BINS]:
        if not resid[b]:
            continue
        m = statistics.median(resid[b])
        rsd = 1.4826 * statistics.median([abs(x - m) for x in resid[b]])
        print(f'| {b} | {len(resid[b]):,} | {100 * outliers[b] / len(resid[b]):.2f} | {rsd:.1f} | {loci[b]:,} | {100 * loci_masked[b] / loci[b]:.2f} |')
    print()

print(f'Sampling noise alone, n aligned columns, T = {T}:')
print(f'| n | Low-identity: P(outlier), expected identity 50 / 70 / 90% | Low-identity: P(2 of 3) at 70% | '
      f'High-identity: chance matches >= 90% identity among {CANDIDATES:,} genomes at 50 / 60 / 70% identity |')
print('|---|---|---|---|')
for n in (10, 20, 30, 40, 50, 75, 100):
    low = [binom_cdf(math.ceil(n * (p - T / 100)) - 1, n, p) for p in (0.5, 0.7, 0.9)]
    high = [CANDIDATES * (1 - binom_cdf(math.ceil(0.9 * n) - 1, n, q)) for q in (0.5, 0.6, 0.7)]
    print(f'| {n} | {" / ".join(f"{x:.1e}" for x in low)} | {3 * low[1] ** 2:.0e} | {" / ".join(f"{max(x, 0):.2g}" for x in high)} |')
