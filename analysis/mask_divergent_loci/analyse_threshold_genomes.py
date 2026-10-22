"""Per-genome effect of lowering both outlier thresholds (T and T') from 25 to 20 on the original data:
genomes gaining the most masked loci, and genomes that drop below 30% of all loci after masking.
Requires WORK with details/orig_k3_n120 (default run with both tests; UPPER=1 KEEP=1 run_job.sh orig 3 120 25).
Masks at other thresholds are recomputed from the comparisons tables, since expected values don't depend on T."""
import collections, csv, os

base = os.path.join(os.environ['WORK'], 'details', 'orig_k3_n120')
THRESHOLDS = (25, 20)
MIN_FRAC = 0.3
TOP = 20


def load(filename, y_col, expected_col):
    by_key = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(base, filename)), delimiter='\t'):
        by_key[(r['genome_id'], r['locus'])].append((float(r[y_col]), float(r[expected_col])))
    return by_key


def masked(by_key, t):
    return {k for k, v in by_key.items() if 2 * sum(y < e - t for y, e in v) > len(v)}


low = load('divergent_loci_comparisons.tsv', 'locus_identity', 'expected_identity')
high = load('divergent_loci_upper_comparisons.tsv', 'genome_aai_loo', 'expected_aai')
genome_rows = list(csv.DictReader(open(os.path.join(base, 'divergent_loci_genome_summary.tsv')), delimiter='\t'))
present = {r['genome_id']: int(r['num_loci']) for r in genome_rows}
num_loci = len({locus for _, locus in low})

per_genome = {}
for t in THRESHOLDS:
    lo, hi = masked(low, t), masked(high, t)
    print(f"T = T' = {t}: low-identity {len(lo):,}, high-identity {len(hi):,}, either {len(lo | hi):,} masked genome loci")
    per_genome[t] = {name: collections.Counter(g for g, _ in s) for name, s in (('low', lo), ('high', hi), ('either', lo | hi))}

a, b = THRESHOLDS
added = {g: per_genome[b]['either'][g] - per_genome[a]['either'][g] for g in present}
dist = collections.Counter(min(v, 6) for v in added.values())
print(f'\nGenomes gaining masked loci: {sum(v > 0 for v in added.values()):,} of {len(present):,}; '
      f'gained loci (6 = 6+): {sorted(dist.items())}')

print(f'\nTop {TOP} genomes by additional masked loci, T = T\' = {a} -> {b} ({num_loci} loci in total):')
print(f'| Genome | Loci present | Masked T={a} | Masked T={b} | Added | Added by low-identity | Added by high-identity | Loci left T={a} | Loci left T={b} |')
print('|---|---|---|---|---|---|---|---|---|')
for g in sorted(present, key=lambda g: (-added[g], -per_genome[b]['either'][g], g))[:TOP]:
    ma, mb = per_genome[a]['either'][g], per_genome[b]['either'][g]
    dl = per_genome[b]['low'][g] - per_genome[a]['low'][g]
    dh = per_genome[b]['high'][g] - per_genome[a]['high'][g]
    left = lambda m: f'{present[g] - m} ({100 * (present[g] - m) / num_loci:.0f}%)'
    print(f'| {g} | {present[g]} | {ma} | {mb} | +{added[g]} | +{dl} | +{dh} | {left(ma)} | {left(mb)} |')

min_loci = MIN_FRAC * num_loci
print(f'\nGenomes with >= {MIN_FRAC:.0%} of loci ({min_loci:.1f}) before masking: {sum(n >= min_loci for n in present.values()):,}')
for t in THRESHOLDS:
    drop = sorted(g for g, n in present.items() if n >= min_loci and n - per_genome[t]['either'][g] < min_loci)
    print(f"T = T' = {t}: {len(drop)} drop below {MIN_FRAC:.0%}: " + ', '.join(f"{g} ({present[g]} -> {present[g] - per_genome[t]['either'][g]})" for g in drop))
