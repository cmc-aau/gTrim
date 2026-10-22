"""Summarise sensitivity sweep 1 of mask-divergent-loci (AAI mode): k in {1,3,5,7} x window size. Requires WORK."""
import collections, csv, itertools, os, re, sys

base = os.environ["WORK"]
KS = [1, 3, 5, 7]
NS = [20, 40, 80, 120, 200, 400]
REF = (3, 120)
CONFIGS = [(k, n) for k in KS for n in NS]

def masked_set(ds, k, n):
    path = os.path.join(base, 'runs', f'{ds}_k{k}_n{n}', 'divergent_loci_masked.tsv')
    with open(path) as f:
        return {(r['locus'], r['genome_id']) for r in csv.DictReader(f, delimiter='\t')}

def outlier_cmps(ds, k, n):
    log = open(os.path.join(base, 'runs', f'{ds}_k{k}_n{n}', 'gtrim.log')).read()
    return int(re.search(r'Identified (\d+) outlier comparisons', log).group(1))

def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0

def table(header, rows):
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
    line = lambda r: '  '.join(str(x).rjust(w) for x, w in zip(r, widths))
    print(line(header)); print('  '.join('-' * w for w in widths))
    for r in rows: print(line(r))
    print()

orig = {c: masked_set('orig', *c) for c in CONFIGS}
planted = {c: masked_set('planted', *c) for c in CONFIGS}
plants = list(csv.DictReader(open(os.path.join(base, 'planted_loci.tsv')), delimiter='\t'))
plant_keys = {(p['locus'], p['genome_id']) for p in plants}
bands = sorted({p['aai_band'] for p in plants})

# consensus: genome loci masked in at least half the configurations on the original data
freq = collections.Counter(x for s in orig.values() for x in s)
consensus = {x for x, c in freq.items() if c >= len(CONFIGS) / 2}

print('== A. Original data: masked genome loci by configuration')
rows = []
for k, n in CONFIGS:
    s = orig[(k, n)]
    per_gene = collections.Counter(l for l, _ in s)
    others = [orig[c] for c in CONFIGS if c != (k, n)]
    rows.append([k, n, outlier_cmps('orig', k, n), len(s), len(per_gene),
                 sum(1 for v in per_gene.values() if v >= 25),
                 f'{jaccard(s, orig[REF]):.3f}',
                 f'{sum(jaccard(s, o) for o in others) / len(others):.3f}',
                 f'{len(s & consensus) / len(consensus):.3f}',
                 f'{len(s & consensus) / len(s):.3f}' if s else 'na'])
table(['k', 'N', 'outlier_cmps', 'masked', 'genes>=1', 'genes>=25', 'J_vs_k3N120', 'mean_J', 'consensus_recall', 'consensus_prec'], rows)
print(f'Consensus set (masked in >= {len(CONFIGS) // 2} of {len(CONFIGS)} configurations): {len(consensus)} genome loci; '
      f'masked in all configurations: {sum(1 for c in freq.values() if c == len(CONFIGS))}; union: {len(freq)}\n')

print('== B. How many configurations mask each genome locus (original data)')
hist = collections.Counter(freq.values())
bins = [(1, 1), (2, 6), (7, 12), (13, 18), (19, 23), (24, 24)]
table(['configs', 'genome_loci'], [[f'{a}-{b}' if a != b else str(a), sum(hist[i] for i in range(a, b + 1))] for a, b in bins])

print('== C. Jaccard similarity between window sizes (original data), per k')
for k in KS:
    rows = [[n] + [f'{jaccard(orig[(k, n)], orig[(k, m)]):.2f}' for m in NS] for n in NS]
    print(f'k = {k}'); table(['N'] + [str(m) for m in NS], rows)

print('== D. Jaccard similarity between k values (original data), per window size')
rows = [[n] + [f'{jaccard(orig[(a, n)], orig[(b, n)]):.2f}' for a, b in itertools.combinations(KS, 2)] for n in NS]
table(['N'] + [f'k{a}v{b}' for a, b in itertools.combinations(KS, 2)], rows)

print('== E. Genes with the most masked genomes (original data, k=3 and N=120), and their rank in other configurations')
ref_counts = collections.Counter(l for l, _ in orig[REF])
top = [l for l, _ in ref_counts.most_common(10)]
rows = []
for l in top:
    row = [l, ref_counts[l]]
    for c in [(1, 120), (5, 120), (7, 120), (3, 20), (3, 400)]:
        row.append(collections.Counter(x for x, _ in orig[c])[l])
    rows.append(row)
table(['gene', 'k3N120', 'k1N120', 'k5N120', 'k7N120', 'k3N20', 'k3N400'], rows)

print('== F. Planted divergent loci: recall by donor AAI band, and background (non-planted) masks')
rows = []
for k, n in CONFIGS:
    s = planted[(k, n)]
    rec = []
    for b in bands:
        keys = {(p['locus'], p['genome_id']) for p in plants if p['aai_band'] == b}
        rec.append(f'{len(s & keys) / len(keys):.2f}')
    rows.append([k, n] + rec + [f'{len(s & plant_keys) / len(plant_keys):.2f}', len(s - plant_keys), len(orig[(k, n)])])
table(['k', 'N'] + [f'recall_{b}' for b in bands] + ['recall_all', 'background', 'orig_masked'], rows)
