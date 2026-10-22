"""Explain k=1 masks relative to k=3 (N=120, threshold 25): mutual nearest neighbours, knock-on masks,
and masks missed at k=1. Requires WORK with details/{orig,planted}_k{1,3}_n120 (run_job.sh with KEEP=1)."""
import collections, csv, os

base = os.environ['WORK']


def load(name):
    d = os.path.join(base, 'details', name)
    rows = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(d, 'divergent_loci_comparisons.tsv')), delimiter='\t'):
        rows[(r['locus'], r['genome_id'])].append((r['neighbour_id'], r['outlier'] == 'true', float(r['identity_deficit'])))
    masked = {(r['locus'], r['genome_id']) for r in csv.DictReader(open(os.path.join(d, 'divergent_loci_masked.tsv')), delimiter='\t')}
    return rows, masked


k1, m1 = load('orig_k1_n120')
k3, m3 = load('orig_k3_n120')
nn = {key: v[0][0] for key, v in k1.items()}

print(f'k=1 masked: {len(m1)}   k=3 masked: {len(m3)}   both: {len(m1 & m3)}   k1 only: {len(m1 - m3)}   k3 only: {len(m3 - m1)}\n')

mutual = other = notmasked = 0
for (l, a) in m1:
    b = nn[(l, a)]
    if (l, b) in m1:
        if nn.get((l, b)) == a:
            mutual += 1
        else:
            other += 1
    else:
        notmasked += 1
print('k=1 masked genome A with nearest neighbour B:')
print(f'  B also masked, mutual nearest neighbours: {mutual}')
print(f'  B also masked, via a different nearest neighbour: {other}')
print(f'  B not masked: {notmasked}\n')

indeg = collections.Counter((l, b) for (l, a), b in nn.items())
print(f'Mean number of genomes using a genome locus as nearest neighbour: {sum(indeg.values()) / len(k1):.2f}')
print(f'  k=3-masked genome loci: {sum(indeg[x] for x in m3) / len(m3):.2f}; nobody\'s nearest neighbour: {sum(1 for x in m3 if indeg[x] == 0) / len(m3):.2f}\n')

coll = sum(1 for (l, a) in m1 - m3 if (l, nn[(l, a)]) in m3)
border = [k1[x][0][2] for x in m1 - m3 if (x[0], nn[x]) not in m3]
print(f'k=1-only masks: {len(m1 - m3)}')
print(f'  knock-on (nearest neighbour masked at k=3): {coll}')
print(f'  other: {len(border)}; within 5 points of the threshold: {sum(1 for d in border if d < 30)}\n')

print(f'k=3-only masks: {len(m3 - m1)}; nearest-neighbour comparison not an outlier: {sum(1 for x in m3 - m1 if not k1[x][0][1])}')
print(f'  nearest neighbour also masked at k=3 (shared divergent copy): {sum(1 for x in m3 - m1 if (x[0], nn[x]) in m3)}\n')

p1, pm1 = load('planted_k1_n120')
pnn = {key: v[0][0] for key, v in p1.items()}
plants = {(p['locus'], p['genome_id']) for p in csv.DictReader(open(os.path.join(base, 'planted_loci.tsv')), delimiter='\t')}
det = plants & pm1
coll = collections.Counter()
for (l, a) in pm1 - plants:
    b = pnn[(l, a)]
    if (l, b) in det:
        coll[(l, b)] += 1
dist = collections.Counter(coll.get(p, 0) for p in det)
print(f'Planted data, k=1: {len(det)} of {len(plants)} plants detected; knock-on masks: {sum(coll.values())}')
print(f'  knock-on masks per detected plant: {dict(sorted(dist.items()))}')
