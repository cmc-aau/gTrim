"""Look for divergent copies shared by a group of close relatives, which the k = 3 tests can't see: a genome locus
whose 3 nearest neighbours agree with it (at most 1 outlier), while at least 4 of its next 7 neighbours are
outliers. Candidates at the same locus that are among each other's 3 nearest neighbours are grouped.
Requires WORK with details/orig_k10_n120 (UPPER=1 KEEP=1 run_job.sh orig 10 120 25) and details/orig_k3_n120.
Only groups of 3-10 genomes are visible this way."""
import collections, csv, os, statistics

work = os.environ['WORK']
K10 = os.path.join(work, 'details', 'orig_k10_n120')
K3 = os.path.join(work, 'details', 'orig_k3_n120')
T = 25
NEAR = 3
MIN_FAR_OUTLIERS = 4
CLEAN_DEFICIT = 10


def load(filename, y_col, expected_col):
    by_key = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(K10, filename)), delimiter='\t'):
        y, e = float(r[y_col]), float(r[expected_col])
        by_key[(r['locus'], r['genome_id'])].append((int(r['neighbour_rank']), r['neighbour_id'], y < e - T, e - y, float(r['genome_aai_loo'])))
    for v in by_key.values():
        v.sort()
    return by_key


def groups_of(cand):
    """Connected groups of candidates at each locus, linked through their nearest neighbours."""
    by_locus = collections.defaultdict(set)
    for locus, g in cand:
        by_locus[locus].add(g)

    groups = []
    for locus, genomes in by_locus.items():
        adj = collections.defaultdict(set)
        for g in genomes:
            for _, n, *_ in cand[(locus, g)][:NEAR]:
                if n in genomes:
                    adj[g].add(n)
                    adj[n].add(g)
        seen = set()
        for g in sorted(genomes):
            if g in seen:
                continue
            stack, comp = [g], []
            seen.add(g)
            while stack:
                x = stack.pop()
                comp.append(x)
                for y in adj[x] - seen:
                    seen.add(y)
                    stack.append(y)
            groups.append((locus, sorted(comp)))
    return groups


masked_k3 = {(r['locus'], r['genome_id']) for r in csv.DictReader(open(os.path.join(K3, 'divergent_loci_masked.tsv')), delimiter='\t')}
total = sum(int(r['num_loci']) for r in csv.DictReader(open(os.path.join(K3, 'divergent_loci_genome_summary.tsv')), delimiter='\t'))
for test, filename, y_col, e_col in (('low-identity', 'divergent_loci_comparisons.tsv', 'locus_identity', 'expected_identity'),
                                     ('high-identity', 'divergent_loci_upper_comparisons.tsv', 'genome_aai_loo', 'expected_aai')):
    by_key = load(filename, y_col, e_col)
    n_full = sum(len(v) == 10 for v in by_key.values())
    cand = {k: v for k, v in by_key.items()
            if len(v) == 10 and sum(o for _, _, o, _, _ in v[:NEAR]) <= 1 and sum(o for _, _, o, _, _ in v[NEAR:]) >= MIN_FAR_OUTLIERS}
    groups = groups_of(cand)
    big = [(locus, c) for locus, c in groups if len(c) >= 3]
    keys = [(locus, g) for locus, c in big for g in c]
    clean = [(locus, c) for locus, c in big
             if all(statistics.mean(d for *_, d, _ in cand[(locus, g)][:NEAR]) < CLEAN_DEFICIT for g in c)]
    clean_keys = [(locus, g) for locus, c in clean for g in c]

    print(f'== {test} test, k = 10, T = {T}: {n_full:,} genome loci with 10 comparisons; {len(cand):,} candidates in {len(groups):,} groups')
    print('   group sizes (6 = 6+):', sorted(collections.Counter(min(len(c), 6) for _, c in groups).items()))
    print(f'   groups of >= 3: {len(big)} groups, {len(keys)} genome loci ({100 * len(keys) / total:.2f}% of genome loci), '
          f'{len({g for _, g in keys})} genomes, {len({l for l, _ in big})} loci; {sum(k in masked_k3 for k in keys)} already masked at k = 3')
    print('   groups of >= 3 by locus:', collections.Counter(l for l, _ in big).most_common(8))
    if clean:
        print(f'   clean groups (every member: mean deficit to its {NEAR} nearest < {CLEAN_DEFICIT}): {len(clean)} groups, '
              f'{len(clean_keys)} genome loci, {len({l for l, _ in clean})} loci; median deficit to neighbours 4-10: '
              f'{statistics.median(x[3] for k in clean_keys for x in cand[k][NEAR:]):.1f}')
        print('   clean groups by locus:', collections.Counter(l for l, _ in clean).most_common(8))
    for locus, c in sorted(big, key=lambda x: (-len(x[1]), x[0]))[:6]:
        v = cand[(locus, c[0])]
        print(f'   {locus} n={len(c)} e.g. {c[0]}: neighbours 1-3 AAI {[round(x[4], 1) for x in v[:NEAR]]} deficit {[round(x[3], 1) for x in v[:NEAR]]}; '
              f'4-10 AAI {[round(x[4], 1) for x in v[NEAR:]]} deficit {[round(x[3], 1) for x in v[NEAR:]]}')
