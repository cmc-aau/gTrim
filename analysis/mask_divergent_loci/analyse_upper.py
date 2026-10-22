"""Evaluate the high-identity (upper) test on planted loci with donor genomes removed. Requires WORK with
planted_nodonor_loci.tsv and details/{control,planted}_nodonor_k3_n120 (default k=3, N=120, T=25)."""
import collections, csv, os, statistics

base = os.environ['WORK']
TESTS = ['low_identity', 'high_identity', 'either']


def load_masked(name):
    by_test = {t: set() for t in TESTS}
    for r in csv.DictReader(open(os.path.join(base, 'details', name, 'divergent_loci_masked.tsv')), delimiter='\t'):
        key = (r['locus'], r['genome_id'])
        by_test['either'].add(key)
        if r['masked_by'] in ('low_identity', 'both'):
            by_test['low_identity'].add(key)
        if r['masked_by'] in ('high_identity', 'both'):
            by_test['high_identity'].add(key)
    return by_test


def table(header, rows):
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
    line = lambda r: '  '.join(str(x).rjust(w) for x, w in zip(r, widths))
    print(line(header))
    print('  '.join('-' * w for w in widths))
    for r in rows:
        print(line(r))
    print()


plants = list(csv.DictReader(open(os.path.join(base, 'planted_nodonor_loci.tsv')), delimiter='\t'))
plant_keys = {(p['locus'], p['genome_id']) for p in plants}
control = load_masked('control_nodonor_k3_n120')
planted = load_masked('planted_nodonor_k3_n120')

# per plant: mean neighbour LOO AAI (low test) and best locus identity of the planted sequence (high test)
nb_aai = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(base, 'details', 'planted_nodonor_k3_n120', 'divergent_loci_comparisons.tsv')), delimiter='\t'):
    key = (r['locus'], r['genome_id'])
    if key in plant_keys:
        nb_aai[key].append(float(r['genome_aai_loo']))
top_identity = {}
for r in csv.DictReader(open(os.path.join(base, 'details', 'planted_nodonor_k3_n120', 'divergent_loci_upper_comparisons.tsv')), delimiter='\t'):
    key = (r['locus'], r['genome_id'])
    if key in plant_keys and r['neighbour_rank'] == '1':
        top_identity[key] = float(r['locus_identity'])

gap = {}
for p in plants:
    key = (p['locus'], p['genome_id'])
    gap[key] = statistics.mean(nb_aai[key]) - float(p['donor_aai']) if nb_aai[key] else float('nan')


def recall(keys, test):
    return f'{len(planted[test] & keys) / len(keys):.2f}' if keys else 'na'


print(f'Plants: {len(plants)} (donor genomes removed from the data)\n')

print('== 1. Recall by donor AAI band and test')
bands = sorted({p['aai_band'] for p in plants})
rows = []
for b in bands + ['all']:
    keys = {(p['locus'], p['genome_id']) for p in plants if b == 'all' or p['aai_band'] == b}
    high_only = len(keys & (planted['high_identity'] - planted['low_identity'])) / len(keys)
    rows.append([b, len(keys)] + [recall(keys, t) for t in TESTS] + [f'{high_only:.2f}'])
table(['band', 'n'] + TESTS + ['high_only'], rows)
print('high_only = fraction of plants in the band found only by the high-identity test\n')

print('== 2. Recall by gap (neighbour AAI - donor AAI)')
GAP_BINS = [(-100, 15), (15, 20), (20, 25), (25, 30), (30, 40), (40, 101)]
bin_name = lambda lo, hi: f'<{hi}' if lo < 0 else (f'{lo}+' if hi > 100 else f'{lo}-{hi}')
rows = []
for lo, hi in GAP_BINS:
    keys = {k for k, v in gap.items() if lo <= v < hi}
    rows.append([bin_name(lo, hi), len(keys)] + [recall(keys, t) for t in TESTS])
table(['gap', 'n'] + TESTS, rows)

print('== 3. Recall of the high-identity test by best locus identity of the planted sequence to a remaining genome')
print('     (how closely the donor lineage is still represented after removing the donor)')
ID_BINS = [(0, 80), (80, 90), (90, 95), (95, 98), (98, 100.01)]
rows = []
for lo, hi in ID_BINS:
    keys = {k for k, v in top_identity.items() if lo <= v < hi}
    rows.append([f'{lo}-{min(hi, 100):.0f}', len(keys)] + [recall(keys, t) for t in TESTS])
table(['best_identity', 'n'] + TESTS, rows)

print('== 4. Masks on the control (donors removed, no plants) and knock-on masks caused by planting')
rows = []
for t in TESTS:
    rows.append([t, len(control[t]), len(planted[t] - plant_keys), len(planted[t] - plant_keys) - len(control[t])])
table(['test', 'control_masked', 'planted_non_plant_masked', 'knock_on'], rows)
