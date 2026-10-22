"""Summarise the sweep of high-identity test parameters (k', N', T') with the low-identity test fixed at
k=3, N=120, T=25. Uses $WORK/runs/{control,planted}_nodonor_uk{k}_un{n}[_ut{t}] (run_job.sh with
UPPER_SWEEP=1), $WORK/planted_nodonor_loci.tsv, and the comparison tables in
$WORK/details/planted_nodonor_k3_n120 (for each plant's best remaining locus identity). Requires WORK."""
import collections, csv, os

base = os.environ['WORK']
KS = [1, 2, 3, 4, 5, 7]
NS = [20, 40, 80, 120, 200, 400]
TS = [10, 15, 20, 25, 30, 35, 40]
DEFAULT = (3, 120, 25)


def run_dir(ds, k, n, t):
    d = os.path.join(base, 'runs', f'{ds}_nodonor_uk{k}_un{n}')
    return d if t == 25 else f'{d}_ut{t}'


def load(ds, k, n, t):
    """Return masked sets for the high-identity test, the low-identity test, and either test."""
    sets = {'high': set(), 'low': set(), 'either': set()}
    for r in csv.DictReader(open(os.path.join(run_dir(ds, k, n, t), 'divergent_loci_masked.tsv')), delimiter='\t'):
        key = (r['locus'], r['genome_id'])
        sets['either'].add(key)
        if r['masked_by'] in ('high_identity', 'both'):
            sets['high'].add(key)
        if r['masked_by'] in ('low_identity', 'both'):
            sets['low'].add(key)
    return sets


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0


def table(title, header, rows):
    print(title)
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
    line = lambda r: '  '.join(str(x).rjust(w) for x, w in zip(r, widths))
    print(line(header))
    print('  '.join('-' * w for w in widths))
    for r in rows:
        print(line(r))
    print()


plants = list(csv.DictReader(open(os.path.join(base, 'planted_nodonor_loci.tsv')), delimiter='\t'))
plant_keys = {(p['locus'], p['genome_id']) for p in plants}
band_keys = {b: {(p['locus'], p['genome_id']) for p in plants if p['aai_band'] == b} for b in sorted({p['aai_band'] for p in plants})}

# best locus identity of each planted sequence to a genome remaining in the data
top_identity = {}
for r in csv.DictReader(open(os.path.join(base, 'details', 'planted_nodonor_k3_n120', 'divergent_loci_upper_comparisons.tsv')), delimiter='\t'):
    key = (r['locus'], r['genome_id'])
    if key in plant_keys and r['neighbour_rank'] == '1':
        top_identity[key] = float(r['locus_identity'])
ID_BINS = [(0, 80), (80, 90), (90, 95), (95, 98), (98, 100.01)]
id_keys = {f'{lo}-{min(hi, 100):.0f}': {k for k, v in top_identity.items() if lo <= v < hi} for lo, hi in ID_BINS}

configs = sorted({(k, n, 25) for k in KS for n in NS} | {(k, 120, t) for k in KS for t in TS})
control = {c: load('control', *c) for c in configs}
planted = {c: load('planted', *c) for c in configs}

def metrics(c):
    ch, ph = control[c]['high'], planted[c]['high']
    return {
        'control_high': len(ch),
        'new_vs_low': len(ch - control[c]['low']),
        'J_default': jaccard(ch, control[DEFAULT]['high']),
        'recall_high': len(ph & plant_keys) / len(plant_keys),
        'recall_either': len(planted[c]['either'] & plant_keys) / len(plant_keys),
        'knock_on': len(ph - plant_keys) - len(ch),
    }


m = {c: metrics(c) for c in configs}

print(f"Low-identity test fixed at k=3, N=120, T=25: {len(control[DEFAULT]['low'])} masks on the control; "
      f"recall of planted loci {len(planted[DEFAULT]['low'] & plant_keys) / len(plant_keys):.2f}\n")

grid = lambda key, fmt: [[k] + [fmt(m[(k, n, 25)][key]) for n in NS] for k in KS]
table("== A. High-identity masks on the control data (T'=25)", ["k'"] + [f"N'={n}" for n in NS], grid('control_high', str))
table("== B. Recall of planted loci by the high-identity test (T'=25)", ["k'"] + [f"N'={n}" for n in NS], grid('recall_high', lambda v: f'{v:.2f}'))
table("== C. Knock-on high-identity masks caused by planting (T'=25)", ["k'"] + [f"N'={n}" for n in NS], grid('knock_on', str))

tgrid = lambda key, fmt: [[t] + [fmt(m[(k, 120, t)][key]) for k in KS] for t in TS]
table("== D. High-identity masks on the control data (N'=120)", ["T'"] + [f"k'={k}" for k in KS], tgrid('control_high', str))
table("== D2. ... of which not masked by the low-identity test (N'=120)", ["T'"] + [f"k'={k}" for k in KS], tgrid('new_vs_low', str))
table("== E. Recall of planted loci by the high-identity test (N'=120)", ["T'"] + [f"k'={k}" for k in KS], tgrid('recall_high', lambda v: f'{v:.2f}'))
table("== F. Recall of planted loci by either test (N'=120; low-identity test fixed)", ["T'"] + [f"k'={k}" for k in KS], tgrid('recall_either', lambda v: f'{v:.2f}'))
table("== G. Knock-on high-identity masks caused by planting (N'=120)", ["T'"] + [f"k'={k}" for k in KS], tgrid('knock_on', str))
table("== H. Jaccard of high-identity masks on the control vs the default (k'=3, N'=120, T'=25)", ["T'"] + [f"k'={k}" for k in KS], tgrid('J_default', lambda v: f'{v:.2f}'))

for k in [2, 3, 5]:
    rows = []
    for t in TS:
        ph = planted[(k, 120, t)]['high']
        rows.append([t] + [f'{len(ph & keys) / len(keys):.2f}' if keys else 'na' for keys in band_keys.values()]
                    + [f'{len(ph & keys) / len(keys):.2f}' if keys else 'na' for keys in id_keys.values()])
    table(f"== I. High-identity recall by donor AAI band and by best remaining locus identity (k'={k}, N'=120)",
          ["T'"] + [f'band {b}' for b in band_keys] + [f'id {b} (n={len(v)})' for b, v in id_keys.items()], rows)

print('== J. Genes with the most high-identity masks on the control, by configuration (N\'=120)')
for c in [(3, 120, 25), (3, 120, 20), (5, 120, 20), (5, 120, 25)]:
    counts = collections.Counter(l for l, _ in control[c]['high']).most_common(5)
    print(f"  k'={c[0]} T'={c[2]}: " + ', '.join(f'{l} {n}' for l, n in counts))
