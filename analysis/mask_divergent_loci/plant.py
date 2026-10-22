"""Plant divergent loci into a copy of the per-locus MSAs.

For each of the AAI bands, picks genomes at random (one plant per genome) and a random locus
the genome has, and replaces its sequence with that of a donor genome whose pooled AAI to the
focal genome falls within the band. Writes the planted MSAs and a TSV of planted loci.

Usage: plant.py <src_dir> <out_dir> <seed> [control_dir]
If control_dir is given, donor genomes are removed from the planted MSAs, and the original MSAs
with donor genomes removed are written to control_dir (for measuring knock-on masks).
"""
import gzip, os, random, sys
import numpy as np

src_dir, out_dir, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
control_dir = sys.argv[4] if len(sys.argv) > 4 else None
BANDS = [(40, 55), (55, 70), (70, 85)]
PER_BAND = 100
AA = np.zeros(256, dtype=bool)
for c in b'ACDEFGHIKLMNPQRSTVWY':
    AA[c] = True

rng = random.Random(seed)
files = sorted(f for f in os.listdir(src_dir) if f.endswith('.faa.gz'))
loci = []  # (name, order, seqs)
for f in files:
    order, seqs, k = [], {}, None
    for line in gzip.open(os.path.join(src_dir, f), 'rt'):
        line = line.strip()
        if line.startswith('>'):
            k = line[1:].split()[0]; order.append(k); seqs[k] = []
        else:
            seqs[k].append(line)
    loci.append((f, order, {k: ''.join(v) for k, v in seqs.items()}))

genomes = sorted(set(g for _, o, _ in loci for g in o))
gidx = {g: i for i, g in enumerate(genomes)}
lens = [len(next(iter(s.values()))) for _, _, s in loci]
offsets = np.cumsum([0] + lens)
mat = np.zeros((len(genomes), offsets[-1]), dtype=np.uint8)
for li, (_, order, seqs) in enumerate(loci):
    for g in order:
        mat[gidx[g], offsets[li]:offsets[li + 1]] = np.frombuffer(seqs[g].upper().encode(), dtype=np.uint8)
valid = AA[mat]

def aai(a, b):
    both = valid[a] & valid[b]
    n = both.sum()
    return 100.0 * (both & (mat[a] == mat[b])).sum() / n if n > 0 else None

used = set()
plants = []  # (locus file, focal, donor, band, donor_aai)
focal_order = genomes[:]
rng.shuffle(focal_order)
fi = 0
for lo, hi in BANDS:
    count = 0
    while count < PER_BAND:
        focal = focal_order[fi]; fi += 1
        locus_choices = [li for li, (_, _, s) in enumerate(loci) if focal in s]
        li = rng.choice(locus_choices)
        cands = [g for g in loci[li][1] if g != focal and g not in used]
        rng.shuffle(cands)
        for donor in cands[:3000]:
            v = aai(gidx[focal], gidx[donor])
            if v is not None and lo <= v < hi:
                plants.append((loci[li][0], focal, donor, f'{lo}-{hi}', v))
                used.update([focal, donor])
                count += 1
                break

os.makedirs(out_dir, exist_ok=True)
planted = {(p[0], p[1]): p[2] for p in plants}
removed = {p[2] for p in plants} if control_dir else set()
for f, order, seqs in loci:
    with gzip.open(os.path.join(out_dir, f), 'wt', compresslevel=1) as out:
        for g in order:
            if g in removed:
                continue
            donor = planted.get((f, g))
            out.write(f'>{g}\n{seqs[donor] if donor else seqs[g]}\n')

if control_dir:
    os.makedirs(control_dir, exist_ok=True)
    for f, order, seqs in loci:
        with gzip.open(os.path.join(control_dir, f), 'wt', compresslevel=1) as out:
            for g in order:
                if g not in removed:
                    out.write(f'>{g}\n{seqs[g]}\n')

with open(os.path.join(os.path.dirname(out_dir.rstrip('/')), os.path.basename(out_dir.rstrip('/')).replace('_msas', '') + '_loci.tsv'), 'w') as out:
    out.write('locus\tgenome_id\tdonor_id\taai_band\tdonor_aai\n')
    for f, g, d, band, v in plants:
        out.write(f'{f.split(".")[0]}\t{g}\t{d}\t{band}\t{v:.2f}\n')
print(f'Planted {len(plants)} loci across {len(BANDS)} AAI bands into {out_dir}')
