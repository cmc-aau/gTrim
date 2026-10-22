# mask-divergent-loci: methodology and evaluation

This document describes how `gtrim mask-divergent-loci` identifies and masks divergent loci, and records how well the method works and how sensitive it is to its parameters. Update it whenever new tests are run.

- **Status:** under evaluation.
  - **Low-identity test:** tested on the example data with planted divergent loci, including sweeps over k, N and T (sections 4.2–4.5).
  - **High-identity test:** added 2026-10-22; evaluated with planted loci whose donor genomes were removed (section 4.8).
  - **Bucketed evaluation:** Experiment 1 (low-identity test, section 4.10) and Experiment 2 (high-identity test, section 4.11) done.
  - **Review follow-ups:** tie handling in the running median fixed, suspect-genome flag added, per-genome effect of T = 20 and a search for divergent copies shared by close relatives (section 4.12); minimum aligned columns added (section 4.13). T = T′ = 25 kept as the default (conclusion 3).
- **Non-decreasing trend:** default since 2026-10-23 (`--method monotone-median`, section 4.16); earlier results used `median`.
- **Last updated:** 2026-10-23 (gTrim 0.1.0, unreleased; working tree after commit cf0ec4e).

## 1. Purpose

Phylogenomic analyses based on concatenated marker genes assume that each marker in a genome reflects that genome's evolutionary history. A gene can break this assumption when it is a paralog, a contaminant, a horizontally transferred copy, or misaligned. Such a gene is far less similar to the same gene in the genome's closest relatives than the genomes are to each other overall. `mask-divergent-loci` finds these cases and replaces the affected sequence with gaps, so the rest of the genome's data can still be used.

## 2. Method

### 2.1 Input

- A directory of per-locus MSAs (`.faa`, `.fna` or `.fasta`, optionally gzipped), one file per gene.
- In each FASTA header, the first whitespace-delimited token is the genome ID, e.g. `GCA_003055145.1`.
- The sequence type is detected automatically. Only the 20 standard amino acids (or A, C, G, T) count as residues. Gaps and ambiguous characters (e.g. `X`, `*`) are treated as gaps.

### 2.2 Pairwise identity

For two aligned sequences, a column is **aligned** if both sequences have a residue there. **Identity** = identical columns / aligned columns.

### 2.3 Genome similarity and neighbour selection

Each locus *L* of each focal genome *g* is evaluated separately:

1. **Leave-one-out AAI:** for every other genome *h*, the AAI between *g* and *h* is computed by pooling identical and aligned counts across all loci **except L**. The locus being evaluated therefore has no influence on which neighbours are chosen or on the x-value it is judged against.
2. **Candidate filters:** *h* is a candidate only if:
   - *h* has a sequence for *L*;
   - the locus overlap is sufficient: aligned(*L*) ≥ *f* × the smaller of the two genomes' residue counts at *L*, and aligned(*L*) ≥ *m* columns;
   - the genome-wide overlap is sufficient: aligned(all loci except *L*) ≥ *f* × the smaller of the two genomes' residue counts outside *L*.

   *f* is `--min-aligned-frac`, default 0.5. *m* is `--min-aligned-columns`, default 50 (section 4.13); identity over fewer columns is too noisy to test. A genome locus with no candidates is not evaluated and is never masked.
3. **Neighbours:** the *k* candidates with the highest leave-one-out AAI (`--neighbours`). Ties are broken by genome ID order.

Each (focal genome, neighbour, locus) **comparison** gives a point with x = leave-one-out AAI (%) and y = locus identity (%).

Genome similarity is always AAI computed from the input MSAs. An earlier ANI mode, which read skani output, was removed (see the change log): it can't leave the evaluated locus out, skani omits distant pairs, and the high-identity test needs AAI between all genomes anyway.

### 2.4 Expected identity and outliers

- **Trend line:** for each locus, points are sorted by x. The expected identity of a point is the **running median** of y over the *N* points nearest to it in x, including itself (`--window-size`). If a locus has fewer than *N* points, all of them are used. Points with the same x are equally near, so the window is extended to include every point tied with the x values at its ends: tied points get the same expected value, whatever their order (section 4.12).
- **Non-decreasing fit (default, `--method monotone-median`):** expected identity should not fall as AAI rises, so the running medians are replaced by their least-squares non-decreasing fit in x (pool adjacent violators; points with the same x are pooled first, so they keep one value). Local dips in the running median are treated as noise (section 4.16). `--method median` uses the running median alone.
- **Outlier:** a comparison is an outlier if y < expected − *T*, where *T* is `--outlier-threshold` in absolute percentage points (default 25). For example, with an expected identity of 50%, a locus identity below 25% is an outlier.

This **low-identity test** is one-sided: identity higher than expected is never flagged by it. Unexpectedly high identity is handled by a separate test (section 2.5).

### 2.5 High-identity test (on by default)

A locus copied from another lineage (contamination or lateral transfer) can be *more* similar to that lineage than the genome's relatedness predicts. A locus with 99% identity to a genome it shares only 65% AAI with is an example. The low-identity test can miss these copies. They are detected by reversing the roles of the two measures:

1. **Neighbours by locus identity:** for each locus *L* of focal genome *g*, take the *k*′ candidates (`--upper-neighbours`; same filters as section 2.3) with the **highest identity at *L***. Ties are broken by the highest leave-one-out AAI, then genome ID order. Breaking ties by AAI means *g*'s own lineage wins when the locus can't tell the lineages apart.
2. **Expected AAI:** for each locus, the expected leave-one-out AAI of a comparison is the running median of AAI over the *N*′ comparisons nearest in locus identity (`--upper-window-size`). Neighbours are chosen by maximal identity, so their identities are biased upwards; using AAI as the tested value avoids that bias. With the default method the trend is made non-decreasing in identity, as in section 2.4.
3. **Outlier:** a comparison is an outlier if AAI < expected AAI − *T*′, where *T*′ is in AAI points (`--upper-outlier-threshold`).
4. **Vote:** majority vote over the *k*′ comparisons, as in section 2.6.

*k*′, *N*′ and *T*′ are set separately from the low-identity test's *k*, *N* and *T*. They default to the same values (3, 120 and 25).

If *g*'s copy at *L* came from lineage *D*, all of its closest matches at *L* are in *D* while *g*'s AAI to them is low, so all *k* comparisons are outliers. Attribution follows from this:
- **Blaming the right genome:** a member of *D* that has *g* among its own closest matches gets at most one outlier vote, so it isn't masked.
- **Limits:** the test needs the source lineage to be represented in the data. A copy whose source lineage is absent, such as a paralog or a very distant contaminant, is left to the low-identity test.

It is disabled with `--no-upper-mask`.

### 2.6 Majority vote and masking

A focal genome's locus is **masked** if **more than half** of its comparisons are outliers:

| k | Outlier comparisons needed |
|---|---|
| 1 | 1 of 1 |
| 2 | 2 of 2 |
| 3 | 2 of 3 |
| 4 | 3 of 4 |
| 5 | 3 of 5 |

If fewer than *k* neighbours pass the filters, the vote is over the comparisons available. Each test votes separately. A genome locus is **masked if either test masks it**. A masked sequence is replaced by gaps of the same length in the output MSA, so row order and alignment length don't change.

Only the focal genome is masked by its own vote. A neighbour is masked only through its own comparisons. This matters for small *k* (section 4.4).

### 2.7 Outputs

| File | Contents |
|---|---|
| `masked_msas/<locus>_masked.fasta[.gz]` | Input MSAs with masked genome loci gapped |
| `divergent_loci_comparisons.tsv` | Low-identity test, every comparison: neighbour, rank, x, locus identity, aligned columns, expected identity, deficit, outlier flag |
| `divergent_loci_upper_comparisons.tsv` | High-identity test, every comparison: neighbour, rank by locus identity, locus identity, aligned columns, leave-one-out AAI, expected AAI, AAI deficit, outlier flag |
| `divergent_loci_masked.tsv` | Masked (locus, genome) pairs with vote counts from each test, `masked_by` (`low_identity`, `high_identity` or `both`), and the source of the masked copy (section 2.8): `best_match_id`, `best_match_identity`, `best_match_aai`, `source_coherent`, `match_taxon`, `shared_taxon` |
| `divergent_loci_summary.tsv` | Per locus: genomes, genomes compared, comparisons, outlier comparisons and masked genomes for each test, total masked |
| `divergent_loci_genome_summary.tsv` | Per genome: loci present (number and % of all loci) before and after masking, loci masked by each test, `suspect` (true if at least `--suspect-masked-perc` % of the genome's loci, and at least 3 loci, are masked; e.g. a contaminated or chimeric assembly), and the main source of its high-identity masks (`main_source_id`, `main_source_masks`; section 2.8). Suspect genomes are also listed in the log |
| `plots/divergent_loci_plots{,_log}.{svg,png}` | With `--plot`, low-identity test: one panel per locus showing points, trend line, cutoff and outliers, on linear and log axes |
| `plots/divergent_loci_upper_plots{,_log}.{svg,png}` | With `--plot`, high-identity test: as above, with locus identity on x and AAI on y |

### 2.8 Source of a masked copy

For each masked genome locus, `divergent_loci_masked.tsv` reports which genomes carry the most similar copies. These come from the neighbours ranked by locus identity (section 2.5, same overlap filters). The ranking is made even with `--no-upper-mask`, in which case only the single best match is kept.

| Column | Meaning |
|---|---|
| `best_match_id` | Genome whose copy of the locus is most identical (ties: highest leave-one-out AAI) |
| `best_match_identity` | Locus identity (%) to the best match |
| `best_match_aai` | Leave-one-out AAI (%) of the masked genome to the best match. High identity with low AAI points to a foreign source in the data; low identity to every genome points to a paralog, misalignment or a source that isn't sampled |
| `source_coherent` | High-identity masks only: true if the neighbours that were outliers are all ≥ 90% genome-wide AAI to each other (one source lineage), false if they come from several lineages. NA for low-identity-only masks or a single outlier neighbour |
| `match_taxon` | With `--taxonomy`: lowest taxon shared by the matching genomes. These are the outlier neighbours for high-identity masks, otherwise the *k*′ neighbours by locus identity (the single best match with `--no-upper-mask`) |
| `shared_taxon` | With `--taxonomy`: lowest taxon shared by the masked genome and the matching genomes, i.e. how far apart the source is |

**Per genome.** `divergent_loci_genome_summary.tsv` reports the main source of a genome's high-identity masks. `main_source_id` is the most common best match among them (ties: first in genome ID order), and `main_source_masks` is how many of them have a best match in that lineage: the main source itself or a genome with ≥ 90% genome-wide AAI to it. Many masks from one lineage point to contamination of the assembly by one organism; masks from several lineages point to separate transfers or a less clear cause. Genomes without high-identity masks have NA and 0; both are NA with `--no-upper-mask`.

`--taxonomy` takes a tab-separated file of genome ID and a taxonomy string with ranks separated by `;`, such as a GTDB taxonomy file (`RS_`/`GB_` prefixes on IDs are removed). The shared lineage ends at the first rank without a name (e.g. `g__`). Taxa are NA when any matching genome has no taxonomy.

### 2.9 Default parameters

| Parameter | Option | Default |
|---|---|---|
| Low-identity test: neighbours (*k*) | `--neighbours` | 3 |
| Low-identity test: running median window (*N*) | `--window-size` | 120 |
| Low-identity test: outlier threshold (*T*, identity points) | `--outlier-threshold` | 25 |
| High-identity test: neighbours (*k*′) | `--upper-neighbours` | 3 |
| High-identity test: running median window (*N*′) | `--upper-window-size` | 120 |
| High-identity test: outlier threshold (*T*′, AAI points) | `--upper-outlier-threshold` | 25 |
| Minimum aligned fraction | `--min-aligned-frac` | 0.5 |
| Minimum aligned columns at the locus (*m*) | `--min-aligned-columns` | 50 |
| Trend method (both tests) | `--method` | `monotone-median` |
| High-identity test | `--no-upper-mask` disables it | on |
| Suspect genome: minimum % of its loci masked (and at least 3) | `--suspect-masked-perc` | 10 |
| Taxonomy file for reporting the source of masked copies | `--taxonomy` | none |

## 3. Evaluation design

### 3.1 Data

- **Example data:** `/srv/projects/donovan/gtrim/data/gtdb_r232_fun77_markers`. 77 fungal marker genes (amino acid MSAs) across 5,813 genomes, with about 5,500 sequences per gene and 36,810 concatenated columns, giving 422,337 genome loci.
- **Planted data:** a copy of the example data with **300 known divergent loci**, made by `analysis/mask_divergent_loci/plant.py` (seed 7):
  - Three bands of 100 plants each, by the donor's pooled AAI to the focal genome: **40–55%, 55–70% and 70–85%**.
  - For each plant, a focal genome and one of its loci are chosen at random, and that locus's sequence is replaced with the sequence of a randomly chosen donor genome in the required AAI band.
  - Each genome is used at most once, as either a focal genome or a donor, so plants don't interact.

**How plants are chosen**
1. **Focal genome:** genomes are taken in a random order.
2. **Locus:** one locus the focal genome has, chosen at random.
3. **Donor:** up to 3,000 randomly ordered genomes that have that locus and haven't been used are tried. The first whose pooled AAI to the focal genome (over all 77 loci) falls in the required band is used. If none does, the next focal genome is tried.
4. **Replacement:** the focal genome's aligned row for that locus is replaced by the donor's row.

The 300 plants span 76 loci, with at most 10 in one locus.

**Divergence of the plants** (medians, with interquartile ranges in brackets). Neighbour AAI and realised deficit come from the k = 3, N = 120 run. The **realised deficit** is the mean of (expected − observed locus identity) over the plant's three comparisons.

| Band | Donor AAI | Focal genome's neighbour AAI | Gap | Realised deficit |
|---|---|---|---|---|
| 40–55% | 51.0 (49.7–52.5), min 45.2 | 92.6 (88.1–96.5) | 41.0 (36.9–45.6) | 39.0 (32.1–47.1) |
| 55–70% | 62.8 (59.1–66.2) | 92.1 (86.0–95.8) | 28.7 (23.0–32.7) | 25.2 (17.9–35.0) |
| 70–85% | 74.0 (71.4–76.0) | 91.8 (87.6–95.5) | 16.8 (11.9–21.2) | 14.8 (7.0–21.9) |
| All | 62.8 (52.5–71.4) | 92.3 (87.4–96.0) | 27.7 (18.5–37.8) | 26.8 (16.7–37.8) |

No donor below 45% AAI was found, so the 40–45% part of the lowest band is empty.

**Why this design.** The bands are a judgement call, not derived from data. Most genomes' closest neighbours have about 90% AAI, and the threshold was T = 25, so the three bands were chosen to give:
- **40–55%:** plants that should be clearly detectable (gap of about 35–50);
- **55–70%:** plants that straddle the threshold (gap of about 20–35);
- **70–85%:** plants mostly below it (gap of about 5–20), to show where detection stops.

The bands are defined by donor AAI because that can be computed when planting. The gap depends on the focal genome's neighbours, which are only known after running the method, so the gap is spread unevenly across bands. Each band has 100 plants so that each band's recall is measured with similar precision, not to reflect how often each kind of divergence occurs.

**Consequences**
- **Overall recall reflects this mix:** it is a property of this deliberately threshold-centred mix (median realised deficit 27 points), not the method's sensitivity on real data. Recall *by gap* (section 4.5) is the transferable result.
- **Real masked loci are more divergent:** on the original data, masked genome loci have a median deficit of 55 points, larger than nearly all plants.
- **What the plants model:** replacing a sequence with another genome's copy models contamination or cross-lineage horizontal transfer. It does not model paralogs, misalignment or partial sequences.

### 3.2 Metrics

| Metric | Definition |
|---|---|
| Masked | Number of masked genome loci on the original data |
| Genes | Number of loci with at least one masked genome |
| Jaccard (J) | Overlap of two masked sets: \|A ∩ B\| / \|A ∪ B\| |
| Consensus | Genome loci masked in at least half of the configurations in a sweep. Consensus recall is the fraction of the consensus a configuration masks; consensus precision is the fraction of its masks that are in the consensus |
| Planted recall | Fraction of the 300 planted loci masked, overall and by AAI band |
| Gap | For a planted locus: mean leave-one-out AAI of the focal genome to its neighbours minus the donor's AAI to the focal genome. This is roughly the expected identity deficit of the planted sequence |
| Excess (knock-on) masks | Non-planted masks in the planted data minus masks in the original data. Estimates the masks that the 300 plants cause in *other* genomes |
| Mirror noise estimate | Comparisons more than *T* points **above** the trend (which divergence can't produce) are used to estimate how many of the comparisons more than *T* points **below** it are noise. Only comparisons where both directions are possible are counted (expected − *T* ≥ 0 and expected + *T* ≤ 100) |

All runs use the `gtrim` binary itself. The metrics are computed from its output files.

### 3.3 Bucketed evaluation

This evaluation measures how often a foreign gene (contamination or lateral transfer) is detected, as a function of how divergent it is and how closely related the genomes involved are. It uses the default settings. Plants are assigned to buckets using values measured *before* planting, so each bucket can be filled to a target size.

**Experiment 1 (low-identity test)**
- **Axes:** *A* = the focal genome's leave-one-out AAI to its 3 nearest relatives at the locus (from a baseline run), and *G* = identity between the planted copy and the native copy it replaces.
- **Bins:** *A*: <70, 70–80, 80–85, 85–90, 90–95, 95–98, 98–100. *G*: <40, 40–50, 50–60, 60–70, 70–80, 80–90, 90–95.
- **Buckets used:** those where the upper bound of *G* is at most the upper bound of *A* minus 5 (36 buckets).
- **Foreign donors:** AAI(focal, donor) ≤ *A* − 5, and *G* ≤ *A* − 5.

**Experiment 2 (high-identity test)**
- **Axes:** *B* = AAI between the focal genome and the best remaining match of the planted copy, and *I* = identity of the planted copy to that best match.
- **Bins:** *B*: <50, 50–60, 60–70, 70–80, 80–90. *I*: <80, 80–90, 90–95, 95–98, 98–100 (25 buckets).
- **Setting *I*:** the donor's sequence, and any sequence at the locus more identical to the planted copy than the upper bound of the *I* bin, are removed so the best match falls in the bin. A plant may remove up to 100 sequences, and a locus loses at most 2% of its sequences per run, so well-sampled donor lineages can be used in the low-*I* buckets without any locus being hollowed out.
- **Foreign donors:** AAI(focal, donor) ≤ *A* − 10.
- **Bucket assignment:** final buckets come from the best match actually reported by the run.

**Targets:** 500 plants per bucket, and 200 for the *A* < 70% row, which has only 136 genomes. Genomes in that row are reused across runs, but at different loci: unused (locus, genome) pairs are always planted first. No locus contributes more than 15 plants to a bucket. 95% confidence intervals come from resampling focal genomes (a cluster bootstrap); with 500 plants the half-width is about ±4.4 points, and wider in the reused rows, which are also reported pooled.

**Design choices made after test runs on the full data.** Early test runs showed that packing plants densely created artificial knock-on masks of donor lineages, so:
- **Low density:** at most 0.25% of the genomes with a locus are planted at that locus per run (13 per locus for this data), one plant per focal genome per run (a planted locus would otherwise lower the genome-wide AAI used at the genome's other plants), and focal and donor genomes are disjoint.
- **Unrelated donors per locus:** donors planted at the same locus in a run come from different lineages (donor copies < 90% identical and donor genomes < 80% AAI to each other), and each donor copy is used once per run.
- **Donor removed in both experiments:** the donor's own sequence at the locus is removed. The source of real contamination is rarely in the dataset, and a donor that keeps its copy sees that copy in a foreign genome.
- **No planted neighbours:** genomes that are low-identity test neighbours of each other at a locus (in the baseline) are not both planted there. Plants whose comparisons in either test still include another planted genome at the same locus are excluded from the analysis.
- **Controls:** each run has a control dataset with the same removals and no plants. Knock-on masks = masks of non-planted genome loci in the planted run minus those in the control.

A one-run test of Experiment 1 at this density gave about 26 high-identity knock-on masks per 1,000 plants; a one-run test of Experiment 2 gave none.

**Reported per bucket:** recall of the evaluated test, the other test and either test, with confidence intervals and the number of distinct genomes and loci; row totals; recall against the realised deficit; recall recomputed at other thresholds (expected values don't depend on T or T′, so no reruns are needed); per-locus recall; and knock-on masks per 1,000 plants.

Three further views address known limits of the design:
- **Identity to neighbours (Experiment 1):** *G* is identity to the replaced native copy, but detection depends on the planted copy's identity to the focal genome's neighbours. A second grid bins plants by *A* and the mean identity to the actual neighbours, taken from the run output.
- **Locus conservation strata:** at a conserved locus even distant genomes exceed 80% identity, so low-*I* buckets can only be filled from variable loci, and the locus mix differs between buckets. Recall grids are also given within tertiles of locus conservation (median neighbour identity per locus in the baseline).
- **Knock-on by donor lineage depth:** knock-on masks are the method's blame ambiguity in sparsely sampled lineages (a lineage member whose best identity match at the locus becomes the planted genome). They are attributed to the responsible plant and reported by the donor lineage's sampling depth (other genomes with ≥ 90% AAI to the donor), where the cost should vanish for lineages sampled two or more deep.

**Scale:** 17,100 plants in 18 runs (Experiment 1) and 12,500 plants in roughly 15–20 runs (Experiment 2; later runs yield fewer plants as buckets fill). Each run also has a control, so about 70 `gtrim` runs in all: roughly 1–1.5 hours on 96 cores and about 11 GB under `WORK`.

## 4. Results

### 4.1 Default run on the example data

k = 3, N = 120, T = 25, 96 threads:

- **Low-identity test:** 1,604 of 422,337 genome loci masked (0.38%), from 6,543 outlier comparisons. 70 of the 77 genes have at least one masked genome.
- **High-identity test** (section 4.8): 112 masked, from 690 outlier comparisons; 21 of these are also masked by the low-identity test.
- **Total:** 1,695 genome loci masked (0.40%).
- **Resources:** with `--plot`, 44 s and 1.7 GB of RAM. The low-identity test alone took 24 s and 1.0 GB without plots, and 31 s and 1.4 GB with them.

Sections 4.1–4.7 below evaluate the **low-identity test alone**. They were produced before the high-identity test existed, and the sweep scripts pass `--no-upper-mask`.

Masked genomes are concentrated in a few genes:

| Gene | Masked genomes |
|---|---|
| 2934341at4751 | 492 (9.1%) |
| 3323237at4751 | 129 |
| 3453803at4751 | 95 |
| 94760at4751 | 75 |
| 3456310at4751 | 73 |
| 520361at4751 | 65 |
| 3459837at4751 | 65 |

Many masked sequences have locus identity of 5–30% to neighbours with 85–98% AAI. For example, GCA_000142805.1 at 2934341at4751 has 5.4% identity over 838 aligned columns to a neighbour with 87.5% AAI. That is below what unrelated proteins share, which suggests paralogs or misalignment rather than gradual divergence.

### 4.2 Sweep 1: number of neighbours (k) × window size (N), T = 25

**Masked genome loci (original data)**

| k | N=20 | N=40 | N=80 | N=120 | N=200 | N=400 |
|---|---|---|---|---|---|---|
| 1 | 1,924 | 1,986 | 2,058 | 2,159 | 2,369 | 2,845 |
| 3 | 1,486 | 1,531 | 1,567 | **1,604** | 1,635 | 1,750 |
| 5 | 1,420 | 1,481 | 1,510 | 1,530 | 1,557 | 1,630 |
| 7 | 1,399 | 1,460 | 1,495 | 1,513 | 1,527 | 1,597 |

**Window size (N) has a modest effect:**
- **Counts:** at k = 3, the masked count rises by 18% from N = 20 to N = 400.
- **Agreement:** the masked sets are very similar for N from 20 to 200 (J = 0.87–0.95 at k = 3). N = 80, 120 and 200 agree most closely with each other (J = 0.93–0.95).
- **N = 400:** inflates the masked count, especially at k = 1.
- **Genes:** 12–13 genes have 25 or more genomes masked in every k ≥ 3 configuration.

**The number of neighbours (k) matters more:**

| k (N=120) | J vs k=3 | Consensus recall | Consensus precision |
|---|---|---|---|
| 1 | 0.49 | 0.78 | 0.57 |
| 3 | 1.00 | 0.92 | 0.90 |
| 5 | 0.73 | 0.91 | 0.93 |
| 7 | 0.65 | 0.85 | 0.88 |

Across all 24 configurations:
- **Core:** 919 genome loci are masked in every configuration; most of these are clear cases.
- **Singletons:** 534 are masked by only one configuration, mostly k = 1.
- **Consensus:** 1,573 genome loci are masked in at least half.

### 4.3 Planted divergent loci (sweep 1)

**Recall and knock-on masks, N = 120, T = 25**

| k | Donor AAI 40–55% | 55–70% | 70–85% | All | Excess (knock-on) masks |
|---|---|---|---|---|---|
| 1 | 0.94 | 0.59 | 0.23 | 0.59 | **+166** |
| 3 | 0.94 | 0.51 | 0.20 | 0.55 | −9 |
| 5 | 0.92 | 0.46 | 0.18 | 0.52 | −2 |
| 7 | 0.88 | 0.47 | 0.18 | 0.51 | +9 |

- **N barely affects recall:** at fixed k, overall recall varies by at most 0.02 across N, and recall within a band by at most 0.04.
- **Detection depends on the gap, not on k or N.** At k = 3, N = 120:

| Donor AAI band | Detected: mean neighbour AAI | Detected: mean gap | Missed: mean neighbour AAI | Missed: mean gap |
|---|---|---|---|---|
| 40–55% | 91.8 | 40.8 | 81.4 | 29.4 |
| 55–70% | 92.8 | 31.0 | 87.2 | 24.1 |
| 70–85% | 95.4 | 22.2 | 88.5 | 13.7 |

Plants are detected when their gap clearly exceeds *T* = 25 and missed when it is near or below *T*. All missed plants received their full set of comparisons. So for k ≥ 3, the outlier threshold is what controls sensitivity.

### 4.4 Why k = 1 performs poorly

With k = 1, a genome locus is masked if its single comparison with its nearest neighbour is an outlier. A single comparison can't show which genome in the pair carries the divergent copy. Its nearest neighbour is also masked whenever that neighbour's own nearest-neighbour comparison is an outlier, which is usually the comparison back to the divergent genome.

Breakdown of the 2,159 k = 1 masks (N = 120):

| Partner of a masked genome | Count |
|---|---|
| Also masked, mutual nearest neighbours | 996 |
| Also masked, through a different nearest neighbour | 391 |
| Not masked | 772 |

**Planted data:** k = 1 detected 176 plants and also masked 173 other genomes only because their nearest neighbour was a detected plant (0–4 per plant). For an isolated divergent locus, k = 1 therefore roughly doubles the number of masks.

**Original data:** k = 1 gives 1.35× as many masks as k = 3, not 2×, because the two sets differ in both directions:
- **k = 1 adds 920 masks that k = 3 doesn't make:**
  - 582 are knock-on masks of genomes whose nearest neighbour is the one masked at k = 3;
  - 338 are single-comparison noise, 255 of them within 5 points of the threshold.
- **k = 1 misses 365 masks that k = 3 makes,** all because the nearest-neighbour comparison is not an outlier. In 279 of these the nearest neighbour is also masked at k = 3. These are close relatives that share the same divergent copy, which only more distant neighbours reveal.
- **Not every divergent locus causes a knock-on mask:** 41% of k = 3-masked genome loci are nobody's nearest neighbour.

### 4.5 Sweep 2: k = 2 and 4, and outlier threshold

Design:
- **k × N grid:** k ∈ {2, 4} × N ∈ {20, 40, 80, 120, 200, 400} at T = 25, which completes the grid for k ∈ {1, 2, 3, 4, 5, 7}.
- **Threshold sweep:** T ∈ {10, 15, 20, 25, 30, 35, 40} for k ∈ {1, 2, 3, 4, 5, 7} at N = 120.
- **Datasets:** both the original and the planted data, 96 runs in total.

#### Even k: k = 2 and k = 4

Masked genome loci on the original data, T = 25:

| k | Votes needed | N=20 | N=40 | N=80 | N=120 | N=200 | N=400 |
|---|---|---|---|---|---|---|---|
| 1 | 1 of 1 | 1,924 | 1,986 | 2,058 | 2,159 | 2,369 | 2,845 |
| 2 | 2 of 2 | 1,058 | 1,105 | 1,119 | 1,154 | 1,210 | 1,396 |
| 3 | 2 of 3 | 1,486 | 1,531 | 1,567 | 1,604 | 1,635 | 1,750 |
| 4 | 3 of 4 | 1,169 | 1,218 | 1,250 | 1,280 | 1,301 | 1,379 |
| 5 | 3 of 5 | 1,420 | 1,481 | 1,510 | 1,530 | 1,557 | 1,630 |
| 7 | 4 of 7 | 1,399 | 1,460 | 1,495 | 1,513 | 1,527 | 1,597 |

The counts zigzag between odd and even k because of the majority rule. Even k needs a larger fraction of outlier comparisons (2 of 2 = 100%, 3 of 4 = 75%) than odd k (2 of 3 = 67%, 3 of 5 = 60%), so even k is stricter.

Comparison at N = 120, T = 25. The consensus here is over all 36 k × N configurations: the 1,399 genome loci masked in at least half of them.

| k | Masked | Genes | J vs k=3 | Consensus recall | Consensus precision | Planted recall | Knock-on masks |
|---|---|---|---|---|---|---|---|
| 1 | 2,159 | 73 | 0.49 | 0.84 | 0.54 | 0.59 | +166 |
| 2 | 1,154 | 69 | 0.70 | 0.79 | 0.95 | 0.54 | −6 |
| 3 | 1,604 | 70 | 1.00 | 0.99 | 0.86 | 0.55 | −9 |
| 4 | 1,280 | 67 | 0.79 | 0.89 | 0.98 | 0.52 | −5 |
| 5 | 1,530 | 67 | 0.73 | 0.93 | 0.85 | 0.52 | −2 |
| 7 | 1,513 | 66 | 0.65 | 0.87 | 0.80 | 0.51 | +9 |

- **No knock-on masks for k ≥ 2:** requiring 2 of 2 is enough to stop the knock-on effect seen with k = 1.
- **k = 2 and k = 4 are the most conservative settings:** 95–98% of their masks are in the consensus. They miss more of it (79–89%) and have slightly lower planted recall than k = 3.
- **k = 3 is the most permissive setting without knock-on masks:** it recovers 99% of the consensus and has the highest planted recall among k ≥ 2.
- **The consensus definition matters:** with 36 configurations, the stricter even-k settings shrink the consensus (1,399, compared with 1,573 in sweep 1), which lowers the consensus precision of k = 3 to 0.86.

#### Outlier threshold (N = 120)

**Masked genome loci on the original data**

| T | k=1 | k=2 | k=3 | k=4 | k=5 | k=7 |
|---|---|---|---|---|---|---|
| 10 | 12,166 | 7,414 | 11,635 | 9,438 | 12,107 | 12,603 |
| 15 | 5,616 | 3,199 | 4,876 | 3,980 | 5,052 | 5,228 |
| 20 | 3,207 | 1,724 | 2,529 | 2,049 | 2,527 | 2,557 |
| 25 | 2,159 | 1,154 | **1,604** | 1,280 | 1,530 | 1,513 |
| 30 | 1,654 | 863 | 1,151 | 932 | 1,095 | 1,051 |
| 35 | 1,397 | 744 | 960 | 798 | 901 | 869 |
| 40 | 1,282 | 688 | 878 | 736 | 823 | 799 |

**Recall of planted loci (all 300)**

| T | k=1 | k=2 | k=3 | k=4 | k=5 | k=7 |
|---|---|---|---|---|---|---|
| 10 | 0.89 | 0.86 | 0.87 | 0.84 | 0.85 | 0.83 |
| 15 | 0.80 | 0.78 | 0.78 | 0.75 | 0.75 | 0.72 |
| 20 | 0.72 | 0.65 | 0.66 | 0.62 | 0.63 | 0.61 |
| 25 | 0.59 | 0.54 | 0.55 | 0.52 | 0.52 | 0.51 |
| 30 | 0.48 | 0.42 | 0.42 | 0.38 | 0.39 | 0.37 |
| 35 | 0.37 | 0.32 | 0.33 | 0.29 | 0.29 | 0.25 |
| 40 | 0.24 | 0.20 | 0.20 | 0.19 | 0.19 | 0.18 |

**Recall by gap (k = 3).** The number of planted loci in each gap bin is shown in brackets.

| T | Gap <15 (44) | 15–20 (40) | 20–25 (49) | 25–30 (31) | 30–35 (38) | 35–40 (38) | 40+ (60) |
|---|---|---|---|---|---|---|---|
| 10 | 0.41 | 0.72 | 0.98 | 0.97 | 0.97 | 1.00 | 1.00 |
| 15 | 0.23 | 0.60 | 0.78 | 0.90 | 0.95 | 1.00 | 1.00 |
| 20 | 0.09 | 0.28 | 0.57 | 0.81 | 0.92 | 0.89 | 1.00 |
| 25 | 0.05 | 0.07 | 0.45 | 0.61 | 0.68 | 0.87 | 1.00 |
| 30 | 0.02 | 0.03 | 0.27 | 0.26 | 0.53 | 0.79 | 0.90 |
| 35 | 0.00 | 0.00 | 0.16 | 0.13 | 0.39 | 0.66 | 0.78 |
| 40 | 0.00 | 0.00 | 0.06 | 0.10 | 0.21 | 0.34 | 0.55 |

k = 2 gives nearly identical values. A planted locus is detected reliably (recall ≥ 0.8) once its gap exceeds about T + 10. About half are detected when the gap is between T and T + 5 (recall 0.53–0.61 in that bin for T = 15–30). The gap is calculated from genome-wide AAI, while the trend line is fitted per locus, so the actual deficit of a planted locus from its expected identity differs somewhat from its gap.

**Mirror noise estimate** (k = 3, original data). The overall robust SD of deviations from the trend line is 2.1 points.

| T | Comparisons > T below trend | Comparisons > T above trend | Genome loci masked (eligible comparisons) | Mirror-masked genome loci | Estimated noise fraction |
|---|---|---|---|---|---|
| 10 | 26,745 | 7,565 | 7,625 | 1,984 | 0.26 |
| 15 | 7,082 | 1,231 | 1,914 | 318 | 0.17 |
| 20 | 1,938 | 262 | 514 | 61 | 0.12 |
| 25 | 584 | 81 | 149 | 23 | 0.15 |
| 30 | 149 | 49 | 32 | 15 | (0.47) |

**How it works.** If random scatter around the trend line is symmetric, a comparison is as likely to be T points above the trend by chance as T points below. Divergence (contamination, horizontal transfer, paralogs) only lowers identity. So the number of comparisons more than T points *above* the trend estimates how many of those more than T points *below* it are chance, and the same mirror logic is applied to the majority vote.

The counts are restricted to comparisons whose expected identity lies between T and 100 − T, where both directions are possible. A comparison expected at 95% can't be 25 points above it.

**What limits it** (k = 3, T = 25)
- **Coverage is small.** Only 5.4% of comparisons qualify (expected identity 25–75%; median expected identity overall is 94.5%). That covers 584 of 6,542 outlier comparisons and 149 of 1,604 masked genome loci. **The other 91% of masks aren't covered.**
- **It measures the noisiest region.** Scatter shrinks as identity rises:

| Expected identity | 25–75% | 75–85% | 85–90% | 90–95% | 95–100% |
|---|---|---|---|---|---|
| Robust SD of deviation (points) | 6.8 | 5.4 | 4.2 | 2.9 | 1.1 |

  At T = 25 an outlier is about 3.7 robust SDs out in the covered region, but more than 6 SDs out where expected identity exceeds 85%. The uncovered masks are also mostly far from the cutoff: the median deficit of outlier comparisons in masked genome loci is 55 points.
- **The "above" side is mostly not noise.** At T = 25 there are 81 above-trend outlier comparisons, against 6,542 below. They are concentrated, full-length and near-identical, which points to contamination, chimeric assemblies or horizontal transfer rather than chance:
  - **Concentrated:** they fall in only 23 genomes, and one genome (GCA_023085455.1) accounts for 40 of them. That genome also has 24 below-trend outliers and a median leave-one-out AAI to its neighbours of only 64%. It is 98–99% identical to GCA_004016305.1 at some loci and to GCA_022204945.1 at others, which looks like a chimeric or contaminated assembly.
  - **Full-length:** the median aligned fraction of above-trend outliers is 1.00, so they aren't artefacts of partial overlap.

  These cases make the mirror estimate overstate noise. See section 4.6.
- **The scatter isn't perfectly symmetric.** Even close to the trend, deviations below outnumber those above (27,261 vs 24,185 beyond ±2 points in the covered region). Noise that can only lower identity, such as local misalignment or short overlaps, isn't measured at all. Much of it is arguably what the method should catch.

**Interpretation.** Among masks in the low-identity region, at most about 12–17% at T = 15–25 (and about 26% at T = 10) are consistent with symmetric chance scatter. Across all masks the fraction is very likely much lower, because most masks lie where scatter is smaller and deficits are larger. The estimate is not reliable above T = 25.

**Other observations**
- **Masked sets are nested across T:** at k = 3, every genome locus masked at a higher T is also masked at a lower T, so changing T only adds or removes the least extreme cases.
- **Knock-on masks:** they stay large for k = 1 at every threshold (+77 to +233). For k ≥ 2 they stay within about ±40 of zero.

#### Trade-off at k = 3

| T | Masked (% of genome loci) | Planted recall | Recall for gap 20–25 | Recall for gap 25–30 | Estimated noise fraction |
|---|---|---|---|---|---|
| 15 | 4,876 (1.15%) | 0.78 | 0.78 | 0.90 | ≤ 0.17 |
| 20 | 2,529 (0.60%) | 0.66 | 0.57 | 0.81 | ≤ 0.12 |
| 25 | 1,604 (0.38%) | 0.55 | 0.45 | 0.61 | ≤ 0.15 |
| 30 | 1,151 (0.27%) | 0.42 | 0.27 | 0.26 | not reliable |

### 4.6 Loci more similar than expected

*Superseded by the high-identity test (sections 2.5 and 4.8). This section records the analysis that motivated it.*

The low-identity test is one-sided. A locus that is far *more* similar than genome-wide relatedness predicts is just as inconsistent with vertical descent: for example, 99% identity between genomes with 65% AAI. It is most likely contamination, a chimeric assembly or horizontal transfer, and it would pull the genome towards the donor lineage in a tree.

Data at k = 3, N = 120, T = 25:

| Rule | Genome loci masked |
|---|---|
| Current: majority of comparisons below by > T | 1,603 |
| Majority above by > T | 23 |
| Majority outlying in either direction | 1,626 |
| At least one comparison above by > T | 45 (22 with 1 of 3, 10 with 2 of 3, 13 with 3 of 3) |

Design issues:
- **Limited reach.** An above-trend outlier is only possible where expected identity is at most 100 − T, which applies to 5.4% of comparisons at T = 25. The test therefore only helps for genomes without close relatives, which is where the below test is also weakest.
- **The majority vote is a poor fit.** A copy acquired from the lineage of *one* neighbour is above-trend only against that neighbour; against the others it usually looks normal. So majority voting misses it: 22 of the 45 genome loci with an above outlier have just 1 of 3.
- **Which genome to blame.** A single above-trend pair (f, n) could mean f carries n's lineage's copy, or n carries f's lineage's copy. In the second case n's locus would normally be a below outlier against n's own close relatives, so this can be resolved by checking whether n's locus is masked.
- **Chimeric genomes.** Genomes like GCA_023085455.1, with many above *and* below outliers, are better excluded entirely than masked locus by locus.

### 4.7 Data-driven thresholds (exploratory)

**Question:** can T be set from the data, and separately for each locus, without planted loci?

**Approach tested:** the threshold becomes T = *c* × *s*, where *s* is the local robust SD of the residuals (deficits). It is estimated for each locus over the same N = 120 nearest points (in genome-wide AAI) used for the running median: *s* = 1.4826 × median |residual|, with a floor of 1.0 point.

The analysis was done in Python on the k = 3, N = 120 comparison tables (`analysis/mask_divergent_loci/analyse_scale_threshold.py`). It is not implemented in `gtrim`.

**The scale differs a lot across the data:**

| Genome-wide AAI | <75 | 75–85 | 85–90 | 90–95 | 95–100 |
|---|---|---|---|---|---|
| Median local robust SD (points) | 5.9 | 4.9 | 3.8 | 2.5 | 1.1 |
| Implied T at *c* = 6 (median) | 35.6 | 29.5 | 22.5 | 15.0 | 6.8 |

Loci differ too: the median implied T at *c* = 6 ranges from 6 to 44 points across the 77 loci (median 15). A fixed T = 25 is about 4 SDs in the noisiest region but more than 20 SDs among close relatives.

**Comparison with a fixed T** (k = 3 majority vote, low-identity outliers only):

| Rule | Masked | J vs T=25 | Planted recall | Gap 15–20 | Gap 20–25 | Gap 25–30 | Gap 30–40 | Knock-on |
|---|---|---|---|---|---|---|---|---|
| T = 25 (fixed) | 1,603 | 1.00 | 0.55 | 0.07 | 0.45 | 0.61 | 0.78 | −8 |
| T = 20 (fixed, section 4.5) | 2,529 | 0.63 | 0.66 | 0.28 | 0.57 | 0.81 | — | +1 |
| *c* = 4 | 5,085 | 0.26 | 0.80 | 0.60 | 0.84 | 0.90 | 0.99 | −39 |
| *c* = 5 | 3,190 | 0.36 | 0.73 | 0.50 | 0.69 | 0.77 | 0.95 | −27 |
| *c* = 6 | 2,238 | 0.43 | 0.69 | 0.38 | 0.59 | 0.77 | 0.93 | +9 |
| *c* = 6, deficit > 10 | 1,921 | 0.48 | 0.68 | 0.33 | 0.59 | 0.77 | 0.92 | +4 |
| *c* = 8 | 1,452 | 0.50 | 0.58 | 0.15 | 0.47 | 0.58 | 0.82 | +17 |
| *c* = 10 | 1,122 | 0.50 | 0.49 | 0.07 | 0.39 | 0.42 | 0.67 | +2 |

The section 4.5 gap bins differ slightly, so the T = 20 values for gap 30–40 aren't directly comparable.

- **Better recall for the same number of masks:** *c* = 8 masks fewer genome loci than T = 25 (1,452 vs 1,603) with higher recall (0.58 vs 0.55). *c* = 6 masks fewer than T = 20 (2,238 vs 2,529) with higher recall (0.69 vs 0.66). *c* = 10 vs T = 30 is 1,122 vs 1,151 masks with recall 0.49 vs 0.42.
- **The masks move:** the masked sets differ substantially from T = 25 (J = 0.43–0.50). Scaled thresholds mask more where genomes are closely related (e.g., 92% identity where 99% is expected) and less among distant relatives, where scatter is largest.
- **Caveat:** the plants are concentrated at neighbour AAI of about 92%, where scaled thresholds are low, which favours this approach.

**What the data can and can't determine**
- **The scale *s* can be estimated reliably:** per locus and as it varies along the AAI axis.
- **The multiplier *c* can't be determined from the data alone.** The deficit distribution has a distinct second mode at 55–70 points (near-unrelated sequences: paralogs or contaminants) with a valley at 40–45 points. Below that, there's a smooth continuum from noise into moderate divergence, with no break between them:

| Deficit (points) | 5–10 | 10–15 | 15–20 | 20–25 | 25–30 | 30–35 | 35–40 | 40–45 | 45–50 | 50–55 | 55–60 | 60–65 | 65–70 | 70–75 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Comparisons | 80,154 | 22,950 | 8,203 | 3,306 | 1,590 | 729 | 319 | 131 | 167 | 267 | 497 | 726 | 1,400 | 225 |

- **Standard statistical cutoffs don't apply.** The tails are far heavier than a normal distribution: 21,208 comparisons have z > 4, where a normal distribution would give 80. An empirical-null fit (fitting the noise distribution to the central part of the data) is hampered by the 100% ceiling, which truncates the upper side: there are only 570 comparisons with z < −3.

### 4.8 High-identity test

**Example data** (k = 3, N = 120, T = 25):
- **Masks:** 112 genome loci, 91 of them not masked by the low-identity test. The outlier comparisons fall 20–60 AAI points below expectation. Of the 690 outlier comparisons, 75% have locus identity ≥ 90% and 82% have AAI of 50–70%; 61% meet both conditions.
- **Concentrated in a few genomes:** GCA_023085455.1, the suspected chimeric assembly from section 4.6, has 21 loci masked by this test and 8 by the low-identity test, so 29 of its 44 loci are masked. GCA_029303655.1 has 10 of 32 loci masked by this test.
- **Concentrated in a few loci:** locus 520361at4751 has 33 masks. For example, GCA_000442015.1 is 77% identical there to three genomes it shares 53% AAI with, against an expected AAI of 82% for that identity.

**Planted loci with donors removed.** The original planted data contains each donor genome, and a planted copy is 100% identical to its donor, so this test would find every plant trivially. For a fair evaluation, the same 300 plants (seed 7) were made with **all donor genomes removed from the data**, leaving the donor's lineage represented only by its relatives. A control dataset (original data minus the donors) gives the knock-on masks.

| Donor AAI band | Low-identity test | High-identity test | Either test | Found only by high-identity test |
|---|---|---|---|---|
| 40–55% | 0.92 | 0.91 | 0.98 | 0.06 |
| 55–70% | 0.49 | 0.76 | 0.85 | 0.36 |
| 70–85% | 0.21 | 0.12 | 0.28 | 0.07 |
| All | 0.54 | 0.60 | **0.70** | 0.16 |

| Gap (points) | <15 | 15–20 | 20–25 | 25–30 | 30–40 | 40+ |
|---|---|---|---|---|---|---|
| n | 46 | 38 | 47 | 32 | 75 | 57 |
| Low-identity test | 0.04 | 0.11 | 0.47 | 0.59 | 0.77 | 1.00 |
| High-identity test | 0.15 | 0.18 | 0.45 | 0.84 | 0.87 | 0.91 |
| Either test | 0.20 | 0.29 | 0.70 | 0.91 | 0.96 | 1.00 |

**Detection depends on how well the donor lineage is represented.** The table below groups plants by the best locus identity of the planted copy to any genome remaining in the data:

| Best remaining identity | 0–80% | 80–90% | 90–95% | 95–98% | 98–100% |
|---|---|---|---|---|---|
| n | 26 | 31 | 54 | 66 | 118 |
| High-identity test recall | 0.15 | 0.42 | 0.65 | 0.64 | 0.72 |

**Knock-on masks:** essentially none. The control has 107 high-identity masks; the planted data has 110 non-plant high-identity masks (+3), and the low-identity test adds 0. So the high-identity test attributes contamination to the right genome.

**Summary:** adding the high-identity test raises recall from 0.54 to 0.70 without knock-on masks. The gain is largest for moderately divergent copies (donor AAI 55–70%: 0.49 → 0.85; gap 25–30: 0.59 → 0.91) whose source lineage is still represented in the data.

### 4.9 High-identity test parameters (k′, N′, T′)

**Design:**
- **Data:** the planted data with donor genomes removed, and its control (section 4.8). The control stands in for real data when counting masks and knock-on masks.
- **Fixed:** the low-identity test stays at k = 3, N = 120, T = 25. On the control it gives 1,521 masks; its planted recall is 0.54.
- **Grid:** k′ ∈ {1, 2, 3, 4, 5, 7} × N′ ∈ {20, 40, 80, 120, 200, 400} at T′ = 25, plus T′ ∈ {10, 15, 20, 30, 35, 40} for every k′ at N′ = 120.

That is 144 runs, which took 37 minutes running 4 at a time.

**Neighbours and window (T′ = 25)**

| k′ | Votes needed | Masks on control (N′=20 / 120 / 400) | Planted recall (N′=120) | Knock-on (N′=120) |
|---|---|---|---|---|
| 1 | 1 of 1 | 273 / 363 / 998 | 0.66 | +170 |
| 2 | 2 of 2 | 85 / 94 / 188 | 0.60 | +4 |
| 3 | 2 of 3 | 116 / **107** / 162 | **0.60** | +3 |
| 4 | 3 of 4 | 94 / 92 / 105 | 0.54 | +1 |
| 5 | 3 of 5 | 104 / 107 / 104 | 0.54 | −3 |
| 7 | 4 of 7 | 93 / 97 / 100 | 0.50 | −1 |

- **k′ = 1 behaves like k = 1 in the low-identity test.** It adds about 170 knock-on masks (genomes blamed for a neighbour's foreign copy), and its masks inflate with N′.
- **k′ ≥ 2 has no knock-on masks.** k′ = 2 and 3 have the highest recall (0.60). k′ ≥ 4 loses recall without reducing masks.
- **N′ barely matters:** for k′ ≥ 2, recall varies by at most 0.02 across N′ (0.04 for k′ = 1). Masks are stable between N′ = 20 and 200; N′ = 400 inflates them for k′ ≤ 3.

**Threshold (N′ = 120)**

| T′ | Masks on control (k′=3) | Not masked by low-identity test | High-identity recall (k′=3) | Either test (k′=3) | Knock-on (k′=3) |
|---|---|---|---|---|---|
| 10 | 5,185 | 5,064 | 0.88 | 0.90 | −33 |
| 15 | 1,027 | 955 | 0.83 | 0.87 | −4 |
| 20 | 269 | 227 | 0.73 | 0.80 | 0 |
| **25** | **107** | 85 | **0.60** | **0.70** | +3 |
| 30 | 57 | 45 | 0.43 | 0.62 | −2 |
| 35 | 29 | 19 | 0.27 | 0.57 | 0 |
| 40 | 16 | 8 | 0.17 | 0.55 | 0 |

Values for other k′ follow the same pattern; k′ = 2 is nearly identical to k′ = 3, and k′ ≥ 4 is lower throughout. For k′ ≥ 2, knock-on masks stay within ±10 of zero for T′ ≥ 15. Masked sets are nested across T′: every locus masked at a higher T′ is also masked at a lower one.

**High-identity recall by how closely the donor lineage is still represented** (best remaining locus identity of the planted copy; k′ = 3, N′ = 120):

| T′ | 0–80% (26) | 80–90% (31) | 90–95% (54) | 95–98% (66) | 98–100% (118) |
|---|---|---|---|---|---|
| 15 | 0.38 | 0.68 | 0.89 | 0.89 | 0.94 |
| 20 | 0.31 | 0.55 | 0.74 | 0.80 | 0.86 |
| 25 | 0.15 | 0.42 | 0.65 | 0.64 | 0.72 |
| 30 | 0.08 | 0.29 | 0.46 | 0.45 | 0.53 |

**What extra masks at lower T′ look like.** The expected AAI doesn't depend on T′, so the kept comparison table (k′ = 3, N′ = 120) shows which genome loci each threshold masks on the control:

| Masks | n | Median locus identity | Median AAI | Median expected AAI | Median deficit | Comparisons with identity ≥ 95% |
|---|---|---|---|---|---|---|
| T′ = 25 | 107 | 96.3 | 59.4 | 93.3 | 30.6 | 54% |
| Added at T′ = 20 | 162 | 81.1 | 59.9 | 81.3 | 21.7 | 13% |
| Added at T′ = 15 | 758 | 79.3 | 62.7 | 79.3 | 16.4 | 11% |

The T′ = 25 masks look like contamination: a near-identical copy in a genome with only about 60% AAI to the source. Masks added at lower T′ are a different kind of case. Their locus identity is moderate (about 80%), at AAI of about 60% where about 80% would be expected. That is consistent with transfer from an intermediate lineage, but also with the gene evolving slowly in that lineage, which doesn't conflict with vertical descent. The planted loci only model the first case, so recall at low T′ is likely to overstate the benefit.

**Genes:** the masks are dominated by 520361at4751 at every threshold (29 masks at T′ = 25, 61 at T′ = 20). Lower thresholds add 3472996at4751, 97101at4751 and 94528at4751.

### 4.10 Bucketed evaluation, Experiment 1: low-identity test

Design in section 3.3. Default settings (k = 3, N = 120, T = 25, both tests on). 17,100 plants in 18 runs, each with a control. 547 plants were excluded because a neighbour in either test was also planted, leaving 16,553. The 36 gtrim runs took 11.5 minutes on 96 cores. Full report: `$WORK/buckets_low/report.txt`; per-bucket recall with 95% CIs (cluster bootstrap over focal genomes): `$WORK/buckets_low/results.tsv`.

**Plants per bucket:** 434–500 (188–193 in the A < 70 row), from 128–455 distinct focal genomes and 63–77 loci per bucket. CI half-widths are about ±0.04, and ±0.05–0.08 in the A < 70 row.

**How real loci are distributed over the same grid** (`real_g_distribution.py`). For each genome locus in the unmodified data, A is computed as for the plants and G is its mean identity to its 3 neighbours, the quantity the test sees. Counts are real genome loci, with those masked by the low-identity test in brackets:

| A | Loci | Median G | 1st pct of G | G <40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 | 90–95 | 95–100 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 98–100 | 49,167 | 99.0 | 89.2 | 50 (50) | 6 (6) | 39 (39) | 25 (16) | 100 (8) | 351 (1) | 2,089 | 46,507 |
| 95–98 | 102,532 | 97.2 | 81.5 | 136 (136) | 11 (11) | 80 (79) | 82 (22) | 567 (18) | 4,407 (2) | 18,297 | 78,952 |
| 90–95 | 121,018 | 94.1 | 72.2 | 200 (200) | 25 (25) | 85 (56) | 565 (71) | 2,978 (28) | 22,631 | 46,108 | 48,426 |
| 85–90 | 69,541 | 89.9 | 63.4 | 103 (102) | 50 (31) | 272 (74) | 1,300 (73) | 6,810 (16) | 26,709 | 25,528 | 8,769 |
| 80–85 | 42,853 | 85.6 | 56.7 | 97 (90) | 116 (40) | 453 (43) | 2,404 (41) | 8,192 (12) | 20,932 | 9,136 | 1,523 |
| 70–80 | 28,806 | 79.9 | 47.3 | 105 (61) | 321 (71) | 1,284 (60) | 4,036 (28) | 8,770 (1) | 11,648 | 2,411 | 231 |
| <70 | 8,420 | 66.7 | 32.4 | 315 (50) | 776 (24) | 1,640 (15) | 2,313 (4) | 2,217 | 1,025 | 120 | 14 |

Cells without brackets have no masks. Example, A 98–100 and G 70–80: plant recall there is 0.49 (below), but only 100 of the row's 49,167 real loci fall there (0.2%), and 8 are masked. Those 100 are concentrated in 2934341at4751 (53), 3459837at4751 (18) and 3323237at4751 (13), and their G sits near the top of the bin: median realised deficit 22.9, against 24.9 for the plants in that bucket (7 vs 240 of 500 above 25). Plant recall describes what happens to a foreign copy that lands in a cell, not how many real loci are there; real masks are concentrated in the far tail, especially G < 40.

**Recall of the low-identity test** (rows A, columns G; about 500 plants per cell, about 190 in the A < 70 row):

| A \ G | <40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 | 90–95 | Row |
|---|---|---|---|---|---|---|---|---|
| 98–100 | 1.00 | 1.00 | 1.00 | 0.99 | 0.49 | 0.00 | 0.00 | 0.64 |
| 95–98 | 1.00 | 1.00 | 1.00 | 0.98 | 0.27 | 0.00 | – | 0.70 |
| 90–95 | 1.00 | 1.00 | 0.98 | 0.81 | 0.10 | 0.00 | – | 0.64 |
| 85–90 | 1.00 | 0.98 | 0.92 | 0.54 | 0.01 | – | – | 0.68 |
| 80–85 | 1.00 | 0.92 | 0.76 | 0.27 | 0.00 | – | – | 0.58 |
| 70–80 | 0.92 | 0.77 | 0.50 | 0.07 | – | – | – | 0.55 |
| <70 | 0.44 | 0.31 | 0.16 | – | – | – | – | 0.30 |

Row totals depend on the bucket mix and are not meaningful on their own.

**Recall of either test:**

| A \ G | <40 | 40–50 | 50–60 | 60–70 | 70–80 | 80–90 | 90–95 |
|---|---|---|---|---|---|---|---|
| 98–100 | 1.00 | 1.00 | 1.00 | 0.99 | 0.60 | 0.08 | 0.01 |
| 95–98 | 1.00 | 1.00 | 1.00 | 0.98 | 0.43 | 0.11 | – |
| 90–95 | 1.00 | 1.00 | 0.98 | 0.85 | 0.36 | 0.10 | – |
| 85–90 | 1.00 | 0.99 | 0.96 | 0.72 | 0.34 | – | – |
| 80–85 | 1.00 | 0.96 | 0.91 | 0.66 | 0.34 | – | – |
| 70–80 | 1.00 | 0.95 | 0.86 | 0.68 | – | – | – |
| <70 | 0.95 | 0.90 | 0.91 | – | – | – | – |

The high-identity test alone masks 0.44–0.80 of plants per row. The donor's own copy is removed, but other members of its lineage usually remain, so the planted copy's best matches are genomes the focal genome is only distantly related to.

**Recall by realised deficit** (mean of expected − observed identity over the plant's comparisons, all buckets pooled):

| Deficit | <0 | 0–10 | 10–15 | 15–20 | 20–25 | 25–30 | 30–35 | 35–40 | 40–50 | ≥50 |
|---|---|---|---|---|---|---|---|---|---|---|
| Plants | 36 | 1,231 | 1,322 | 1,756 | 1,645 | 1,837 | 1,638 | 1,673 | 2,747 | 2,668 |
| Recall | 0.00 | 0.00 | 0.00 | 0.003 | 0.06 | 0.85 | 0.98 | 0.996 | 0.999 | 1.00 |

Detection is a near step function of the realised deficit at T = 25. The grid's pattern is therefore a statement about which plants end up with small deficits:
- **Diagonal:** recall falls from about 1 to 0 over 10–20 points of G, with the edge at roughly G = A − 20 to A − 30.
- **Distant relatives:** for A < 70 the trend already expects low identity, so even G < 40 gives a small deficit (recall 0.16–0.44).
- **Locus conservation:** within a bucket, recall at the most variable third of loci (median neighbour identity 80–93%) is far lower than at the most conserved third (95–99%). Examples: A 80–85, G 50–60: 0.41 vs 0.98; A 85–90, G 60–70: 0.14 vs 0.84; A < 70, G < 40: 0.20 vs 0.64.
- **G as the axis:** binning by the planted copy's mean identity to the actual neighbours, instead of to the native copy, gives an almost identical grid.

**Other thresholds** (recomputed from the same runs; expected values don't depend on T). Lowering T from 25 to 20 moves the detection edge about 5 points closer to the diagonal:

| Bucket | T = 15 | T = 20 | T = 25 | T = 30 |
|---|---|---|---|---|
| A 98–100, G 70–80 | 0.99 | 0.89 | 0.49 | 0.00 |
| A 90–95, G 70–80 | 0.87 | 0.49 | 0.10 | 0.00 |
| A 85–90, G 60–70 | 0.93 | 0.81 | 0.54 | 0.16 |
| A 80–85, G 60–70 | 0.77 | 0.59 | 0.27 | 0.02 |
| A 70–80, G 50–60 | 0.84 | 0.70 | 0.50 | 0.18 |
| A < 70, G < 40 | 0.82 | 0.63 | 0.44 | 0.25 |
| Masks on the original data (section 4.5) | 4,876 | 2,529 | 1,604 | 1,151 |

**Knock-on masks** (non-planted masks in the planted run minus those in its control): +255 low-identity, +69 high-identity, +338 either over 17,100 plants, i.e. **+19.8 per 1,000 plants**. Per run the net ranges from −3 to +66. Gross changes are larger and go both ways: a planted run gains about 90 and loses about 45 masks against its control, concentrated at 2934341at4751, the gene with the most masks, so most of this is borderline calls flipping. Masks tied to a plant (the planted genome is an outlier neighbour), by the donor lineage's sampling depth (other genomes with AAI ≥ 90 to the donor):

| Donor relatives | Plants | Low / 1,000 plants | High / 1,000 plants |
|---|---|---|---|
| 0 | 5,428 | 27.8 | 11.6 |
| 1 | 2,146 | 12.1 | 9.8 |
| 2–4 | 2,932 | 18.4 | 2.0 |
| 5+ | 6,047 | 11.6 | 1.7 |

772 low-identity and 35 high-identity gained masks could not be tied to a plant. As expected, high-identity knock-on nearly vanishes for donor lineages sampled two or more deep.

**Hardest loci:** pooled low-identity recall is lowest at 333279at4751 (0.45), 351421at4751, 3500106at4751, 3487871at4751 and 3453803at4751 (0.51–0.52). These are variable loci whose plants fall mostly in the variable stratum.

**Limits:** see section 3.3. Plants model contamination and cross-lineage transfer, not paralogs, misalignment or partial sequences. Recall is per bucket, so overall recall on real data depends on the real mix of foreign copies, which this design doesn't estimate.

### 4.11 Bucketed evaluation, Experiment 2: high-identity test

Design in section 3.3. Default settings. 12,500 plants in 14 runs, each with a control; all 25 buckets were filled to 500. 62 plants were excluded (38 whose reported best match fell outside the bucket ranges, 24 with a planted neighbour), leaving 12,438. Plant generation took 21 minutes; the 28 gtrim runs and the analysis took 10. Full report: `$WORK/buckets_high/report.txt`; per-bucket CIs: `$WORK/buckets_high/results.tsv`.

**Plants per bucket:** 474–535, from 448–500 distinct focal genomes and 70–77 loci per bucket; CI half-widths about ±0.04. Holdouts per plant (median, 90th percentile) by I bin: <80: 7, 48; 80–90: 5, 33; 90–95: 3, 19; 95–98: 1, 7; 98–100: 1, 1. Low-I buckets were filled from sparsely sampled donor lineages (median donor relatives 1 for I < 80, 7 for I ≥ 95).

**Recall of the high-identity test** (rows B = AAI of the focal genome to the planted copy's reported best match; columns I = identity of that match):

| B \ I | <80 | 80–90 | 90–95 | 95–98 | 98–100 | Row |
|---|---|---|---|---|---|---|
| <50 | 0.70 | 0.95 | 0.99 | 0.99 | 0.97 | 0.92 |
| 50–60 | 0.52 | 0.88 | 0.96 | 0.98 | 0.99 | 0.86 |
| 60–70 | 0.02 | 0.12 | 0.41 | 0.72 | 0.87 | 0.43 |
| 70–80 | 0.00 | 0.00 | 0.00 | 0.06 | 0.21 | 0.05 |
| 80–90 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

**Recall of either test:**

| B \ I | <80 | 80–90 | 90–95 | 95–98 | 98–100 |
|---|---|---|---|---|---|
| <50 | 0.97 | 0.99 | 1.00 | 1.00 | 1.00 |
| 50–60 | 0.90 | 0.94 | 0.97 | 0.99 | 1.00 |
| 60–70 | 0.65 | 0.52 | 0.63 | 0.82 | 0.92 |
| 70–80 | 0.34 | 0.18 | 0.16 | 0.23 | 0.33 |
| 80–90 | 0.26 | 0.06 | 0.05 | 0.04 | 0.03 |

The low-identity test alone masks 0.81–0.89 of plants when B < 60, 0.47–0.64 at B 60–70 and 0.03–0.34 above that. Its recall is highest in the I < 80 column, probably because there the best remaining match is often outside the donor's lineage, so the donor is more distant than B suggests.

**Recall by realised AAI deficit:** none of the 4,099 plants below a deficit of 15 and about 1 of 1,316 at 15–20; 0.045 at 20–25; 0.91 at 25–30; over 0.99 at 30 and above. As in Experiment 1, detection is a step function of the deficit at the threshold.

**The detection ceiling is structural.** A comparison is an outlier only if its AAI is more than T′ below the AAI expected at its locus identity. Expected AAI of the rank-1 match in the baseline run, by identity:

| I | <80 | 80–90 | 90–95 | 95–98 | 98–100 |
|---|---|---|---|---|---|
| Median expected AAI (IQR) | 78.5 (73–83) | 85.9 (83–90) | 91.0 (88–93) | 94.9 (93–96) | 97.9 (97–99) |
| B must be below about | 53 | 61 | 66 | 70 | 73 |

These limits match where recall falls in each column. At T′ = 25, copies from lineages above about 70% AAI to the focal genome look like ordinary relatives.

**Locus conservation** has the opposite effect from Experiment 1: weak matches are harder at conserved loci, where only distant genomes fall below 90% identity, so the expected AAI is already low. Examples (variable / middle / conserved): B < 50, I < 80: 0.85 / 0.60 / 0.47; B 50–60, I < 80: 0.72 / 0.35 / 0.18; B 60–70, I 90–95: 0.69 / 0.43 / 0.19. The effect fades for near-identical matches (B 60–70, I 98–100: 0.92 / 0.83 / 0.87).

**Other thresholds** (recomputed):

| Bucket | T′ = 15 | T′ = 20 | T′ = 25 | T′ = 30 |
|---|---|---|---|---|
| B 50–60, I < 80 | 0.89 | 0.75 | 0.52 | 0.25 |
| B 60–70, I 90–95 | 0.94 | 0.80 | 0.41 | 0.11 |
| B 60–70, I 95–98 | 0.98 | 0.92 | 0.72 | 0.30 |
| B 70–80, I 95–98 | 0.73 | 0.41 | 0.06 | 0.00 |
| B 70–80, I 98–100 | 0.87 | 0.62 | 0.21 | 0.00 |
| B 80–90, I 98–100 | 0.38 | 0.00 | 0.00 | 0.00 |
| Masks on control_nodonor (section 4.9) | 1,027 | 269 | 107 | 57 |

**Knock-on masks:** high-identity −13 over 12,500 plants (−1.0 per 1,000), low-identity +140 (+11.2), either +133 (+10.6); per run −13 to +37. Gross churn against the controls is about 58 gained and 49 lost masks per run, led again by 2934341at4751. Masks tied to a plant fall with donor lineage depth (high-identity: 3.8, 2.7, 1.8 and 0.5 per 1,000 plants for 0, 1, 2–4 and 5+ relatives).

**Hardest loci:** pooled high-identity recall is lowest at 101857at4751, 3470561at4751, 94760at4751, 3494367at4751 and 3483416at4751 (0.37–0.39); pooled recall depends on each locus's bucket mix.

**Limits:** every plant keeps a match in the data by construction (copies whose source lineage is absent are left to the low-identity test). Weak-match buckets were filled from sparsely sampled donor lineages with more holdouts. Recall is per bucket.

### 4.12 Review follow-ups

These analyses followed a review of the method (2026-10-23). They use the current code, default settings and the original example data. Masks at T = 20 are recomputed from the comparisons tables, since expected values don't depend on T; the recomputation gives 1,692 masks at T = 25 against gtrim's 1,693, because one comparison's deficit equals T to the 4 decimals written in the TSV.

**Tied x values in the running median (code change).** Previously each point's window was the *N* nearest points in sorted order, so where many points share the same x, which tied points fell in the window depended on genome order. In the high-identity test, 21% of comparisons sit in runs of more than 120 identical locus identities (mostly 100%), and every point in such a run got the median AAI of the first 120 genomes in accession order. That differed from the median over the whole run by at most 2.1 AAI points. The window now includes every point tied with the x values at its ends (section 2.4). Effect on the example data:
- **Low-identity test:** outlier comparisons unchanged (6,542).
- **High-identity test:** outlier comparisons 690 → 683, masks 112 → 111.
- **Total:** 1,695 → 1,693 masked genome loci. Two masks within 0.1 point of the threshold were removed and none added.
- **Integration test:** outputs unchanged.

**Suspect genomes (new output).** `divergent_loci_genome_summary.tsv` has a `suspect` column, true when at least `--suspect-masked-perc` (default 10) percent of a genome's loci, and at least 3 loci, are masked. 10% is 25 times the overall masked rate (0.40%). On the example data, 11 genomes are flagged at the default. At 15%, the flag would keep only the 3 clearest cases: GCA_023085455.1 (29 of 44 loci masked), GCA_029303655.1 (10 of 32) and GCA_029300685.1 (7 of 29). The flag only reports; nothing is dropped.

**Per-genome effect of lowering T and T′ from 25 to 20.** Masks rise from 1,692 to 2,765 (low-identity 1,602 → 2,530, high-identity 111 → 286).
- **Spread:** 845 of 5,813 genomes gain at least one masked locus, and 690 of them gain exactly one.
- **Concentrated in suspect genomes:** the largest gains are in genomes already flagged. GCA_029300685.1 gains 13 high-identity masks (7 → 20 of its 29 loci), and GCA_023085455.1 gains 7 (29 → 36 of 44).

| Genome | Loci present | Masked T=25 | Masked T=20 | Added | Added by low-identity | Added by high-identity | Loci left T=25 | Loci left T=20 |
|---|---|---|---|---|---|---|---|---|
| GCA_029300685.1 | 29 | 7 | 20 | +13 | +1 | +13 | 22 (29%) | 9 (12%) |
| GCA_023085455.1 | 44 | 29 | 36 | +7 | +7 | +0 | 15 (19%) | 8 (10%) |
| GCA_000760515.2 | 60 | 5 | 12 | +7 | +6 | +3 | 55 (71%) | 48 (62%) |
| GCA_022605425.2 | 43 | 3 | 10 | +7 | +3 | +5 | 40 (52%) | 33 (43%) |
| GCA_018985225.1 | 58 | 0 | 6 | +6 | +3 | +4 | 58 (75%) | 52 (68%) |
| GCA_003600395.1 | 35 | 5 | 10 | +5 | +5 | +0 | 30 (39%) | 25 (32%) |
| GCA_001875675.1 | 47 | 4 | 9 | +5 | +3 | +4 | 43 (56%) | 38 (49%) |
| GCA_026539955.1 | 49 | 4 | 9 | +5 | +4 | +1 | 45 (58%) | 40 (52%) |
| GCA_022664585.2 | 76 | 0 | 5 | +5 | +5 | +0 | 76 (99%) | 71 (92%) |
| GCA_029303655.1 | 32 | 10 | 14 | +4 | +0 | +4 | 22 (29%) | 18 (23%) |
| GCA_002087915.1 | 34 | 0 | 4 | +4 | +2 | +2 | 34 (44%) | 30 (39%) |
| GCA_010405405.1 | 77 | 0 | 4 | +4 | +3 | +1 | 77 (100%) | 73 (95%) |
| GCA_026539545.1 | 56 | 2 | 5 | +3 | +3 | +0 | 54 (70%) | 51 (66%) |
| GCA_000385875.2 | 38 | 1 | 4 | +3 | +1 | +2 | 37 (48%) | 34 (44%) |
| GCA_001432165.1 | 33 | 1 | 4 | +3 | +0 | +3 | 32 (42%) | 29 (38%) |
| GCA_001661335.1 | 76 | 1 | 4 | +3 | +1 | +2 | 75 (97%) | 72 (94%) |
| GCA_024243835.1 | 40 | 1 | 4 | +3 | +2 | +1 | 39 (51%) | 36 (47%) |
| GCA_025603075.1 | 73 | 1 | 4 | +3 | +0 | +3 | 72 (94%) | 69 (90%) |
| GCA_026025075.1 | 76 | 1 | 4 | +3 | +3 | +0 | 75 (97%) | 72 (94%) |
| GCA_027600645.1 | 63 | 1 | 4 | +3 | +0 | +3 | 62 (81%) | 59 (77%) |

**Genomes dropping below 30% of the 77 loci** (fewer than 24 loci kept). All 5,813 genomes have at least 24 loci before masking.
- **T = T′ = 25:** 7 genomes drop below 30%: GCA_011428085.1 (25 → 22), GCA_023085455.1 (44 → 15), GCA_029300685.1 (29 → 22), GCA_029303655.1 (32 → 22), GCA_036937195.1 (24 → 22), GCA_910592325.1 (26 → 23), GCA_958449995.1 (24 → 22).
- **T = T′ = 20:** 9 genomes drop below 30%, adding GCA_000383075.1 and GCA_964256305.1 (each 26 → 23).

**Divergent copies shared by close relatives.** Both tests compare a genome with its 3 nearest relatives, so a copy shared by a genome and 2 or more of them passes. Such sharing can come from a transfer into a common ancestor, a contaminant common to several assemblies, or a lineage carrying a different paralog. To look for it, gtrim was run with 10 neighbours in both tests (`--neighbours 10 --upper-neighbours 10`). A candidate is a genome locus with at most 1 outlier among its 3 nearest neighbours and at least 4 among the next 7. Candidates at the same locus that are among each other's 3 nearest neighbours were grouped (`analyse_shared_divergence.py`).

| Test | Groups of ≥ 3 genomes | Genome loci | Loci | Already masked at k = 3 | Clean groups (mean deficit to the 3 nearest < 10) |
|---|---|---|---|---|---|
| Low-identity | 96 | 424 (0.10%) | 30 | 6 | 26 groups, 108 genome loci, 17 loci |
| High-identity | 6 | 28 (0.01%) | 2 | 1 | 2 groups, 9 genome loci |

- **Shared low identity exists, mostly in problem loci.** The groups concentrate in loci that already have many low-identity masks: 2934341at4751 (18 groups), 3494367at4751 (11), 3505740at4751 (8), 3472996at4751 (7), 520361at4751 (7) and 3323237at4751 (6). This fits a lineage carrying a different member of the gene family (hidden paralogy) better than lateral transfer. A clear example: at 3479320at4751, GCA_007182795.1 agrees with its 3 nearest relatives (95.6–97.6% AAI, deficits under 3 points), but sits about 30 points below the trend against 6 of the next 7, which are almost as closely related (93.3–95.2% AAI).
- **Shared copies from sampled lineages are rare.** The high-identity version of the pattern, a group whose locus is closest to an unrelated lineage, gives 28 genome loci, 5 of the 6 groups at 520361at4751.
- **Limits:** with 10 neighbours, only groups of 3–10 genomes are visible. The pattern doesn't say which side is divergent (the group or the more distant genomes). The trend is fitted over all 10 ranks.

**Assessment:** a dedicated clade-level test isn't justified by these numbers. Shared divergence affects about 0.1% of genome loci, a quarter as many as the k = 3 masks, and sits mostly in a handful of loci that are already problematic. Reviewing or dropping those loci deals with most of it.

### 4.13 Minimum aligned columns

`--min-aligned-columns` (*m*, default 50) excludes a neighbour from both tests when fewer than *m* columns are aligned at the locus. `--min-aligned-frac` alone is relative: a 20-residue fragment passes it with 10 aligned columns.

**The example data has no short overlaps.** Every comparison has at least 50 aligned columns, so the data can't show short-overlap failures directly. Residual scatter does widen as overlap shrinks (default run, low-identity test; `analyse_min_aligned.py`):

| Aligned | Comparisons | Outlier % | Robust SD of residual | Genome loci | Masked % |
|---|---|---|---|---|---|
| 50-74 | 103 | 0.00 | 8.8 | 55 | 0.00 |
| 75-99 | 143 | 1.40 | 4.8 | 68 | 1.47 |
| 100-149 | 2,278 | 0.75 | 5.3 | 1,156 | 0.43 |
| 150-199 | 73,540 | 0.20 | 3.1 | 25,018 | 0.15 |
| 200-299 | 216,216 | 0.24 | 2.2 | 72,807 | 0.19 |
| 300+ | 974,731 | 0.60 | 2.0 | 323,233 | 0.44 |

**Sampling noise sets the value.** For a vertical copy compared over *n* columns, identity is binomial. The table gives the chance that noise alone makes a comparison an outlier at T = 25, and the expected number of unrelated genomes that reach ≥ 90% identity by chance among 5,000 candidates. The second is what fools the high-identity test, which picks neighbours by highest identity:

| n | Low-identity: P(outlier), expected identity 50 / 70 / 90% | Low-identity: P(2 of 3) at 70% | High-identity: chance matches >= 90% identity among 5,000 genomes at 50 / 60 / 70% identity |
|---|---|---|---|
| 10 | 5.5e-02 / 4.7e-02 / 1.3e-02 | 7e-03 | 54 / 2.3e+02 / 7.5e+02 |
| 20 | 5.9e-03 / 5.1e-03 / 4.2e-04 | 8e-05 | 1 / 18 / 1.8e+02 |
| 30 | 2.6e-03 / 2.1e-03 / 8.9e-05 | 1e-05 | 0.021 / 1.6 / 47 |
| 40 | 3.4e-04 / 2.7e-04 / 3.5e-06 | 2e-07 | 0.00046 / 0.14 / 13 |
| 50 | 1.5e-04 / 1.2e-04 / 7.6e-07 | 4e-08 | 1.1e-05 / 0.013 / 3.6 |
| 75 | 3.6e-06 / 2.5e-06 / 1.4e-09 | 2e-11 | 2.9e-10 / 1.6e-05 / 0.082 |
| 100 | 9.1e-08 / 5.9e-08 / 2.9e-12 | 1e-14 | 0 / 1.2e-07 / 0.0078 |

- **Low-identity test:** the 2 of 3 vote already protects it. At *n* = 20, about 1 in 10,000 vertical genome loci would be masked by noise; at 50, about 1 in 10 million.
- **High-identity test:** this is the tighter constraint. At *n* = 20, about 18 of 5,000 genomes sharing 60% identity reach ≥ 90% by chance and could become the "closest" matches. At 50, the expected number is 0.01.

**Effect on the example data** (default settings otherwise): *m* = 50 gives output identical to no minimum (1,693 masks). *m* = 100 leaves 47 genome loci unevaluated, with masks unchanged. *m* = 150 leaves 523 unevaluated and loses 2 masks. 50 is the smallest round value at which both tests are safe from sampling noise, and it costs nothing here.

### 4.14 Source of masked copies

Default run on the example data (no taxonomy is available for these genomes). The new columns don't change which loci are masked; runtime (33 s on 96 cores) and memory (1.2 GB) are unchanged.
- **High-identity masks (111):** `source_coherent` is true for 67 and false for 44, so about 60% point to a single source lineage. All 21 high-identity masks of GCA_023085455.1 have the same best match, GCA_004016305.1, the signature of contamination from one organism.
- **Low-identity-only masks (1,582):** best-match identity is below 50% for 201 (no similar copy anywhere: paralog, misalignment or unsampled source), 50–80% for 618, 80–95% for 464 and ≥ 95% for 299. In the last group the best match is a close relative (median leave-one-out AAI 96.9%) that shares the odd copy, not a foreign lineage.
- **Per genome:** 60 genomes have high-identity masks, 12 of them two or more:

| Genome | Loci | High-identity masks | Main source | Masks from its lineage | Suspect |
|---|---|---|---|---|---|
| GCA_023085455.1 | 44 | 21 | GCA_004016305.1 | 21 | true |
| GCA_029303655.1 | 32 | 9 | GCA_016859295.1 | 7 | true |
| GCA_004114315.1 | 69 | 5 | GCA_036320915.1 | 4 | false |
| GCA_029300685.1 | 29 | 5 | GCA_016859295.1 | 2 | true |
| GCA_004114245.1 | 58 | 4 | GCA_000507065.1 | 3 | true |
| GCA_011033665.1 | 68 | 4 | GCA_011426335.1 | 4 | false |
| GCA_024708745.1 | 30 | 4 | GCA_000256225.1 | 4 | true |
| GCA_029300755.1 | 36 | 3 | GCA_016859295.1 | 3 | true |
| GCA_002994575.1 | 44 | 2 | GCA_000149925.1 | 1 | false |
| GCA_006912115.1 | 64 | 2 | GCA_000517465.1 | 2 | false |
| GCA_041381035.1 | 68 | 2 | GCA_039778275.1 | 1 | false |
| GCA_910592055.1 | 61 | 2 | GCA_003073295.1 | 2 | false |

  Most genomes with several high-identity masks have one main source lineage, consistent with contamination by one organism. Three suspect genomes with neighbouring accessions (GCA_029300685.1, GCA_029300755.1 and GCA_029303655.1) share the same main source, GCA_016859295.1, which suggests the same contaminant across a batch of assemblies.

### 4.15 Marker gene review: paralogy and transfer per gene

Question: should some of the 77 markers be dropped from phylogenetic inference? Masking handles isolated problems (one genome's copy unlike its nearest relatives' copies) but cannot see a problem shared by a whole lineage, whose members are each other's nearest relatives. `gene_quality.py` combines two kinds of evidence per gene:
- **Masks** from the default run on the unmodified data (low- and high-identity tests); low-identity masks with mean identity < 50% to the neighbours are counted as severe (paralog or misalignment-like).
- **Pairwise discordance** across the whole AAI range: 99,197 genome pairs sampled evenly over 8 AAI strata; for each gene, identity is compared with the gene's own median identity at the pair's leave-one-out AAI, and a pair is discordant if it is more than 25 points off. Counted among pairs with leave-one-out AAI ≥ 60, before and after removing masked copies. **Problem genomes** have ≥ 4 sampled pairs after masking, at least half discordant; they are grouped by single linkage at 80% and 65% AAI (no taxonomy file is available). Each gene file has one sequence per genome, so extra copies can't be seen directly.

Results (identical ranking with a second random sample of pairs):
- The median gene has 0.13% of genomes masked and 0.009% discordant pairs after masking; 56 genes are under 0.1% and 6 are above 0.5%.

| Gene | Masked | Severe | High-identity masks | Discordant before → after masking | Problem genomes (largest lineage at 65%) | Within lineage / vs others | Reading |
|---|---|---|---|---|---|---|---|
| 2934341at4751 | 9.08% | 369 | 2 | 14.8% → 3.0% | 56 (53) | 0 / −60 | lineage carries a different gene |
| 520361at4751 | 1.79% | 10 | 33 | 2.3% → 2.1% | 427 (427) | −31 / −28 | large lineage diverges inconsistently |
| 3472996at4751 | 0.76% | 9 | 5 | 3.2% → 3.0% (1.8% upward) | 45 (23) | +27 / +25 | copies too similar across lineages: transfer |
| 3505740at4751 | 0.86% | 18 | 0 | 1.1% → 0.9% | 141 (141) | −27 / −25 | as 520361at4751, smaller lineage |
| 3494367at4751 | 1.12% | 16 | 3 | 1.5% → 1.2% (0.7% upward) | 28 (19) | +33 / +23 | transfer-like |
| 3453803at4751 | 1.71% | 89 | 1 | 1.9% → 0.7% | 16 (12) | −5 / −27 | small lineage with a different copy |
| 3323237at4751 | 2.30% | 115 | 1 | 3.1% → 0.3% | 15 (9) | −4 / −45 | small lineage with a different copy |
| 3496677at4751 | 0.15% | 5 | 3 | 0.3% → 0.3% | 35 (29) | −26 / −26 | lineage-wide divergence |

Within / vs others: median residual (points) for pairs among the problem genomes and between them and other genomes.

**Recommendations.**
- **Remove** 2934341at4751, 520361at4751 and 3472996at4751: their signals are lineage-wide or cross-lineage, so masking can't fix them. Cost: 2,173 of 36,810 columns (5.9%), leaving 74 genes.
- **Review** with gene trees: 3505740at4751, 3494367at4751, 3453803at4751 (another 2,192 columns; 11.9% in total with the first three).
- **Keep** the other 71 with masking (discordance after masking ≤ 0.4%, no large problem lineage). For 3323237at4751 and 3453803at4751 (and 2934341at4751 if kept), masking the gene in the affected lineage only is an alternative to dropping it; gtrim doesn't support this yet.

**Order of the tests and repeated masking** (runs in `$WORK/sequential/`). The two tests are evaluated together and independently on the unmodified data, and their masks are combined. Running one test on the other's masked output instead changes only borderline calls:

| Variant | Masks | Added vs together | Dropped |
|---|---|---|---|
| Together (current; this binary) | 1,693 | | |
| Low-identity first, then high-identity on its output | 1,675 | 2 | 20 |
| High-identity first, then low-identity on its output | 1,678 | 7 | 22 |
| Together, repeated on the masked output until stable | 1,838 | 145 | 0 |

The dropped masks were borderline (median second-largest deficit 28.8 → 22.9 for the 20 high-identity masks) and only 1 of the 42 had an outlier neighbour masked by the other test, so sequencing doesn't correct blame errors; it moves borderline calls. Repeating adds 122, 15 and 8 masks in rounds 2–4. These concentrate in the flagged genes: the 6 remove or review genes take 67 of the 145 (46%), led by 520361at4751 (24) and 2934341at4751 (19); the other 78 are spread over 36 genes (94760at4751 has 9). The mechanism is lineage erosion: genomes whose neighbours share a divergent copy look normal until those neighbours are masked, after which normal copies replace them as neighbours (median second-largest deficit 22.0 → 26.8). This supports the gene review: repeated masking chips at lineage-wide problems without clearing them. Simultaneous, single-pass masking is kept as the default.

**Caveats.** Identity alone can't separate a lineage-wide paralog from a strong lineage-specific rate change; gene trees for the flagged genes, named taxonomy for the problem lineages, and concatenated trees with and without the three genes (topology and support) are the next steps.

### 4.16 Non-decreasing trend lines (default method since 2026-10-23)

The expected value should not fall as relatedness rises: expected identity should be non-decreasing in leave-one-out AAI (low-identity test), and expected AAI non-decreasing in locus identity (high-identity test). The running median dips locally, by a median of 2.7–2.8 points per gene and up to 7.4 (low-identity, 97101at4751) or 15.2 (high-identity, 2934341at4751).

**What-if** (`monotone_trend_whatif.py`, on the baseline comparison tables): replacing each locus's running medians by their isotonic fit and re-applying the threshold and vote changed 46 low-identity and 10 high-identity masks. The median gene's expected values moved by only about 0.2 points; the changes concentrate in the flagged genes (3505740at4751, 520361at4751, 3496677at4751, 3472996at4751, 3494367at4751, 2934341at4751), where problem lineages pull the median down within their AAI range.

**Implemented** as `--method monotone-median`, now the default for both tests (section 2.4); `--method median` keeps the plain running median. Full example data, current binary:

| Method | Masked | Low-identity | High-identity | Added vs median | Removed |
|---|---|---|---|---|---|
| `median` | 1,693 | 1,603 | 111 | | |
| `monotone-median` (default) | 1,691 | 1,599 | 109 | 23 | 25 |

The run took 31 s and 1.2 GB on 96 threads, as before. The integration test outputs are unchanged. Results in sections 4.1–4.15 were produced with `median`; the bucketed experiments have not been rerun with the new default.

## 5. Conclusions so far

1. **Use k = 3 by default.**
   - **k = 1:** never use it. It can't tell which genome in a pair is divergent: about 27% of its masks on the real data are knock-on masks, and it adds about one false mask per planted locus.
   - **k = 3:** among settings without knock-on masks, it has the highest recall and recovers 99% of the consensus.
   - **Conservative alternatives:** k = 2 (2 of 2) and k = 4 (3 of 4) give the most conservative masked sets, with 95–98% consensus precision, at a small cost in recall.
2. **The window size is not critical between 40 and 200.** N = 120 is on the plateau of mutual agreement. Avoid N = 400.
3. **The outlier threshold T is the main lever, and is a choice between sensitivity and how much gets masked.**
   - **What T determines:** a planted locus is detected reliably once its identity falls about T + 10 points below its neighbours' AAI. Masked sets are nested across T.
   - **T = 25** (current default) masks 0.38% of genome loci and finds 55% of the planted loci, including all with a gap of 40 points or more.
   - **T = 20** masks 58% more (0.60%) and raises planted recall to 0.66, with a similar or lower estimated noise fraction (≤ 12% vs ≤ 15%). It is a reasonable alternative if moderately divergent copies matter.
   - **T = 15** triples the masked set for a recall of 0.78.
   - **T = 10** masks 2.8% of genome loci with up to about 26% noise, which is too permissive.
   - **Decision (2026-10-23): keep T = T′ = 25.** False positives cost more than masking 1/77 of a complete genome's data suggests: (1) a partial genome has fewer loci to lose, and (2) removing the most divergent genuine genes shortens the branch lengths inferred in the tree. Section 4.12 shows the per-genome effect of T = 20 that was weighed.
4. **The high-identity test is complementary and on by default.** It finds copies acquired from lineages present in the data (contamination, lateral transfer), including many the low-identity test misses. Together the two tests raise planted recall from 0.54 to 0.70 with no knock-on masks. The high-identity test has its own parameters (*k*′, *N*′ and *T*′). The sweep in section 4.9 supports keeping the defaults (3, 120, 25):
   - **Neighbours:** k′ = 1 produces knock-on masks, and k′ ≥ 4 loses recall.
   - **Window:** N′ hardly matters.
   - **Threshold:** T′ = 25 isolates clear contamination-like cases. T′ = 20 raises combined recall from 0.70 to 0.80 and still masks only 269 genome loci on the control (0.07%). It's a reasonable choice when sensitivity matters more, but the added masks are less clearly foreign copies. T′ ≤ 15 is too permissive.

5. **The low-identity test does what its threshold says (Experiment 1, section 4.10).** A foreign copy is masked when its realised deficit exceeds T: recall is 0.85 at deficits of 25–30 and over 0.99 beyond 30. Its blind spots are copies within about 20 points of the focal lineage, genomes whose nearest relatives are below 70% AAI, and fast-evolving loci. The high-identity test recovers much of that gap when the source lineage is in the data. Knock-on costs about 2 masks per 100 plants.

6. **The high-identity test also does what its threshold says (Experiment 2, section 4.11).** It masks a borrowed copy when the best match's AAI is at least T′ below the AAI expected at that identity: recall is 0.88–0.99 for source lineages below 60% AAI with matches of 80% identity or more. At T′ = 25 it cannot see copies from lineages above about 70% AAI (they look like ordinary relatives); T′ = 20 raises recall near that edge by 0.2–0.4. It causes no knock-on masks. With the low-identity test, recall is 0.90–1.00 when the source lineage is below 60% AAI.

7. **Divergent copies shared by groups of close relatives are rare and concentrated (section 4.12).** About 0.1% of genome loci show the pattern, mostly at loci that already have many masks, and copies shared from sampled foreign lineages are rarer still (28 genome loci). Handling the few problematic loci is enough; a clade-level test is not needed. Lowering T and T′ to 20 adds 1,073 masks, spread thinly (690 genomes gain one); the largest gains are in genomes already flagged as suspect.

8. **Drop three markers (section 4.15).** 2934341at4751, 520361at4751 and 3472996at4751 carry lineage-wide or cross-lineage paralogy or transfer signals that masking cannot fix (5.9% of columns). 3505740at4751, 3494367at4751 and 3453803at4751 need a gene-tree check. The other 71 genes are sound once masked.

9. **Use a non-decreasing trend (section 4.16).** Making each locus's running median non-decreasing changes about 3% of masks, mostly in the flagged genes, and removes trend dips that have no biological basis. It is the default (`--method monotone-median`).

## 6. Reproducing the analysis

The scripts are in `analysis/mask_divergent_loci/`. Results are written to a working directory given by `WORK`, which shouldn't be inside the repository. The runs need a release build (`cargo build --release`).

```bash
export WORK=/path/to/work_dir
A=analysis/mask_divergent_loci
DATA=/srv/projects/donovan/gtrim/data/gtdb_r232_fun77_markers

# planted data: writes $WORK/planted_msas/ and $WORK/planted_loci.tsv
python3 $A/plant.py $DATA $WORK/planted_msas 7

# planted data with donor genomes removed, plus a control (original data minus donors):
# writes $WORK/planted_nodonor_msas/, $WORK/planted_nodonor_loci.tsv and $WORK/control_nodonor_msas/
python3 $A/plant.py $DATA $WORK/planted_nodonor_msas 7 $WORK/control_nodonor_msas

# sweeps of the low-identity test (run_job.sh passes --no-upper-mask unless UPPER=1)
# each line of a job file is "<dataset> <k> <N> <T>"; completed runs are skipped
xargs -P 4 -L 1 $A/run_job.sh < $A/jobs_sweep1.txt
xargs -P 4 -L 1 $A/run_job.sh < $A/jobs_sweep2.txt

# runs that keep the full comparisons table (used for the gap, k=1 and mirror analyses; written to $WORK/details/)
KEEP=1 $A/run_job.sh orig 1 120 25
KEEP=1 $A/run_job.sh orig 3 120 25
KEEP=1 $A/run_job.sh planted 1 120 25
KEEP=1 $A/run_job.sh planted 3 120 25

# high-identity test evaluation (section 4.8)
UPPER=1 KEEP=1 $A/run_job.sh control_nodonor 3 120 25
UPPER=1 KEEP=1 $A/run_job.sh planted_nodonor 3 120 25

python3 $A/analyse_sweep1.py
python3 $A/analyse_sweep2.py
python3 $A/knock_on_k1.py
python3 $A/analyse_scale_threshold.py   # data-driven thresholds (section 4.7); about 2 minutes
python3 $A/analyse_upper.py             # high-identity test (section 4.8)

# high-identity test parameter sweep (section 4.9): the low-identity test stays at 3/120/25,
# and k, N and T apply to the high-identity test; results go to $WORK/runs/<dataset>_uk<k>_un<N>[_ut<T>]
UPPER_SWEEP=1 xargs -P 4 -L 1 $A/run_job.sh < $A/jobs_upper_sweep.txt
python3 $A/analyse_upper_sweep.py
```

**Bucketed evaluation (section 3.3).** Run these from the repository root. The timings are for 96 cores.

```bash
export WORK=/path/to/bucket_eval        # about 11 GB
export DATA=/srv/projects/donovan/gtrim/data/gtdb_r232_fun77_markers
A=analysis/mask_divergent_loci

cargo build --release --bins --examples                               # gtrim and examples/aai_matrix.rs
target/release/examples/aai_matrix $DATA $WORK/aai 96                  # genome AAI matrix (~25 s)
UPPER=1 KEEP=1 THREADS=96 $A/run_job.sh orig 3 120 25                  # baseline run (~30 s)

python3 $A/make_bucket_plants.py --experiment low --msa-dir $DATA --work $WORK --target 500 --target-overrides "<70=200"   # ~9 min
python3 $A/make_bucket_plants.py --experiment high --msa-dir $DATA --work $WORK --target 500                                # ~21 min (was ~19 h before the early exit added 2026-10-22)

for e in low high; do                                                   # ~70 runs, ~30 min in total
  UPPER=1 KEEP=1 xargs -P 4 -L 1 $A/run_job.sh < $WORK/buckets_$e/jobs.txt
  UPPER=1 xargs -P 4 -L 1 $A/run_job.sh < $WORK/buckets_$e/jobs_control.txt
done

python3 $A/analyse_buckets.py --experiment low --work $WORK | tee $WORK/buckets_low/report.txt     # ~1 min each
python3 $A/analyse_buckets.py --experiment high --work $WORK | tee $WORK/buckets_high/report.txt
python3 $A/real_g_distribution.py                                     # real loci by A and G (section 4.10), ~1 min
python3 $A/gene_quality.py --msa-dir $DATA                            # marker gene review (section 4.15), ~30 s
python3 $A/gene_grid_plot.py grid low grid_a.png 1 42                 # low-identity test grids, ranked by its masks (needs matplotlib)
python3 $A/gene_grid_plot.py grid low grid_b.png 43 77
python3 $A/gene_grid_plot.py panels gene_panels/                      # one plot per gene and test plus genes.json, for the scrollable page
# the scrollable page uses the monotone-median run (section 4.16):
python3 $A/gene_quality.py --msa-dir $DATA --run-dir $WORK/monotone/default --out-dir $WORK/gene_quality_monotone
RUN_DIR=$WORK/monotone/default GENE_QUALITY=$WORK/gene_quality_monotone/gene_quality.tsv python3 $A/gene_grid_plot.py panels gene_panels/
```

`$WORK/run_all.sh` runs all of the above unattended and resumably: Experiment 1's runs and analysis first, while Experiment 2's plants are generated in the background, then Experiment 2. Launch it detached so it survives the session: `cd $WORK && setsid nohup ./run_all.sh >> run_all.log 2>&1 < /dev/null &`. Experiment 1's 36 runs took 11.5 minutes and its analysis about 1 minute; Experiment 2's 28 runs and analysis took 10 minutes. The whole driver took 31 minutes.

**Review follow-ups (section 4.12).** These use a separate working directory (about 1.4 GB).

```bash
export WORK=/path/to/review_dir
UPPER=1 KEEP=1 THREADS=96 $A/run_job.sh orig 3 120 25     # default run (~35 s)
UPPER=1 KEEP=1 THREADS=96 $A/run_job.sh orig 10 120 25    # 10 neighbours in both tests (~50 s)
python3 $A/analyse_threshold_genomes.py                   # T = T' = 25 vs 20 per genome; 30% counts
python3 $A/analyse_shared_divergence.py                   # divergent copies shared by close relatives
python3 $A/analyse_min_aligned.py                         # choice of --min-aligned-columns (section 4.13)
```


On 96 cores each run takes about 25 s. Sweep 1 (48 runs) took 21 minutes run one at a time; sweep 2 (96 runs) took 28 minutes running 4 jobs in parallel. `analyse_sweep2.py` takes about 30 s.

## 7. Change log

| Date | Change |
|---|---|
| 2026-10-22 | Initial method (full-AAI neighbours, N = 40). Neighbour selection changed to leave-one-out AAI; the default N changed from 40 to 120. Sweep 1 (k × N), planted-locus evaluation and k = 1 analysis added. |
| 2026-10-22 | Sweep 2 added: k = 2 and 4, outlier threshold T = 10–40, recall by gap, mirror noise estimate. Conclusions updated: T is the main lever; T = 20 is an alternative to the default of 25. |
| 2026-10-22 | Documented how plants are chosen, their divergence distribution and the rationale for it, and the coverage and caveats of the mirror noise estimate (scatter by identity, non-noise above-trend cases, asymmetry). |
| 2026-10-22 | Characterised above-trend outliers (section 4.6): rare, concentrated in a few genomes, full-length, consistent with contamination or chimeras. Two-sided detection raised as an open design question. |
| 2026-10-22 | Explored data-driven thresholds (section 4.7): T = c × local robust SD per locus. Gives higher recall for a similar number of masks; c itself is not determinable from the data alone. |
| 2026-10-22 | Added the high-identity test (sections 2.5, 4.8; on by default, `--no-upper-mask` disables it) and the per-genome summary output. Evaluated on planted loci with donors removed: recall 0.54 → 0.70 with no knock-on masks. `plant.py` can remove donors and write a control dataset; `run_job.sh` disables the new test unless UPPER=1. |
| 2026-10-22 | Removed the ANI (skani) mode: genome similarity is always leave-one-out AAI from the MSAs. Results in this document were all produced in AAI mode and are unaffected. |
| 2026-10-22 | The high-identity test has separate parameters: `--upper-neighbours`, `--upper-window-size` and `--upper-outlier-threshold` (defaults 3, 120 and 25, so results are unchanged). `run_job.sh` with UPPER=1 sets them to the same k, N and T. Chimeric genomes are deliberately not dropped. |
| 2026-10-22 | High-identity test parameter sweep (section 4.9): keep k′ = 3, N′ = 120, T′ = 25. T′ = 20 is a more sensitive option, but its added masks are less clearly foreign copies. `run_job.sh` gained UPPER_SWEEP=1. |
| 2026-10-22 | Designed the bucketed evaluation (section 3.3) and wrote its scripts (`examples/aai_matrix.rs`, `make_bucket_plants.py`, `analyse_buckets.py`). Test runs on the full data led to low planting density, unrelated donors per locus, donor removal and control runs. Not yet run. |
| 2026-10-22 | Review of the bucketed evaluation scripts: knock-on count corrected (planted keys removed from both sides), planted-neighbour exclusion extended to both tests, distinct loci per bucket reported, one plant per focal genome per run. Estimates for Experiment 2 revised. |
| 2026-10-22 | Bucketed evaluation design revised after review: holdout cap 100 with a 2% per-locus budget, per-bucket locus quota, unused genome loci first, neighbour-identity grid, locus-conservation strata, knock-on by donor lineage depth. Full run started (WORK = /srv/projects/donovan/gtrim/bucket_eval). |
| 2026-10-22 | First full attempt stopped during Experiment 2 plant generation (process killed with its session). Found that `make_bucket_plants.py --experiment high` spent ~80 min per run searching donors at loci whose holdout budget was spent; added an early exit (~84 s per run). Relaunched detached, Experiment 1 first. |
| 2026-10-22 | Experiment 1 results (section 4.10) and conclusion 5. |
| 2026-10-22 | Experiment 2 results (section 4.11) and conclusion 6. Both bucketed experiments complete. |
| 2026-10-23 | Review follow-ups (section 4.12, conclusion 7): ties in the running median now share one window that includes every tied point (1,695 → 1,693 masks); `suspect` flag and `--suspect-masked-perc` added to the genome summary; per-genome effect of T = T′ = 20 and genomes dropping below 30% of loci; search for divergent copies shared by close relatives (`analyse_threshold_genomes.py`, `analyse_shared_divergence.py`). |
| 2026-10-23 | Added `--min-aligned-columns` (default 50; section 4.13, `analyse_min_aligned.py`). No effect on the example data. Decision recorded: keep T = T′ = 25, because false masks hurt partial genomes and shorten inferred branch lengths (conclusion 3). |
| 2026-10-23 | Masked TSV reports the source of each masked copy (section 2.8): best match by locus identity, `source_coherent` (outlier neighbours ≥ 90% AAI to each other) and, with the new `--taxonomy` option, the lowest common taxon of the matching genomes and of the masked genome with them. Neighbours by locus identity are now always ranked (only the best match without the high-identity test). Masks unchanged (section 4.14). Default parameters moved to section 2.9. |
| 2026-10-22 | Section 4.10: distribution of real loci over the Experiment 1 grid (`real_g_distribution.py`), to put per-bucket recall in context. |
| 2026-10-22 | Marker gene review (section 4.15, `gene_quality.py`): masks plus pairwise discordance over the whole AAI range, problem lineages; recommendation to drop three genes and review three. Conclusion 7. |
| 2026-10-22 | Section 4.15: order of the two tests and repeated masking. Order changes ~1% of masks (borderline); repeating converges after 4 rounds and concentrates in the flagged genes. Default unchanged. Note: the current binary masks 1,693 genome loci on the default run (1,695 recorded in section 4.1). |
| 2026-10-23 | Genome summary reports the main source of each genome's high-identity masks (`main_source_id`, `main_source_masks`; sections 2.8, 4.14). Three suspect genomes from neighbouring accessions share one main source. |
| 2026-10-22 | Section 4.16: what-if of non-decreasing trend lines (isotonic fit of the running medians); about 3% of masks change, mostly in the flagged genes. Not implemented. |
| 2026-10-23 | New default `--method monotone-median`: running medians made non-decreasing by isotonic regression, for both tests (sections 2.4, 4.16). Full data: 1,693 → 1,691 masks (23 added, 25 removed). Renumbered the marker gene review to section 4.15 (conclusion 8) and the trend section to 4.16 to avoid clashes with the review follow-ups. |
