# pf-variant-priority

**Designing and testing antimalarial resistance surveillance panels for *Plasmodium
falciparum*, on Indian parasites.**

Two questions, two answers, both measured on open data and both reported here whichever
way they came out.

**1. Does a resistance-marker panel chosen on Indian parasites beat one chosen globally?**
Yes, and the margin survives resampling. → [Result 1](#result-1-an-india-selected-panel-beats-a-globally-selected-one)

**2. Can an unsupervised protein language model rank variants of unknown significance?**
No. Scored against WHO-validated markers it is indistinguishable from chance, and that
negative is kept on the front page rather than buried. → [Result 2](#result-2-constraint-scoring-does-not-recover-resistance-markers)

The repository also holds the workflow those questions were asked with: a reproducible
Nextflow + Python pipeline that takes targeted Nanopore amplicon data (or a pre-computed
variant table), annotates protein-level variants across a 15-gene antimalarial resistance
panel, scores missense variants, and calibrates any score against known markers before
anyone acts on the ranking.

---

## Why this exists

Targeted amplicon surveillance of antimalarial resistance is now cheap enough to run at
national scale. Khilari *et al.* (*Open Forum Infectious Diseases*, 2026;13(3):ofag106)
sequenced 15 resistance-associated genes from dried blood spots across six Indian states
using Oxford Nanopore and the PfMDR15 panel, and found, among other things, **Pfaat1 S258L** —
a mutation previously reported only from Africa and associated with reduced lumefantrine
susceptibility.

That finding illustrates the bottleneck. Panels of this kind return many non-synonymous
variants per gene. A handful are WHO-validated. The rest are *variants of unknown
significance* (VUS), and today they are triaged by hand against the literature. As panels
widen and sample numbers grow, manual triage does not scale, and genuinely novel resistance
alleles can sit unremarked in a results table for years.

**This workflow adds a prioritisation layer, and a calibration step that decides whether
that layer is trustworthy.** It does not claim to predict resistance phenotype. It
produces a reproducible ordering over VUS *and* measures, against known markers, whether
the ordering carries any signal before anyone acts on it.

Two strands of work came out of that. The first asks a narrower, answerable question —
which markers belong in a panel aimed at Indian parasites — and gets a clear positive
answer. The second asks whether a protein language model can rank the unknowns, and gets a
clear negative. Both are below. A prioritisation tool that cannot say when it is not
working is worse than no tool at all, so the negative stays on the front page.

## Result 1: an India-selected panel beats a globally-selected one

**Question.** Surveillance panels are assembled from markers validated on global data. If
the markers that carry information in India are not the same ones, a panel inherited
wholesale is spending reads on the wrong positions. Is it?

**Data.** MalariaGEN Pf7, 300 Indian samples, scored against the project's drug-resistance
marker genotypes and inferred-resistance classifications.

**Answer.** At a panel size of seven markers, scored on Indian samples across 13 resistance
classifications:

| panel selected on | balanced accuracy |
|---|---|
| **Indian samples** | **0.9338** |
| global samples | 0.8540 |

Resampling the Indian samples 200 times gives a mean advantage of **+0.0410, 95% CI +0.0225
to +0.0609**, with **100% of resamples favouring the India-selected panel**. Quote that
bootstrap mean rather than the single-run gap of +0.0798, which is optimistic.

**Why.** The informativeness of individual markers differs, and the differences are
legible:

| marker | entropy, India | entropy, global |
|---|---|---|
| dhfr_108 | 1.298 | 0.635 |
| dhfr_59 | 1.306 | 0.822 |
| dhps_613 | 0.000 | 0.414 |

dhfr_108 carries about twice the entropy in India. dhps_613 is monomorphic across all 299
Indian calls, so on these samples it separates nothing.

**What this is not.** This is panel design — feature selection over MalariaGEN's
rule-based resistance classifications. It is not phenotype prediction, because those
classifications are rules applied to genotypes, not measured IC50 or survival assays.
Monomorphic also does not mean droppable: surveillance exists to catch what arrives next,
and a marker invariant today is the one you want to be watching when it changes.

Full method, code and caveats: [`india-panel/`](india-panel/).

## Approach

```
FASTQ ──minimap2──> BAM ──Clair3──> VCF ──bcftools norm/csq──> protein variants
                                                                      │
                                       ┌──────────────────────────────┤
                                       ▼                              ▼
                            ESM-2 masked-marginal              population metrics
                            substitution score                 (frequency, spread,
                                       │                        sample count)
                                       └──────────────┬───────────────┘
                                                      ▼
                                      calibration on WHO-validated markers
                                                      ▼
                                             ranked VUS watchlist
```

**Variant effect score.** Each missense variant is scored with the masked-marginal
log-odds from ESM-2 (Lin *et al.*, 2023):

```
score(wt → mut, position i) = log P(mut | x_masked_i) − log P(wt | x_masked_i)
```

A strongly negative score means the language model considers the substitution unlikely
given the protein family context — a proxy for functional constraint. This is an
unsupervised score: no resistance labels are used to produce it, which is what makes the
benchmark below meaningful.

**Calibration, not training.** The 13 WHO-validated *kelch13* markers are held out as
positive controls and a matched set of synonymous-adjacent / observed-neutral positions as
negatives. The workflow reports ROC-AUC, precision at k, and the rank of each validated
marker. If the score cannot recover known markers, the ranking over VUS is not trustworthy,
and the report says so rather than hiding it.

**Prioritisation.** The final rank combines the (calibrated) constraint score with
population-genetic evidence — allele frequency, number of states/sites where the variant is
seen, and whether frequency is rising across collection years where dates are available.
Weights live in `conf/priority.yml` and are deliberately explicit rather than learned, since
no labelled Indian VUS set exists to learn them from.

---

## Quick start

Two entry points. Use the second one if you do not have raw reads.

### A. From Nanopore FASTQ

```bash
nextflow run main.nf \
  --reads      'data/fastq/*.fastq.gz' \
  --reference  data/ref/Pfalciparum3D7_Genome.fasta \
  --gff        data/ref/Pfalciparum3D7.gff \
  --panel      data/pfmdr15_panel.tsv \
  --outdir     results/ \
  -profile     conda
```

### B. From an existing variant table (no sequencing data needed)

```bash
conda env create -f environment.yml && conda activate pf-vp

python bin/score_variants.py \
  --variants data/example_variant_table.tsv \
  --fasta    data/panel_proteins.faa \
  --out      results/scored.tsv

python bin/benchmark.py \
  --scored    results/scored.tsv \
  --validated data/who_validated_k13.tsv \
  --out       results/benchmark/

python bin/prioritise.py \
  --scored results/scored.tsv \
  --config conf/priority.yml \
  --out    results/watchlist.tsv
```

Entry point B runs on a laptop in minutes with the 35M-parameter ESM-2 checkpoint, or on a
single GPU in under an hour with the 650M checkpoint.

### Where to get real input data

| Source | What it gives you | Access |
|---|---|---|
| MalariaGEN **Pf7 / Pf8** | Open genotype calls for thousands of *P. falciparum* samples, including Indian sites | `malariagen_data` Python package |
| ENA **PRJNA956048** | Nanopore reads, targeted resistance panel (de Cesare *et al.*, *Nat Commun* 2024) | ENA / SRA |
| **PlasmoDB** release ≥ 52 | 3D7 reference genome, GFF, protein FASTA | plasmodb.org |
| WHO malaria resistance report | Current validated/candidate *kelch13* marker lists | who.int |

> `data/example_variant_table.tsv` in this repo is a **small synthetic file for smoke-testing
> the code path only**. It is clearly labelled as such and must not be used to draw any
> biological conclusion.

---

## Result 2: constraint scoring does not recover resistance markers

Scored every possible missense substitution across the kelch13 propeller domain
(positions 350–726) with ESM-2 650M on a T4 GPU. Sequence retrieved from UniProt
(`A0A077LQB4_PLAFA`, Kelch protein, *P. falciparum*); all 13 WHO-validated markers matched
the reference residue, confirming the right isoform.

```
substitutions scored : 7163
positives (validated): 13
negatives            : 7137   (WHO candidates excluded)
ROC-AUC              : 0.524   95% CI 0.383 – 0.664
median rank of validated markers: 2662 of 7163
```

| Marker | ESM score | Rank | Percentile |
|---|---|---|---|
| R622I | −8.913 | 683 / 7163 | top 9.5% |
| C580Y | −6.354 | 1542 | top 21.5% |
| I543T | −5.815 | 1749 | top 24.4% |
| N458Y | −5.451 | 1905 | top 26.6% |
| Y493H | −4.507 | 2412 | top 33.7% |
| P574L | −4.491 | 2428 | top 33.9% |
| P553L | −4.131 | 2662 | top 37.2% |
| M476I | −2.756 | 3690 | top 51.5% |
| F446I | −1.768 | 4584 | top 64.0% |
| R561H | −1.158 | 5151 | top 71.9% |
| C469Y | −0.767 | 5509 | top 76.9% |
| R539T | −0.382 | 5844 | top 81.6% |
| A675V | −0.062 | 6131 | top 85.6% |

### What this means

**The constraint score does not recover WHO-validated artemisinin resistance markers.**
AUC 0.524 with a confidence interval spanning 0.5 is indistinguishable from chance. Per the
calibration rule stated above, the prioritised watchlist built on this score alone is
**not usable for triage**, and this README reports that rather than burying it.

This is not a bug, and on reflection it is the biologically expected answer. It is the
caveat in the limitations section, confirmed by measurement:

> A highly constrained position may be essential for parasite viability and therefore
> *less* likely to be a resistance allele.

Resistance mutations are, by definition, substitutions the parasite **tolerates** — they
arise, survive and spread under drug pressure. They sit under *positive* selection, not
purifying selection. An unsupervised model of evolutionary constraint measures the
opposite quantity. The markers scatter across the whole distribution (R622I in the top
10%, A675V in the bottom 15%), which is what that mismatch looks like.

### What would be needed instead

1. **Signatures of positive selection** from population data — iHS, XP-EHH, Tajima's D at
   the variant's locus in MalariaGEN Pf7/Pf8. Resistance alleles should show recent
   selective sweeps; that is a signal pointing the right way.
2. **Structural context** — AlphaFold3 models of the panel proteins, with proximity to
   drug-binding or substrate pockets as a feature. For kelch13, position within the
   propeller and the BTB/POZ interface.
3. **Supervised learning on phenotype**, once *in vitro* IC50 or ring-stage survival assay
   data are linked. The constraint score then becomes one weak feature among several
   rather than the whole model.
4. **A tolerance score rather than a constraint score** — asking "is this substitution
   compatible with viability?" instead of "is it unusual?" is closer to the right question.

The pipeline, the calibration harness and the reporting are unchanged by this; only the
input score needs replacing. That is the point of keeping calibration as a separate,
mandatory step.

---

## Repository layout

```
pf-variant-priority/
├── main.nf                       Nextflow workflow (reads → annotated variants)
├── nextflow.config               executor, container and conda profiles
├── conf/priority.yml             prioritisation weights, documented and editable
├── bin/
│   ├── collect_variants.py       parse bcftools csq output → tidy variant table
│   ├── score_variants.py         ESM-2 masked-marginal scoring
│   ├── benchmark.py              calibration against WHO-validated markers
│   └── prioritise.py             combine evidence → ranked watchlist
├── data/
│   ├── pfmdr15_panel.tsv         panel gene list (VERIFY before use — see below)
│   ├── who_validated_k13.tsv     WHO validated + candidate kelch13 markers
│   └── example_variant_table.tsv synthetic smoke-test input
├── india-panel/                  Result 1: India vs global panel selection on Pf7
│   ├── india_panel.py            entropy, greedy selection, bootstrap
│   ├── test_synthetic.py         selection recovers a planted answer
│   ├── india_panel_colab.ipynb   the run, as executed
│   └── README.md                 method, full results and caveats
└── notebooks/01_walkthrough.ipynb
```

---

## Honest limitations

These are stated up front because they determine how the output should be read.

1. **Low-complexity regions.** *P. falciparum* proteins carry long asparagine and lysine
   repeats arising from an AT-rich genome. Protein language models are trained mostly on
   non-repetitive sequence and are poorly calibrated inside such regions. The workflow
   flags variants falling in low-complexity stretches (`--mask-lcr`, on by default) and
   reports them separately rather than ranking them alongside the rest.
2. **A constraint score is not a resistance prediction.** ESM measures evolutionary
   plausibility, not drug phenotype. A highly constrained position may be essential for
   parasite viability and therefore *less* likely to be a resistance allele. The output is
   a triage ordering for follow-up, nothing more.
3. **Small positive set.** Thirteen validated *kelch13* markers is a thin benchmark, and
   they all sit in one gene and one structural domain. AUC from such a set is indicative,
   not definitive. Confidence intervals are bootstrapped and reported.
4. **Panel gene list must be verified.** `data/pfmdr15_panel.tsv` lists the canonical
   antimalarial resistance genes with PlasmoDB IDs, but the exact 15-gene composition of
   the published PfMDR15 panel should be confirmed against the Supplementary Table of
   Khilari *et al.* (2026) before any comparison with that study is claimed.
5. **No phenotype data.** Linking these rankings to *in vitro* IC50 or clinical treatment
   outcome is the obvious next step and is not attempted here.

## Possible extensions

- Swap ESM-2 for a *Plasmodium*-adapted model fine-tuned on apicomplexan proteomes.
- Add structural context: AlphaFold3 models of panel proteins, with relative solvent
  accessibility and proximity to known binding pockets as additional features.
- Train a supervised genotype→phenotype model once labelled IC50 data are linked, using the
  unsupervised score as one feature.
- Extend the same framework to insecticide-resistance loci in *Anopheles stephensi*.

## Citation / provenance

Built as an independent methods exercise. Method and motivation follow:

- Khilari A, Sharma S, Bajpai M, Viswan AK, Chaturvedi R, Mirdha BR, Rahi M, Sharma A,
  Shanmugam D. Targeted genomic surveillance unveils genetic variations linked to regional
  malaria drug resistance dynamics in India. *Open Forum Infect Dis.* 2026;13(3):ofag106.
- de Cesare M, *et al.* Flexible and cost-effective genomic surveillance of *P. falciparum*
  malaria with targeted nanopore sequencing. *Nat Commun.* 2024;15:1413.
- Lin Z, *et al.* Evolutionary-scale prediction of atomic-level protein structure with a
  language model. *Science.* 2023;379:1123–1130.

Author: Anjali Chaudhary

## Copyright

Copyright (c) 2026 Anjali Chaudhary. All rights reserved. Contact the repository owner for permission to reuse the original project material.
