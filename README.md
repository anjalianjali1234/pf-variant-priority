# pf-variant-priority

**Prioritising variants of unknown significance in *Plasmodium falciparum* drug-resistance
surveillance panels, using protein language models.**

A reproducible Nextflow + Python workflow that takes targeted Nanopore amplicon data (or a
pre-computed variant table), annotates protein-level variants across a 15-gene antimalarial
resistance panel, scores every missense variant with a protein language model, calibrates
those scores against WHO-validated *kelch13* markers, and emits a ranked watchlist of
unvalidated variants circulating in a given set of samples.

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

**This workflow adds a prioritisation layer.** It does not claim to predict resistance
phenotype. It produces a defensible, reproducible ordering over VUS, so that limited
follow-up effort (*in vitro* assays, phenotyping, targeted resequencing) goes to the
variants most likely to matter.

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

Author: Anjali Chaudhary · MIT licence
