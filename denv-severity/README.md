# DENV genome-to-severity, with the controls

Companion to the main `pf-variant-priority` workflow, aimed at the dengue/Zika side of
the lab's work.

## The question

Ravi *et al.* (*PLoS Negl Trop Dis*, 2025) sequenced 3,254 DENV-2 genomes from the 2023
Delhi outbreak on Oxford Nanopore and linked SNV hotspots in **E, NS4B and NS5** to
disease severity using Fisher's exact test and phi-correlation.

Site-by-site association testing on viral genomes has two standing problems:

- **Phylogenetic structure.** Genomes are related by descent, so samples are not
  independent. If one lineage happens to circulate in a sicker population, every SNV that
  marks that lineage looks associated with severity.
- **Multiple testing across linked sites.** Thousands of correlated positions are tested
  one at a time, and linkage means a single causal variant drags its neighbours along.

## What this script adds

| Control | Why |
|---|---|
| Multivariate L1-logistic + gradient boosting | Correlated positions compete instead of each claiming the signal |
| Nested cross-validation | Reported AUC is not inflated by hyperparameter tuning on the test set |
| Label-permutation null | Gives the AUC distribution expected by chance for *this* dataset |
| Lineage-aware grouped CV (`--clade-col`) | A clade never spans train and test, so the model cannot learn lineage instead of biology |
| Stability selection | Positions ranked by how often they survive across resamples, not one lucky fit |

The script refuses to report associated positions when the observed AUC does not beat its
own permutation null.

## Validation

Tested on synthetic alignments with signal planted at known alignment columns:

- **Signal present** (3 causal sites, 120 genomes, 43 SNV features): nested-CV AUC 0.962,
  permutation *p* = 0.048 at 20 permutations, and all three planted sites recovered at
  stability-selection frequency 1.0.
- **No signal** (monomorphic "causal" columns, by accident during development): AUC 0.567,
  *p* = 0.29, correctly reported as no signal. This was a genuine negative control — the
  test data were broken and the script said so rather than producing a hit list.

Reproduce both with the generator in `make_synthetic.py`.

## Usage

```bash
python denv_severity_ml.py \
  --alignment data/denv2_aligned.fasta \
  --metadata  data/metadata.tsv \
  --id-col sample --label-col severity --positive-label severe \
  --clade-col lineage \
  --n-permutations 100 \
  --out results/
```

Outputs `model_report.json` (AUC, fold AUCs, null distribution, permutation *p*, verdict)
and `stability_selection.tsv` (positions by selection frequency).

## Getting data

Public DENV-2 genomes are in GenBank and NextStrain. **Clinical severity labels for the
Delhi cohort are not public** — supply labels from your own cohort, or run the machinery
against a non-clinical surrogate label (collection year, region) to demonstrate the
method. Any such run must be reported as a methods demonstration, not an epidemiological
finding.

## Limitations

- Binary severity collapses a graded clinical picture; an ordinal model would be better.
- Grouped CV controls lineage crudely. A proper phylogenetic mixed model (e.g. treeWAS,
  or a kinship-corrected LMM) is the rigorous alternative.
- Host genetics, prior serotype exposure and age drive dengue severity strongly and are
  not modelled here. Viral genotype is at best a partial explanation.
