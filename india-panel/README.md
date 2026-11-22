# An India-optimised antimalarial resistance marker panel

## The problem

Resistance marker panels in common use were designed around the Southeast Asian and
African resistance landscape. India's is different in ways that matter operationally:

- **Artesunate + sulfadoxine-pyrimethamine is first-line** for uncomplicated *P. falciparum*
  malaria in most of India, except six northeastern states where treatment failure rates
  were high.
- ***kelch13* is largely wild type.** A Madhya Pradesh study found stable high frequencies
  of SP-resistance mutations and an *absence* of K13 mutations three and four years after
  AS+SP was introduced.
- **Artemisinin resistance in India is not always *kelch13*-mediated.** A central India
  study reported day-3 positive parasites showing the artemisinin-resistance phenotype
  with no pfkelch13 mutation.
- **SP resistance is the live operational problem.** AS-SP failure rates above 10% and
  high prevalence of resistant SP parasites, including quintuple mutants, are reported in
  northeastern India, with zones of predicted PfDHPS 540E prevalence above 15% in central,
  eastern and northeastern India.

A panel inherited from elsewhere therefore spends sequencing capacity on positions that
are invariant in Indian parasites, while the markers that drive Indian treatment policy
compete for the same read budget.

## The question

Which markers carry information in Indian parasites? What is the smallest panel that
recovers the resistance picture there? And does a panel selected on Indian data differ
from one selected on global data?

## What this is, and is not

MalariaGEN's resistance classifications are **heuristic rules applied to genotypes**, not
measured phenotypes — a parasite is called pyrimethamine resistant if dhfr 108N is
present, irrespective of other alleles. So this is a **panel design and feature selection**
exercise: how few markers reproduce the agreed interpretation in Indian parasites.

It is **not** independent phenotype prediction. Training a model to predict rule-derived
labels from the same genotypes that generated them and calling it phenotype prediction
would be circular, and a reviewer would say so immediately.

## Method

1. **Informativeness.** Per marker, allele frequency and Shannon entropy in Indian samples
   versus globally. The `uninformative_in_india_only` flag marks markers that vary
   worldwide but are invariant in India.
2. **Panel selection.** Greedy forward selection, adding the marker that most improves
   mean balanced accuracy across the resistance classifications, with cross-validation at
   each step. Run twice: selected on Indian samples, and selected on all samples.
3. **Comparison.** Both panels scored on Indian samples, giving an accuracy-versus-panel-size
   curve. A gap at small panel sizes is the finding.
4. **Regional breakdown.** Marker frequencies by country across South Asia.

---

## Results — MalariaGEN Pf7, run 23 Nov 2026

**Scope.** 16,203 QC-pass samples with metadata, marker genotypes and resistance
classifications. **India: n = 300.** 15 marker columns, 13 resistance labels. India is
modestly sampled in Pf7 (Bangladesh 1,310; Myanmar 985), and every claim below is bounded
by that 300.

### The India-selected panel beats the globally-selected one on Indian parasites

Mean balanced accuracy across the 13 resistance classifications, both panels scored on
Indian samples:

| markers | selected on India | selected globally | difference |
|---|---|---|---|
| 1 | 0.6310 | 0.6310 | 0.0000 |
| 2 | 0.7465 | 0.7285 | +0.0180 |
| 3 | 0.7921 | 0.7725 | +0.0196 |
| 4 | 0.8369 | 0.8040 | +0.0329 |
| 5 | 0.8684 | 0.8523 | +0.0161 |
| 6 | 0.9007 | 0.8383 | +0.0624 |
| **7** | **0.9338** | **0.8540** | **+0.0798** |
| 8 | 0.9370 | 0.9180 | +0.0190 |
| 9 | 0.9375 | 0.9304 | +0.0071 |
| 10 | 0.9378 | 0.9284 | +0.0094 |

**Bootstrap at k = 7** (200 resamples of the Indian samples): mean difference **+0.0410,
95% CI +0.0225 to +0.0609, and 100% of resamples favour the India-selected panel.**

The interval excludes zero, so the gap is not noise. Note the bootstrap mean (+0.041) is
about half the single-run point estimate (+0.080) — the point estimate is optimistic, and
the bootstrap is the number to quote.

### The two panels

| rank | selected on India | selected globally |
|---|---|---|
| 1 | pm2_dup_call | pm2_dup_call |
| 2 | **dhps_437[G]** | dhfr_51[N] |
| 3 | mdr1_dup_call | crt_76[K] |
| 4 | **dhfr_108[S]** | mdr1_dup_call |
| 5 | crt_76[K] | dhps_437[G] |
| 6 | kelch13_349-726_ns_changes | **dhps_540[K]** |
| 7 | dhfr_51[N] | **dhfr_59[C]** |

Same universe of markers, different ordering. India pulls the SP markers *dhps_437* and
*dhfr_108* forward; the global ordering spends slots 6 and 7 on *dhps_540* and *dhfr_59*,
which buys less in Indian parasites.

### Which markers carry information in India

Shannon entropy over called genotypes, India versus global:

| marker | India | global | reading |
|---|---|---|---|
| dhfr_108[S] | 1.298 | 0.635 | **twice as informative in India** |
| dhfr_59[C] | 1.306 | 0.822 | more informative in India |
| crt_76[K] | 1.400 | 1.299 | slightly more |
| dhps_437[G] | 1.126 | 1.029 | slightly more |
| dhps_540[K] | 1.010 | 1.598 | **less informative in India** |
| dhfr_51[N] | 0.825 | 1.191 | less informative in India |
| dhfr_164[I] | 0.123 | 0.891 | **near-invariant in India** (98.7% wild type) |
| dhps_613[A] | 0.000 | 0.414 | **monomorphic in India** (299/299 calls `A`) |

The *dhfr_108* and *dhfr_59* contrast is the clearest signal and matches the known Indian
SP-resistance picture. *dhfr_164* and *dhps_613* carry essentially nothing for
distinguishing Indian parasites.

### A correction to the framing

Before running this, the project argued that *kelch13* would be uninformative in India.
**The data does not support that as stated.** The `kelch13_349-726_ns_changes` column has
*higher* entropy in Indian samples (3.459) than globally (3.160), and the India-selected
panel includes it at rank 6.

That column lists every non-synonymous change found in the propeller domain, so its
entropy reflects the diversity of all such strings, not the frequency of WHO-validated
resistance markers. The literature point — that validated *kelch13* markers are rare in
India — concerns a different quantity and is not tested here. Either measure it directly
or drop the claim. It is not evidence either way as it stands.

### A bug this run exposed

`pm2_breakpoint` was initially flagged as monomorphic in India. It has **exactly one**
Indian call. "Invariant" there is an artefact of missingness, not a fact about Indian
parasites. The code now requires at least 50 calls before calling a marker monomorphic
and reports a `too_few_calls` flag separately. *dhps_613[A]*, with 299 calls, survives
that check and is a real finding.

---

## Validation

`test_synthetic.py` builds data with a known answer: 20 markers where Indian resistance
status is driven by dhfr/dhps, and global status additionally by *kelch13* and two *crt*
positions that are monomorphic in the Indian subset. All three tests pass:

```
TEST 1  flags exactly the 6 markers that are invariant in India but variable globally
TEST 2  India panel  : dhfr_108, crt_76, dhps_437, dhps_540, dhfr_51          (5 markers)
        Global panel : dhfr_108, kelch13_580, crt_76, dhps_437, dhps_540, dhfr_51  (6)
        India panel excludes kelch13; global panel spends a slot on it
TEST 3  on Indian samples, India-selected beats globally-selected
        k=2  0.8366 vs 0.7075   (+0.129)
        k=3  0.9371 vs 0.8366   (+0.101)
        k=4  0.9871 vs 0.9371   (+0.050)
```

Run it with `python test_synthetic.py`. These are synthetic numbers that validate the
machinery, not results about real parasites.

## Running on real data

Open `india_panel_colab.ipynb` in Google Colab and run top to bottom. No GPU needed.
It installs `malariagen_data`, shows the South Asian sample counts and the resistance
table columns before analysing anything, then runs the panel selection.

```bash
python india_panel.py --release Pf7 --country India --max-markers 12 --out results_india
```

Outputs: `marker_informativeness.tsv`, `panel_comparison.tsv`,
`accuracy_vs_panel_size.png`, `regional_heterogeneity.tsv`, `summary.json`.

## Two things that must survive into any write-up

**A monomorphic marker is not a useless marker.** A surveillance panel both measures the
resistance already present and catches what arrives next. A *kelch13* marker invariant in
Indian samples today is uninformative for burden estimation and still worth keeping as a
sentinel — especially given reports of artemisinin-resistance phenotypes with wild-type
*kelch13* in India. The recommendation separates "uninformative for burden" from "safe to
remove", and never collapses the two.

**Sample size governs the claim.** Pf7 is not evenly sampled. If this rests on a few
hundred Indian genomes, that belongs in the first paragraph, not the limitations. The
script prints the count and refuses to pretend otherwise; if India is under ~150 samples
it tells you to widen the scope to South Asia and relabel the work.

## Sources

- MalariaGEN. Pf7: an open dataset of *Plasmodium falciparum* genome variation in 20,000
  worldwide samples. *Wellcome Open Res.* 2023;8:22.
- MalariaGEN. Mapping genetic markers to resistance status classification (Pf7), 2021.
- Khilari A, *et al.* Targeted genomic surveillance unveils genetic variations linked to
  regional malaria drug resistance dynamics in India. *Open Forum Infect Dis.* 2026;13(3):ofag106.
- Systematic review and geospatial modeling of molecular markers of resistance to
  artemisinins and sulfadoxine-pyrimethamine in *P. falciparum* in India. *Am J Trop Med Hyg.* 2024.
