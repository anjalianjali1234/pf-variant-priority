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
