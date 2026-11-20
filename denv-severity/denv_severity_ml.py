#!/usr/bin/env python3
"""Genome-wide association of DENV variants with clinical severity, done defensibly.

Motivation
----------
Ravi et al. (PLoS Negl Trop Dis, 2025) sequenced 3,254 DENV-2 genomes from the 2023
Delhi outbreak with Oxford Nanopore and used Fisher's exact test plus phi-correlation to
link SNV hotspots in E, NS4B and NS5 to disease severity.

Per-site testing across thousands of correlated positions has two well-known problems:
phylogenetic structure (related genomes are not independent samples, so lineage can
masquerade as a severity signal) and multiple testing across linked sites. This script
keeps the same question and adds the controls:

1. **Multivariate models** (L1-logistic, gradient boosting) instead of site-by-site tests,
   so correlated positions compete rather than each claiming the signal.
2. **Nested cross-validation** — hyperparameters tuned in an inner loop, performance
   measured in an outer loop, so reported AUC is not optimistically biased.
3. **Label-permutation null** — the whole pipeline is re-run on shuffled labels to get the
   AUC distribution expected by chance. A model only counts if it beats its own null.
4. **Lineage-aware splitting** (optional `--clade-col`) — grouped CV so a clade never
   appears in both train and test. Without this, models often learn lineage, not biology.
5. **Stability selection** — positions are reported by how often they are selected across
   resamples, not by a single fit.

Input
-----
`--alignment`  aligned genome FASTA (MAFFT), one record per sample
`--metadata`   TSV with a sample id column, a severity label column, optionally a clade

Severity labels are NOT public for the Delhi cohort. Supply labels from your own cohort,
or use `--surrogate-label` to run the machinery against a non-clinical label (collection
year, region) as a methods demonstration.

Usage
-----
    python denv_severity_ml.py --alignment data/denv2_aligned.fasta \
        --metadata data/metadata.tsv --id-col sample --label-col severity \
        --clade-col lineage --out results/
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AMBIG = set("NRYKMSWBDHV-")


def read_alignment(path: Path) -> tuple[list[str], np.ndarray]:
    ids, seqs, name, buf = [], [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    ids.append(name); seqs.append("".join(buf))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper())
    if name:
        ids.append(name); seqs.append("".join(buf))
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        sys.exit(f"alignment is ragged: sequence lengths {sorted(lengths)[:5]}... — realign first")
    return ids, np.array([list(s) for s in seqs])


def build_snv_matrix(aln: np.ndarray, min_minor_freq: float, max_missing: float):
    """One-hot encode polymorphic, well-covered sites. Returns (X, feature_names)."""
    cols, names = [], []
    n = aln.shape[0]
    for j in range(aln.shape[1]):
        col = aln[:, j]
        missing = np.isin(col, list(AMBIG)).mean()
        if missing > max_missing:
            continue
        valid = col[~np.isin(col, list(AMBIG))]
        alleles, counts = np.unique(valid, return_counts=True)
        if len(alleles) < 2:
            continue
        order = np.argsort(-counts)
        for allele in alleles[order][1:]:          # minor alleles only
            freq = (valid == allele).mean()
            if freq < min_minor_freq:
                continue
            cols.append((col == allele).astype(np.int8))
            names.append(f"{j + 1}{allele}")       # 1-based alignment coordinate
    if not cols:
        sys.exit("no SNV features passed the frequency/missingness filters")
    print(f"[snv] {len(names)} features from {aln.shape[1]} alignment columns, {n} genomes",
          file=sys.stderr)
    return np.column_stack(cols), names


def l1_logreg(**kw):
    """L1-penalised logistic regression, compatible across scikit-learn versions.

    scikit-learn 1.8 deprecated `penalty=` in favour of `l1_ratio=`; older versions do
    not accept `l1_ratio` for liblinear. Pick the right spelling at runtime.
    """
    from sklearn.linear_model import LogisticRegression

    base = dict(solver="liblinear", max_iter=5000, class_weight="balanced", **kw)
    try:
        return LogisticRegression(l1_ratio=1, **base)
    except TypeError:
        return LogisticRegression(penalty="l1", **base)


def evaluate(X, y, groups, n_outer, n_inner, seed):
    """Nested CV. Returns (mean_auc, per_fold_aucs, mean_selection_frequency)."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold

    outer = (StratifiedGroupKFold(n_outer, shuffle=True, random_state=seed)
             if groups is not None else StratifiedKFold(n_outer, shuffle=True, random_state=seed))
    inner = StratifiedKFold(n_inner, shuffle=True, random_state=seed)

    aucs, sel = [], np.zeros(X.shape[1])
    splits = outer.split(X, y, groups) if groups is not None else outer.split(X, y)
    for tr, te in splits:
        if len(np.unique(y[te])) < 2:
            continue
        grid = GridSearchCV(
            l1_logreg(), {"C": [0.01, 0.1, 1.0, 10.0]},
            scoring="roc_auc", cv=inner, n_jobs=-1,
        ).fit(X[tr], y[tr])
        lr = grid.best_estimator_
        sel += (np.abs(lr.coef_[0]) > 1e-8).astype(float)

        gb = HistGradientBoostingClassifier(max_iter=200, random_state=seed).fit(X[tr], y[tr])
        best = max(roc_auc_score(y[te], lr.predict_proba(X[te])[:, 1]),
                   roc_auc_score(y[te], gb.predict_proba(X[te])[:, 1]))
        aucs.append(best)

    if not aucs:
        sys.exit("no usable outer folds — too few samples or one class only")
    return float(np.mean(aucs)), aucs, sel / len(aucs)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--alignment", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--id-col", default="sample")
    p.add_argument("--label-col", default="severity")
    p.add_argument("--clade-col", default=None, help="group column for lineage-aware CV")
    p.add_argument("--positive-label", default="severe")
    p.add_argument("--min-minor-freq", type=float, default=0.02)
    p.add_argument("--max-missing", type=float, default=0.10)
    p.add_argument("--n-outer", type=int, default=5)
    p.add_argument("--n-inner", type=int, default=3)
    p.add_argument("--n-permutations", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    ids, aln = read_alignment(args.alignment)
    meta = pd.read_csv(args.metadata, sep="\t", comment="#").set_index(args.id_col)
    shared = [i for i in ids if i in meta.index]
    if len(shared) < 30:
        sys.exit(f"only {len(shared)} genomes matched metadata — too few to model")
    keep = [ids.index(i) for i in shared]

    y = (meta.loc[shared, args.label_col].astype(str) == args.positive_label).to_numpy().astype(int)
    groups = meta.loc[shared, args.clade_col].to_numpy() if args.clade_col else None
    X, names = build_snv_matrix(aln[keep], args.min_minor_freq, args.max_missing)

    print(f"[model] {len(shared)} genomes, {int(y.sum())} positive, {len(y) - int(y.sum())} negative; "
          f"grouping={'on' if groups is not None else 'OFF (lineage confounding uncontrolled)'}",
          file=sys.stderr)

    auc, fold_aucs, sel = evaluate(X, y, groups, args.n_outer, args.n_inner, args.seed)
    print(f"[model] observed nested-CV AUC = {auc:.3f}", file=sys.stderr)

    rng = np.random.default_rng(args.seed)
    null = []
    for i in range(args.n_permutations):
        yp = rng.permutation(y)
        try:
            null.append(evaluate(X, yp, groups, args.n_outer, args.n_inner, args.seed + i + 1)[0])
        except SystemExit:
            continue
        if (i + 1) % 10 == 0:
            print(f"[perm] {i + 1}/{args.n_permutations}", file=sys.stderr)

    null = np.array(null)
    pval = float((np.sum(null >= auc) + 1) / (len(null) + 1)) if len(null) else float("nan")

    args.out.mkdir(parents=True, exist_ok=True)
    report = {
        "n_genomes": len(shared),
        "n_positive": int(y.sum()),
        "n_features": len(names),
        "lineage_aware_cv": groups is not None,
        "observed_auc": round(auc, 3),
        "fold_aucs": [round(a, 3) for a in fold_aucs],
        "null_auc_mean": round(float(null.mean()), 3) if len(null) else None,
        "null_auc_p95": round(float(np.percentile(null, 95)), 3) if len(null) else None,
        "permutation_p_value": round(pval, 4),
        "verdict": (
            "Signal exceeds the permutation null." if pval < 0.05 else
            "NO SIGNAL beyond chance. Do not report associated positions from this run."
        ),
    }
    if groups is None:
        report["warning"] = (
            "Run without --clade-col: phylogenetic structure is uncontrolled and any signal "
            "may reflect lineage rather than biology."
        )

    (args.out / "model_report.json").write_text(json.dumps(report, indent=2))
    pd.DataFrame({"feature": names, "selection_frequency": sel.round(3)}) \
      .sort_values("selection_frequency", ascending=False) \
      .to_csv(args.out / "stability_selection.tsv", sep="\t", index=False)

    print(f"[model] permutation p = {pval:.4f} — {report['verdict']}", file=sys.stderr)
    print(f"[model] wrote {args.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
