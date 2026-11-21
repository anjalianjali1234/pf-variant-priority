#!/usr/bin/env python3
"""Design a drug-resistance marker panel optimised for Indian P. falciparum.

The problem
-----------
Resistance marker panels in common use were designed around the Southeast Asian and
African resistance landscape. India's is different: kelch13 is largely wild type, while
sulfadoxine-pyrimethamine resistance is widespread and artesunate+SP is first-line
treatment in most of the country. A panel inherited from elsewhere spends sequencing
capacity on markers that are invariant in Indian parasites.

The question
------------
Which markers actually carry information in Indian parasites, and what is the smallest
panel that recovers the resistance picture there? Does a panel selected on Indian data
differ from one selected on global data?

What this is, and is not
------------------------
MalariaGEN's resistance classifications are HEURISTIC RULES applied to genotypes, not
measured phenotypes (see Pf7_resistance_classification.pdf: a parasite is called
pyrimethamine resistant if dhfr 108N is present, irrespective of other alleles). So this
is a PANEL DESIGN and FEATURE SELECTION exercise: which markers are needed to reproduce
the agreed interpretation in Indian parasites. It is NOT independent phenotype
prediction, and reporting it as such would be circular.

Output
------
  marker_informativeness.tsv   per-marker allele frequency, entropy, monomorphic flag
  panel_comparison.tsv         India-selected vs globally-selected panel, ranked
  accuracy_vs_panel_size.png   saturation curve for both panels, scored on India
  regional_heterogeneity.tsv   marker frequencies by region within South Asia
  summary.json                 headline numbers
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Genes in the lab's targeted panel (Khilari et al., OFID 2026). Used for comparison only.
PANEL_15 = ["crt", "mdr1", "dhfr", "dhps", "kelch13", "aat1", "mdr2", "mrp1",
            "ubp1", "ap2mu", "coronin", "arps10", "fd", "exo", "plasmepsin2-3"]

SOUTH_ASIA = ["India", "Bangladesh", "Myanmar", "Nepal", "Sri Lanka", "Pakistan", "Bhutan"]


# ---------------------------------------------------------------- data loading

def load_release(prefer: str = "Pf8"):
    """Return (sample_metadata, release_name).

    Tries the preferred release then the other. Column names differ between releases,
    so everything downstream detects columns rather than hard-coding them.
    """
    import malariagen_data

    for name in ([prefer] + [r for r in ("Pf8", "Pf7") if r != prefer]):
        try:
            rel = getattr(malariagen_data, name)()
            meta = rel.sample_metadata()
            print(f"[data] loaded {name}: {meta.shape[0]} samples, {meta.shape[1]} columns")
            return meta, name
        except Exception as exc:
            print(f"[data] {name} unavailable: {exc}", file=sys.stderr)
    sys.exit("could not load any MalariaGEN release")


DR_URLS = {
    "Pf7": {
        "class": "https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_inferred_resistance_status_classification.txt",
        "geno": "https://www.malariagen.net/wp-content/uploads/2023/11/Pf7_drug_resistance_marker_genotypes.txt",
    },
}


def load_dr_tables(release: str, class_url: str | None, geno_url: str | None):
    urls = DR_URLS.get(release, {})
    cu = class_url or urls.get("class")
    gu = geno_url or urls.get("geno")
    if not cu or not gu:
        sys.exit(
            f"No built-in resistance-table URLs for {release}. Pass --class-url and "
            "--geno-url. Find them in the release README on malariagen.net."
        )
    cls = pd.read_csv(cu, sep="\t")
    gen = pd.read_csv(gu, sep="\t")
    print(f"[data] classifications {cls.shape}, marker genotypes {gen.shape}")
    return cls, gen


def find_id_column(df: pd.DataFrame) -> str:
    for c in ("Sample", "sample", "Sample ID", "sample_id"):
        if c in df.columns:
            return c
    return df.columns[0]


def find_qc_column(df: pd.DataFrame) -> str | None:
    for c in df.columns:
        if c.lower().replace("_", " ").strip() in ("qc pass", "qcpass"):
            return c
    return None


def find_country_column(df: pd.DataFrame) -> str:
    for c in ("Country", "country"):
        if c in df.columns:
            return c
    sys.exit(f"no country column found in metadata: {list(df.columns)[:20]}")


# ------------------------------------------------------- marker informativeness

def shannon_entropy(series: pd.Series) -> float:
    """Entropy in bits over non-missing genotype calls."""
    v = series.dropna()
    v = v[~v.astype(str).isin(["-", "", "nan", "NaN", "undetermined"])]
    if len(v) == 0:
        return float("nan")
    p = v.value_counts(normalize=True).to_numpy()
    return float(-(p * np.log2(p)).sum())


def marker_table(gen: pd.DataFrame, marker_cols: list[str], label: str) -> pd.DataFrame:
    rows = []
    for c in marker_cols:
        v = gen[c].dropna()
        v = v[~v.astype(str).isin(["-", "", "nan", "NaN", "undetermined"])]
        counts = v.value_counts(normalize=True)
        rows.append({
            "marker": c,
            f"{label}_n_called": int(len(v)),
            f"{label}_n_alleles": int(v.nunique()),
            f"{label}_major_allele": counts.index[0] if len(counts) else None,
            f"{label}_major_freq": round(float(counts.iloc[0]), 4) if len(counts) else None,
            f"{label}_entropy_bits": round(shannon_entropy(gen[c]), 4),
            f"{label}_monomorphic": bool(v.nunique() <= 1),
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------- panel selection

def encode(gen: pd.DataFrame, marker_cols: list[str]) -> pd.DataFrame:
    """Ordinal-encode each marker column. Missing becomes its own category."""
    out = {}
    for c in marker_cols:
        out[c] = gen[c].astype(str).fillna("NA").astype("category").cat.codes
    return pd.DataFrame(out, index=gen.index)


def score_subset(X: pd.DataFrame, Y: pd.DataFrame, cols: list[str], cv: int, seed: int) -> float:
    """Mean balanced accuracy across resistance labels, using a shallow tree."""
    from sklearn.model_selection import cross_val_score
    from sklearn.tree import DecisionTreeClassifier

    if not cols:
        return 0.0
    scores = []
    for lab in Y.columns:
        y = Y[lab].to_numpy()
        if len(np.unique(y)) < 2:
            continue  # label invariant here: carries no information to recover
        k = min(cv, int(pd.Series(y).value_counts().min()))
        if k < 2:
            continue
        clf = DecisionTreeClassifier(max_depth=min(len(cols) + 1, 8),
                                     class_weight="balanced", random_state=seed)
        s = cross_val_score(clf, X[cols], y, cv=k, scoring="balanced_accuracy")
        scores.append(float(np.mean(s)))
    return float(np.mean(scores)) if scores else 0.0


def greedy_panel(X: pd.DataFrame, Y: pd.DataFrame, max_k: int, cv: int, seed: int,
                 tag: str) -> list[tuple[str, float]]:
    """Forward selection: repeatedly add the marker that most improves mean accuracy."""
    chosen: list[tuple[str, float]] = []
    remaining = list(X.columns)
    current = 0.0
    for step in range(1, min(max_k, len(remaining)) + 1):
        best, best_score = None, current
        for c in remaining:
            s = score_subset(X, Y, [m for m, _ in chosen] + [c], cv, seed)
            if s > best_score + 1e-6:
                best, best_score = c, s
        if best is None:
            print(f"[{tag}] no further marker improves accuracy; stopping at {step - 1}")
            break
        chosen.append((best, round(best_score, 4)))
        remaining.remove(best)
        current = best_score
        print(f"[{tag}] {step:>2}. {best:<28} cumulative balanced accuracy {best_score:.4f}")
    return chosen


# --------------------------------------------------------------------- reporting

def plot_curve(curves: dict[str, list[float]], out: Path, n_india: int) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[plot] skipped: {exc}", file=sys.stderr)
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    for name, ys in curves.items():
        ax.plot(range(1, len(ys) + 1), ys, marker="o", lw=2, label=name)
    ax.set_xlabel("markers in panel")
    ax.set_ylabel("mean balanced accuracy (scored on Indian samples)")
    ax.set_title(f"Recovering resistance status in Indian parasites\n(n = {n_india} QC-pass samples)",
                 fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[plot] wrote {out}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--release", default="Pf7",
                   help="Pf7 (built-in resistance-table URLs) or Pf8 (pass --class-url/--geno-url)")
    p.add_argument("--class-url", default=None)
    p.add_argument("--geno-url", default=None)
    p.add_argument("--country", default="India")
    p.add_argument("--max-markers", type=int, default=12)
    p.add_argument("--cv", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=Path, default=Path("results_india"))
    args = p.parse_args()

    meta, release = load_release(args.release)
    cls, gen = load_dr_tables(release, args.class_url, args.geno_url)

    mid, cid, gid = (find_id_column(meta), find_id_column(cls), find_id_column(gen))
    qc = find_qc_column(meta)
    country_col = find_country_column(meta)

    if qc:
        meta = meta[meta[qc].astype(bool)]
        print(f"[qc] {len(meta)} QC-pass samples")

    meta = meta.set_index(mid)
    cls = cls.set_index(cid)
    gen = gen.set_index(gid)
    shared = meta.index.intersection(cls.index).intersection(gen.index)
    meta, cls, gen = meta.loc[shared], cls.loc[shared], gen.loc[shared]
    print(f"[join] {len(shared)} samples with metadata, genotypes and classifications")

    is_india = meta[country_col] == args.country
    n_india = int(is_india.sum())
    print(f"[scope] {args.country}: {n_india} samples")
    print(meta.loc[meta[country_col].isin(SOUTH_ASIA), country_col].value_counts().to_string())
    if n_india < 50:
        print(f"\n!! Only {n_india} {args.country} samples. Too few for a country-specific "
              "panel. Re-run with --country covering South Asia, and report the wider "
              "scope honestly rather than over-reading a small national sample.\n",
              file=sys.stderr)

    marker_cols = [c for c in gen.columns if gen[c].dtype == object or gen[c].nunique() < 50]
    label_cols = [c for c in cls.columns if cls[c].nunique() <= 6 and cls[c].nunique() >= 2]
    print(f"[features] {len(marker_cols)} marker columns, {len(label_cols)} resistance labels")
    print(f"[labels] {label_cols}")

    args.out.mkdir(parents=True, exist_ok=True)

    # 1. informativeness, India vs global
    info = marker_table(gen[is_india], marker_cols, "india").merge(
        marker_table(gen, marker_cols, "global"), on="marker")
    info["uninformative_in_india_only"] = info.india_monomorphic & ~info.global_monomorphic
    info.sort_values("india_entropy_bits", ascending=False).to_csv(
        args.out / "marker_informativeness.tsv", sep="\t", index=False)
    n_mono = int(info.india_monomorphic.sum())
    n_mono_only = int(info.uninformative_in_india_only.sum())
    print(f"[info] {n_mono}/{len(info)} markers monomorphic in {args.country}; "
          f"{n_mono_only} of those ARE variable globally")

    # 2. panels: selected on India vs selected globally, both scored on India
    X_all, Y_all = encode(gen, marker_cols), cls[label_cols].astype(str)
    Xi, Yi = X_all[is_india.to_numpy()], Y_all[is_india.to_numpy()]

    print("\n--- forward selection on Indian samples ---")
    panel_india = greedy_panel(Xi, Yi, args.max_markers, args.cv, args.seed, "india")
    print("\n--- forward selection on all samples ---")
    panel_global = greedy_panel(X_all, Y_all, args.max_markers, args.cv, args.seed, "global")

    curves = {
        f"selected on {args.country}": [score_subset(Xi, Yi, [m for m, _ in panel_india[:k]],
                                                     args.cv, args.seed)
                                        for k in range(1, len(panel_india) + 1)],
        "selected globally": [score_subset(Xi, Yi, [m for m, _ in panel_global[:k]],
                                           args.cv, args.seed)
                              for k in range(1, len(panel_global) + 1)],
    }
    plot_curve(curves, args.out / "accuracy_vs_panel_size.png", n_india)

    pd.DataFrame({
        "rank": range(1, max(len(panel_india), len(panel_global)) + 1),
        f"panel_{args.country.lower()}": [m for m, _ in panel_india] + [None] * (
            max(len(panel_india), len(panel_global)) - len(panel_india)),
        "panel_global": [m for m, _ in panel_global] + [None] * (
            max(len(panel_india), len(panel_global)) - len(panel_global)),
    }).to_csv(args.out / "panel_comparison.tsv", sep="\t", index=False)

    # 3. regional heterogeneity within South Asia
    sa = meta[country_col].isin(SOUTH_ASIA)
    if sa.sum() > 0:
        reg = []
        for country, idx in meta[sa].groupby(country_col).groups.items():
            t = marker_table(gen.loc[idx], marker_cols, "x")
            t.insert(0, "country", country)
            t["n"] = len(idx)
            reg.append(t)
        pd.concat(reg).to_csv(args.out / "regional_heterogeneity.tsv", sep="\t", index=False)

    summary = {
        "release": release,
        "country": args.country,
        "n_samples_country": n_india,
        "n_samples_total": int(len(shared)),
        "n_markers_considered": len(marker_cols),
        "n_markers_monomorphic_in_country": n_mono,
        "n_markers_monomorphic_in_country_but_variable_globally": n_mono_only,
        "panel_country": [m for m, _ in panel_india],
        "panel_global": [m for m, _ in panel_global],
        "accuracy_curve_country_selected": [round(v, 4) for v in curves[f"selected on {args.country}"]],
        "accuracy_curve_global_selected": [round(v, 4) for v in curves["selected globally"]],
        "caveat": (
            "MalariaGEN resistance classifications are heuristic rules applied to "
            "genotypes, not measured phenotypes. This is panel design and feature "
            "selection, not independent phenotype prediction."
        ),
        "surveillance_note": (
            "A marker that is monomorphic in this country is uninformative for estimating "
            "current burden, but may still be worth retaining as a sentinel for imported "
            "or newly emerging resistance. Do not drop markers on these numbers alone."
        ),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[done] wrote {args.out}/")


if __name__ == "__main__":
    main()
