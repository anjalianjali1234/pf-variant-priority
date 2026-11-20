#!/usr/bin/env python3
"""Combine constraint and population evidence into a ranked VUS watchlist.

Four signals, each scaled to [0, 1] across the scored variants, then combined with the
hand-set weights in conf/priority.yml:

  constraint  -esm_score      more constrained substitutions rank higher
  frequency   allele_freq     common variants are easier to follow up
  spread      n_sites         seen at many sites is harder to explain as drift
  trend       recent rise     frequency rising across collection years

Rank stability is estimated by perturbing the weights (Dirichlet noise) and recording the
spread of each variant's rank. A variant whose rank swings wildly under small weight
changes is not a robust call, and the output says so.

Usage:
    python bin/prioritise.py --scored results/scored.tsv \
        --config conf/priority.yml --out results/watchlist.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def unit_scale(s: pd.Series) -> pd.Series:
    """Min-max to [0,1]; constant input maps to 0.5 (carries no information)."""
    s = s.astype(float)
    lo, hi = s.min(), s.max()
    return pd.Series(0.5, index=s.index) if not np.isfinite(lo) or hi == lo else (s - lo) / (hi - lo)


def temporal_trend(df: pd.DataFrame) -> pd.Series:
    """Crude rise indicator: first observed recently and spanning >1 year.

    Real trend estimation needs per-year counts. Where only first/last year are present
    this returns a coarse proxy, and 0 when dates are absent.
    """
    if not {"first_year", "last_year"}.issubset(df.columns):
        return pd.Series(0.0, index=df.index)
    first = pd.to_numeric(df["first_year"], errors="coerce")
    last = pd.to_numeric(df["last_year"], errors="coerce")
    span = (last - first).clip(lower=0)
    recency = first - first.min()
    return (unit_scale(recency.fillna(0)) * (span > 0).astype(float)).fillna(0.0)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scored", type=Path, required=True)
    p.add_argument("--config", type=Path, default=Path("conf/priority.yml"))
    p.add_argument("--validated", type=Path, default=Path("data/who_validated_k13.tsv"))
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    w = cfg["weights"]
    filt = cfg.get("filters", {})
    boot = cfg.get("bootstrap", {})

    df = pd.read_csv(args.scored, sep="\t", comment="#")
    n_in = len(df)

    in_lcr = df.get("in_lcr", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    lcr = df[in_lcr].copy()
    if filt.get("separate_lcr", True):
        df = df[~in_lcr].copy()
    df = df[df.esm_score.notna()].copy()
    if df.empty:
        sys.exit("no scored variants to prioritise — check score_variants.py output")

    if filt.get("exclude_validated", True) and args.validated.exists():
        v = pd.read_csv(args.validated, sep="\t", comment="#")
        known = set(zip(v.aa_ref, v.aa_pos, v.aa_alt)) if "status" in v else set()
        before = len(df)
        df = df[[k not in known for k in zip(df.aa_ref, df.aa_pos, df.aa_alt)]].copy()
        print(f"[prioritise] removed {before - len(df)} already-catalogued markers", file=sys.stderr)

    min_s = int(filt.get("min_samples", 1))
    if "n_samples" in df.columns:
        before = len(df)
        df = df[df.n_samples >= min_s].copy()
        print(f"[prioritise] removed {before - len(df)} variants seen in < {min_s} samples", file=sys.stderr)

    if df.empty:
        sys.exit("all variants filtered out — loosen conf/priority.yml")

    comp = pd.DataFrame(index=df.index)
    comp["constraint"] = unit_scale(-df.esm_score)
    comp["frequency"] = unit_scale(df.allele_freq) if "allele_freq" in df else 0.0
    comp["spread"] = unit_scale(df.n_sites) if "n_sites" in df else 0.0
    comp["trend"] = temporal_trend(df)

    keys = [c for c in ("constraint", "frequency", "spread", "trend") if c in w]
    base = np.array([w[c] for c in keys], dtype=float)
    base = base / base.sum()
    df["priority_score"] = (comp[keys].to_numpy() @ base).round(4)
    df = df.sort_values("priority_score", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))

    # Rank stability: perturb the weights and record how far each rank moves.
    comp_sorted = pd.DataFrame({
        "constraint": unit_scale(-df.esm_score),
        "frequency": unit_scale(df.allele_freq) if "allele_freq" in df else 0.0,
        "spread": unit_scale(df.n_sites) if "n_sites" in df else 0.0,
        "trend": temporal_trend(df),
    })[keys].to_numpy()

    rng = np.random.default_rng(boot.get("seed", 42))
    n_iter = int(boot.get("n_iterations", 1000))
    ranks = np.empty((n_iter, len(df)), dtype=int)
    for i in range(n_iter):
        perturbed = rng.dirichlet(base * 20.0)
        ranks[i] = (-(comp_sorted @ perturbed)).argsort().argsort() + 1

    df["rank_p05"] = np.percentile(ranks, 5, axis=0).astype(int)
    df["rank_p95"] = np.percentile(ranks, 95, axis=0).astype(int)
    df["rank_stable"] = (df.rank_p95 - df.rank_p05) <= max(2, len(df) // 10)

    for k in keys:
        df[f"component_{k}"] = comp_sorted[:, keys.index(k)].round(3)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# Ranked watchlist of variants of unknown significance.\n"
        "# This is a TRIAGE ORDERING FOR FOLLOW-UP, not a prediction of drug resistance.\n"
        f"# Weights: {', '.join(f'{k}={v:.3f}' for k, v in zip(keys, base))}"
        f"  (from {args.config})\n"
        f"# {n_in} input variants -> {len(df)} ranked; {len(lcr)} in low-complexity regions "
        "reported separately.\n"
        "# Read benchmark.json first: if the score is uncalibrated, this ordering is exploratory.\n"
    )
    with open(args.out, "w") as fh:
        fh.write(header)
        df.to_csv(fh, sep="\t", index=False)
    print(f"[prioritise] wrote {args.out} ({len(df)} variants, "
          f"{int(df.rank_stable.sum())} with stable ranks)", file=sys.stderr)

    if len(lcr):
        lcr_path = args.out.with_name(args.out.stem + "_low_complexity.tsv")
        lcr.to_csv(lcr_path, sep="\t", index=False)
        print(f"[prioritise] wrote {lcr_path} ({len(lcr)} unranked LCR variants)", file=sys.stderr)


if __name__ == "__main__":
    main()
