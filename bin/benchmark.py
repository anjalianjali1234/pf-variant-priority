#!/usr/bin/env python3
"""Calibrate the ESM score against WHO-validated kelch13 markers.

The question this answers: does an unsupervised constraint score recover mutations we
already know matter? If it cannot, the ranking it produces over variants of unknown
significance should not be trusted, and the report says so explicitly.

Design:
  positives  WHO-validated kelch13 markers present in the scored table
  negatives  other missense kelch13 variants observed in the same data that are neither
             validated nor WHO candidates
  metric     ROC-AUC with a bootstrap confidence interval, precision@k, and the rank of
             each individual validated marker

Caveat, stated in the output: the negative set is "not yet known to matter", not "known
not to matter". Some negatives may be undiscovered resistance alleles, which biases AUC
downwards. With ~13 positives in one structural domain the estimate is indicative only.

Usage:
    python bin/benchmark.py --scored results/scored.tsv \
        --validated data/who_validated_k13.tsv --out results/benchmark/
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load_markers(path: Path) -> tuple[set, set]:
    df = pd.read_csv(path, sep="\t", comment="#")
    key = lambda d: set(zip(d["aa_ref"], d["aa_pos"], d["aa_alt"]))
    return key(df[df.status == "validated"]), key(df[df.status == "candidate"])


def bootstrap_auc(y: np.ndarray, s: np.ndarray, n_iter: int, seed: int) -> tuple[float, float]:
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_iter):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) == 2:
            vals.append(roc_auc_score(y[idx], s[idx]))
    if not vals:
        return (float("nan"), float("nan"))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scored", type=Path, required=True)
    p.add_argument("--validated", type=Path, default=Path("data/who_validated_k13.tsv"))
    p.add_argument("--gene", default="pfk13", help="gene to benchmark within")
    p.add_argument("--n-bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    from sklearn.metrics import roc_auc_score, roc_curve

    df = pd.read_csv(args.scored, sep="\t", comment="#")
    df = df[(df.gene == args.gene) & df.esm_score.notna()].copy()
    if df.empty:
        sys.exit(f"no scored {args.gene} variants found in {args.scored}")

    validated, candidate = load_markers(args.validated)
    keys = list(zip(df.aa_ref, df.aa_pos, df.aa_alt))
    df["is_validated"] = [k in validated for k in keys]
    df["is_candidate"] = [k in candidate for k in keys]

    bench = df[df.is_validated | ~df.is_candidate].copy()   # drop candidates: status unclear
    y = bench.is_validated.to_numpy().astype(int)
    # More constrained (more negative score) should rank higher, so negate.
    s = -bench.esm_score.to_numpy()

    args.out.mkdir(parents=True, exist_ok=True)
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())

    report: dict = {
        "gene": args.gene,
        "n_positives": n_pos,
        "n_negatives": n_neg,
        "n_candidates_excluded": int(df.is_candidate.sum()),
    }

    if n_pos < 2 or n_neg < 2:
        report["status"] = "INSUFFICIENT_DATA"
        report["message"] = (
            f"Only {n_pos} positives and {n_neg} negatives were scored. No AUC is reported. "
            "The prioritised watchlist is UNCALIBRATED and must be read as exploratory."
        )
        print(f"[benchmark] {report['message']}", file=sys.stderr)
    else:
        auc = float(roc_auc_score(y, s))
        lo, hi = bootstrap_auc(y, s, args.n_bootstrap, args.seed)
        bench["rank"] = bench.esm_score.rank(method="min").astype(int)
        k = max(1, n_pos)
        prec_at_k = float(bench.nsmallest(k, "esm_score").is_validated.mean())

        report.update({
            "status": "ok",
            "roc_auc": round(auc, 3),
            "roc_auc_ci95": [round(lo, 3), round(hi, 3)],
            f"precision_at_{k}": round(prec_at_k, 3),
            "validated_marker_ranks": {
                f"{r.aa_ref}{r.aa_pos}{r.aa_alt}": {"rank": int(r.rank), "of": len(bench),
                                                    "esm_score": float(r.esm_score)}
                for r in bench[bench.is_validated].itertuples()
            },
            "interpretation": (
                "AUC near 0.5 means the constraint score carries no signal for this gene and "
                "the watchlist should not be used for triage. AUC above ~0.7 with a CI "
                "excluding 0.5 supports cautious use."
            ),
            "caveat": (
                "Negatives are 'not yet known to matter', not 'known not to matter'. "
                f"With {n_pos} positives confined to one domain this estimate is indicative."
            ),
        })
        print(f"[benchmark] AUC={auc:.3f} (95% CI {lo:.3f}-{hi:.3f}), "
              f"{n_pos} positives vs {n_neg} negatives", file=sys.stderr)

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fpr, tpr, _ = roc_curve(y, s)
            fig, ax = plt.subplots(figsize=(4.5, 4.5))
            ax.plot(fpr, tpr, lw=2, label=f"ESM-2  AUC={auc:.2f}")
            ax.plot([0, 1], [0, 1], "k--", lw=1, label="chance")
            ax.set_xlabel("False positive rate")
            ax.set_ylabel("True positive rate")
            ax.set_title(f"Recovery of WHO-validated {args.gene} markers\n"
                         f"({n_pos} positives, {n_neg} negatives)", fontsize=10)
            ax.legend(frameon=False, fontsize=9)
            fig.tight_layout()
            fig.savefig(args.out / "roc.png", dpi=200)
            print(f"[benchmark] wrote {args.out / 'roc.png'}", file=sys.stderr)
        except Exception as exc:  # plotting must never break the run
            print(f"[benchmark] plot skipped: {exc}", file=sys.stderr)

    (args.out / "benchmark.json").write_text(json.dumps(report, indent=2))
    bench.to_csv(args.out / "benchmark_table.tsv", sep="\t", index=False)
    print(f"[benchmark] wrote {args.out / 'benchmark.json'}", file=sys.stderr)


if __name__ == "__main__":
    main()
