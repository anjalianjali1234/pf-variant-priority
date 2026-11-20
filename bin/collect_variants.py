#!/usr/bin/env python3
"""Parse bcftools csq-annotated VCFs into a tidy protein-variant table.

Reads one or more VCF files annotated with `bcftools csq` (which writes a BCSQ INFO
field), extracts missense consequences falling inside the resistance panel, and
aggregates them across samples into the table consumed by score_variants.py.

The BCSQ format is pipe-delimited:
    consequence|gene|transcript|biotype|strand|amino_acid_change|nucleotide_change
e.g. missense|PF3D7_1343700|PF3D7_1343700.1|protein_coding|+|580C>Y|1738T>A

Usage:
    python bin/collect_variants.py --vcf 'results/vcf/*.csq.vcf.gz' \
        --panel data/pfmdr15_panel.tsv --out results/variant_table.tsv
"""
from __future__ import annotations

import argparse
import glob
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# 580C>Y  or  Met1? style; we want ref-pos-alt for simple substitutions
AA_CHANGE = re.compile(r"^(\d+)([A-Z*])>([A-Z*])$")


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def load_panel(path: Path) -> dict[str, str]:
    """Return {plasmodb_id: gene_symbol} from the panel TSV (comment lines start with #)."""
    df = pd.read_csv(path, sep="\t", comment="#")
    missing = {"gene_symbol", "plasmodb_id"} - set(df.columns)
    if missing:
        sys.exit(f"panel file is missing column(s): {sorted(missing)}")
    return dict(zip(df["plasmodb_id"], df["gene_symbol"]))


def parse_bcsq(info: str) -> list[tuple[str, str, int, str]]:
    """Extract (gene_id, aa_ref, aa_pos, aa_alt) tuples from a VCF INFO string."""
    out = []
    for field in info.split(";"):
        if not field.startswith("BCSQ="):
            continue
        for entry in field[len("BCSQ="):].split(","):
            parts = entry.split("|")
            if len(parts) < 6 or "missense" not in parts[0]:
                continue
            gene_id = parts[1]
            m = AA_CHANGE.match(parts[5])
            if m:
                out.append((gene_id, m.group(2), int(m.group(1)), m.group(3)))
    return out


def collect(vcf_paths: list[str], panel: dict[str, str], min_depth: int) -> pd.DataFrame:
    counts: dict[tuple, set[str]] = defaultdict(set)
    samples_seen: set[str] = set()
    skipped_depth = 0

    for path in vcf_paths:
        sample = Path(path).name.split(".")[0]
        samples_seen.add(sample)
        with open_maybe_gzip(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 8 or cols[6] not in ("PASS", "."):
                    continue
                info = cols[7]
                dp = re.search(r"(?:^|;)DP=(\d+)", info)
                if dp and int(dp.group(1)) < min_depth:
                    skipped_depth += 1
                    continue
                for gene_id, ref, pos, alt in parse_bcsq(info):
                    if gene_id in panel:
                        counts[(gene_id, ref, pos, alt)].add(sample)

    n_total = len(samples_seen)
    if n_total == 0:
        sys.exit("no VCFs read — check the --vcf glob")

    rows = []
    for (gene_id, ref, pos, alt), samples in counts.items():
        rows.append(
            {
                "gene": panel[gene_id],
                "plasmodb_id": gene_id,
                "aa_ref": ref,
                "aa_pos": pos,
                "aa_alt": alt,
                "n_samples": len(samples),
                "n_total": n_total,
                "allele_freq": round(len(samples) / n_total, 4),
            }
        )

    df = pd.DataFrame(rows).sort_values(["gene", "aa_pos"]).reset_index(drop=True)
    print(
        f"[collect] {len(vcf_paths)} VCFs, {n_total} samples, "
        f"{len(df)} distinct missense variants, {skipped_depth} calls below depth {min_depth}",
        file=sys.stderr,
    )
    return df


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="glob for bcftools csq-annotated VCF(s)")
    p.add_argument("--panel", type=Path, default=Path("data/pfmdr15_panel.tsv"))
    p.add_argument("--min-depth", type=int, default=20, help="minimum DP to keep a call")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    paths = sorted(glob.glob(args.vcf))
    if not paths:
        sys.exit(f"no files matched {args.vcf!r}")

    df = collect(paths, load_panel(args.panel), args.min_depth)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[collect] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
