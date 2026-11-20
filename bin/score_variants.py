#!/usr/bin/env python3
"""Score missense variants with ESM-2 masked-marginal log-odds.

For each variant (wt, pos, mut) the residue at `pos` is masked and the model's
log-probabilities at that position are read off:

    score = log P(mut | x_masked) - log P(wt | x_masked)

Negative scores mean the substitution is less plausible than wild type given the
sequence context, which is used here as a proxy for evolutionary constraint. This is
unsupervised: no resistance labels enter the calculation.

Two details matter for Plasmodium and are handled explicitly:

* **Long proteins.** ESM-2 accepts ~1022 residues. Several panel proteins (pfmdr1 ~1419 aa,
  pfubp1 ~3800 aa) exceed this, so a window is centred on the variant position and the
  score computed within it. Window size and the variant's offset are recorded.
* **Low-complexity regions.** P. falciparum proteins carry long asparagine/lysine repeats
  that language models handle poorly. Variants inside them are flagged, not silently
  ranked.

Usage:
    python bin/score_variants.py --variants data/example_variant_table.tsv \
        --fasta data/panel_proteins.faa --out results/scored.tsv
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

MAX_LEN = 1022  # ESM-2 context minus BOS/EOS


def read_fasta(path: Path) -> dict[str, str]:
    """Minimal FASTA reader. Key is the first whitespace-delimited token of the header."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    if not seqs:
        sys.exit(f"no sequences read from {path}")
    return seqs


def is_low_complexity(seq: str, pos: int, window: int, threshold: float) -> bool:
    """True if a window centred on pos (1-based) is dominated by a single residue."""
    half = window // 2
    lo, hi = max(0, pos - 1 - half), min(len(seq), pos - 1 + half + 1)
    sub = seq[lo:hi]
    if len(sub) < 5:
        return False
    return Counter(sub).most_common(1)[0][1] / len(sub) >= threshold


def make_window(seq: str, pos: int) -> tuple[str, int]:
    """Return (subsequence, 1-based position of the variant inside it)."""
    if len(seq) <= MAX_LEN:
        return seq, pos
    half = MAX_LEN // 2
    start = max(0, min(pos - 1 - half, len(seq) - MAX_LEN))
    return seq[start : start + MAX_LEN], pos - start


def score_all(df: pd.DataFrame, seqs: dict[str, str], checkpoint: str,
              mask_lcr: bool, lcr_window: int, lcr_threshold: float) -> pd.DataFrame:
    import torch
    from transformers import AutoTokenizer, EsmForMaskedLM

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[score] loading {checkpoint} on {device}", file=sys.stderr)
    tok = AutoTokenizer.from_pretrained(checkpoint)
    model = EsmForMaskedLM.from_pretrained(checkpoint).to(device).eval()

    records = []
    for row in df.itertuples(index=False):
        key = row.plasmodb_id
        seq = seqs.get(key)
        rec = {
            "esm_score": None, "in_lcr": None, "window_len": None, "qc": "ok",
        }

        if seq is None:
            rec["qc"] = "sequence_not_found"
        elif row.aa_pos > len(seq):
            rec["qc"] = f"position_beyond_sequence_len_{len(seq)}"
        elif seq[row.aa_pos - 1] != row.aa_ref:
            # A mismatch means the variant table and the FASTA disagree — usually a
            # transcript-version or coordinate problem. Never silently score it.
            rec["qc"] = f"ref_mismatch_fasta_has_{seq[row.aa_pos - 1]}"
        else:
            rec["in_lcr"] = is_low_complexity(seq, row.aa_pos, lcr_window, lcr_threshold)
            if mask_lcr and rec["in_lcr"]:
                rec["qc"] = "low_complexity_not_scored"
            else:
                sub, p = make_window(seq, row.aa_pos)
                rec["window_len"] = len(sub)
                enc = tok(sub, return_tensors="pt").to(device)
                enc["input_ids"][0, p] = tok.mask_token_id  # +1 offset for BOS
                with torch.no_grad():
                    logits = model(**enc).logits
                logprobs = torch.log_softmax(logits[0, p], dim=-1)
                rec["esm_score"] = round(
                    (logprobs[tok.convert_tokens_to_ids(row.aa_alt)]
                     - logprobs[tok.convert_tokens_to_ids(row.aa_ref)]).item(), 4
                )
        records.append(rec)

    out = pd.concat([df.reset_index(drop=True), pd.DataFrame(records)], axis=1)
    n_ok = out["esm_score"].notna().sum()
    print(f"[score] scored {n_ok}/{len(out)} variants", file=sys.stderr)
    for reason, n in out.loc[out["esm_score"].isna(), "qc"].value_counts().items():
        print(f"[score]   skipped {n}: {reason}", file=sys.stderr)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variants", type=Path, required=True)
    p.add_argument("--fasta", type=Path, required=True, help="panel protein FASTA from PlasmoDB")
    p.add_argument("--checkpoint", default="facebook/esm2_t33_650M_UR50D",
                   help="use facebook/esm2_t12_35M_UR50D for a quick laptop run")
    p.add_argument("--mask-lcr", action="store_true", default=True)
    p.add_argument("--no-mask-lcr", dest="mask_lcr", action="store_false")
    p.add_argument("--lcr-window", type=int, default=20)
    p.add_argument("--lcr-threshold", type=float, default=0.6)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    df = pd.read_csv(args.variants, sep="\t", comment="#")
    required = {"plasmodb_id", "aa_ref", "aa_pos", "aa_alt"}
    if missing := required - set(df.columns):
        sys.exit(f"variant table is missing column(s): {sorted(missing)}")

    scored = score_all(df, read_fasta(args.fasta), args.checkpoint,
                       args.mask_lcr, args.lcr_window, args.lcr_threshold)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.out, sep="\t", index=False)
    print(f"[score] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
