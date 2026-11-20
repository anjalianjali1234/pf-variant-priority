#!/usr/bin/env python3
"""Generate synthetic alignments for validating denv_severity_ml.py.

    python make_synthetic.py --out /tmp/denvtest            # signal present
    python make_synthetic.py --out /tmp/denvnull --no-signal # negative control
"""
import argparse
from pathlib import Path

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--out", type=Path, required=True)
p.add_argument("--n", type=int, default=120)
p.add_argument("--length", type=int, default=300)
p.add_argument("--causal", type=int, nargs="*", default=[50, 120, 200],
               help="0-based alignment columns carrying signal")
p.add_argument("--effect", type=float, default=0.8, help="alt frequency in cases")
p.add_argument("--no-signal", action="store_true")
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()

rng = np.random.default_rng(a.seed)
bases = np.array(list("ACGT"))
aln = np.tile(rng.choice(bases, a.length), (a.n, 1))

# background polymorphism; alt is always forced to differ from the consensus base
for j in rng.choice(a.length, max(1, a.length // 8), replace=False):
    alt = rng.choice([b for b in bases if b != aln[0, j]])
    aln[rng.random(a.n) < 0.3, j] = alt

y = rng.integers(0, 2, a.n)
if not a.no_signal:
    for j in a.causal:
        alt = rng.choice([b for b in bases if b != aln[0, j]])
        aln[(y == 1) & (rng.random(a.n) < a.effect), j] = alt
        aln[(y == 0) & (rng.random(a.n) < 0.1), j] = alt

a.out.mkdir(parents=True, exist_ok=True)
with open(a.out / "aln.fasta", "w") as fh:
    for i in range(a.n):
        fh.write(f">s{i}\n{''.join(aln[i])}\n")
with open(a.out / "meta.tsv", "w") as fh:
    fh.write("sample\tseverity\n")
    for i in range(a.n):
        fh.write(f"s{i}\t{'severe' if y[i] else 'mild'}\n")

label = "no signal" if a.no_signal else f"signal at 1-based columns {[c + 1 for c in a.causal]}"
print(f"wrote {a.out}/aln.fasta and meta.tsv — {a.n} genomes, {label}")
