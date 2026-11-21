#!/usr/bin/env python3
"""Validate the panel-selection logic on synthetic data with a known answer.

Construction: 20 markers. Indian resistance status is driven by dhfr/dhps markers;
global status is additionally driven by kelch13 and crt markers that are MONOMORPHIC
in the Indian subset. A correct implementation should therefore
  (a) flag the kelch13/crt markers as monomorphic in India but variable globally, and
  (b) select a different, smaller panel for India than for the world.
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from india_panel import encode, greedy_panel, marker_table, score_subset

rng = np.random.default_rng(0)
N, N_IND = 2400, 700
AA = list("ACDEFGHIKLMNPQRSTVWY")

cols = ([f"dhfr_{p}" for p in (51, 59, 108, 164)]
        + [f"dhps_{p}" for p in (436, 437, 540, 581)]
        + [f"crt_{p}" for p in (72, 74, 76, 97, 220)]
        + [f"kelch13_{p}" for p in (446, 539, 561, 580)]
        + ["mdr1_86", "mdr1_184", "plasmepsin2-3_CN"])

is_india = np.zeros(N, dtype=bool)
is_india[:N_IND] = True

gen = pd.DataFrame(index=pd.Index([f"S{i}" for i in range(N)], name="Sample"))
for c in cols:
    wt, mut = rng.choice(AA, 2, replace=False)
    if c.startswith(("kelch13", "crt_97", "crt_220")):
        # variable globally, invariant in India  -> the thing we want detected
        g = np.where(rng.random(N) < 0.35, mut, wt)
        g[is_india] = wt
    else:
        g = np.where(rng.random(N) < 0.45, mut, wt)
    gen[c] = g

# Resistance rules. Pyrimethamine/sulfadoxine depend on dhfr/dhps (variable in India);
# artemisinin depends on kelch13 (invariant in India).
d108 = gen["dhfr_108"] != gen["dhfr_108"].mode()[0]
d51 = gen["dhfr_51"] != gen["dhfr_51"].mode()[0]
s437 = gen["dhps_437"] != gen["dhps_437"].mode()[0]
s540 = gen["dhps_540"] != gen["dhps_540"].mode()[0]
k580 = gen["kelch13_580"] != gen["kelch13_580"].mode()[0]
c76 = gen["crt_76"] != gen["crt_76"].mode()[0]

cls = pd.DataFrame({
    "Pyrimethamine": np.where(d108, "Resistant", "Sensitive"),
    "Sulfadoxine": np.where(s437 | s540, "Resistant", "Sensitive"),
    "SP_IPTp": np.where(d108 & d51 & s437, "Resistant", "Sensitive"),
    "Artemisinin": np.where(k580, "Resistant", "Sensitive"),
    "Chloroquine": np.where(c76, "Resistant", "Sensitive"),
}, index=gen.index)

print("=" * 66)
print("TEST 1: does it flag markers invariant in India but variable globally?")
print("=" * 66)
info = marker_table(gen[is_india], cols, "india").merge(marker_table(gen, cols, "global"), on="marker")
info["india_only_mono"] = info.india_monomorphic & ~info.global_monomorphic
flagged = sorted(info.loc[info.india_only_mono, "marker"])
expected = sorted([c for c in cols if c.startswith(("kelch13", "crt_97", "crt_220"))])
print(f"flagged : {flagged}")
print(f"expected: {expected}")
assert flagged == expected, "FAILED: wrong markers flagged"
print("PASS\n")

print("=" * 66)
print("TEST 2: does the India panel differ from the global panel?")
print("=" * 66)
X, Y = encode(gen, cols), cls.astype(str)
Xi, Yi = X[is_india], Y[is_india]
print("-- India --")
pi = greedy_panel(Xi, Yi, 6, 4, 0, "india")
print("-- Global --")
pg = greedy_panel(X, Y, 6, 4, 0, "global")
si, sg = [m for m, _ in pi], [m for m, _ in pg]
print(f"\nIndia  panel: {si}")
print(f"Global panel: {sg}")
assert si != sg, "FAILED: panels identical"
assert not any(m.startswith("kelch13") for m in si), "FAILED: picked a marker invariant in India"
assert any(m.startswith("kelch13") for m in sg), "FAILED: global panel missed kelch13"
print("PASS: India panel excludes kelch13; global panel includes it\n")

print("=" * 66)
print("TEST 3: does the India-selected panel beat the global one on Indian samples?")
print("=" * 66)
for k in (2, 3, 4):
    a = score_subset(Xi, Yi, si[:k], 4, 0)
    b = score_subset(Xi, Yi, sg[:k], 4, 0)
    print(f"k={k}:  India-selected {a:.4f}   global-selected {b:.4f}   delta {a - b:+.4f}")
assert score_subset(Xi, Yi, si[:3], 4, 0) >= score_subset(Xi, Yi, sg[:3], 4, 0), "FAILED"
print("PASS\n")
print("ALL TESTS PASSED")
