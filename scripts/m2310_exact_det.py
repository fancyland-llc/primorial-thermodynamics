"""
m2310_exact_det.py — Exact integer determinant for m=2310
=========================================================

The sanity check showed determinants are ~154 digits (not ~900).
This might finish in minutes, not hours. Let's try.

Prediction (to be tested): ratio = -3^5 * 23 = -5589
Numerical estimate: ratio ≈ -108927 (float, unreliable)

If this takes too long, Ctrl+C — the blocks are saved as .npy files.
"""

import numpy as np
from sympy import Matrix, factorint
import time
import sys

print("Loading 240x240 integer blocks...")
# Rebuild from saved pairs (npy with object dtype may not load cleanly)
from math import gcd

m = 2310
coprimes = sorted([x for x in range(1, m) if gcd(x, m) == 1])
pairs = [(x, m - x) for x in coprimes if x < m - x]
idx = {x: i for i, x in enumerate(coprimes)}

def cyclic_distance(a, b, n):
    d = abs(a - b) % n
    return min(d, n - d)

print("Building D_sym(2310) [480x480]...")
D = [[0]*480 for _ in range(480)]
for i in range(480):
    for j in range(i+1, 480):
        d = cyclic_distance(coprimes[i], coprimes[j], m)
        D[i][j] = d
        D[j][i] = d

print("Building 240x240 integer blocks...")
D_even_int = []
D_odd_int = []
for k in range(240):
    row_e, row_o = [], []
    xk, mxk = pairs[k]
    ik, jk = idx[xk], idx[mxk]
    for l in range(240):
        xl, mxl = pairs[l]
        il, jl = idx[xl], idx[mxl]
        se = D[ik][il] + D[ik][jl] + D[jk][il] + D[jk][jl]
        so = D[ik][il] - D[ik][jl] - D[jk][il] + D[jk][jl]
        row_e.append(se // 2)
        row_o.append(so // 2)
    D_even_int.append(row_e)
    D_odd_int.append(row_o)

print("Blocks built. Starting exact determinant computation...")
print("=" * 60)

# ── D_odd first (it's negative definite → likely faster/simpler) ──
print("\nComputing det(D_odd) [240x240 exact integer]...")
print(f"  Started at: {time.strftime('%H:%M:%S')}")
sys.stdout.flush()

t0 = time.time()
M_odd = Matrix(D_odd_int)
det_odd = M_odd.det()
t1 = time.time()

print(f"  Finished in {t1-t0:.1f}s")
print(f"  det(D_odd) = ...{str(det_odd)[-40:]} ({len(str(abs(det_odd)))} digits)")
print(f"  sign = {'+'  if det_odd > 0 else '-'}")
sys.stdout.flush()

# ── D_even ──
print(f"\nComputing det(D_even) [240x240 exact integer]...")
print(f"  Started at: {time.strftime('%H:%M:%S')}")
sys.stdout.flush()

t0 = time.time()
M_even = Matrix(D_even_int)
det_even = M_even.det()
t1 = time.time()

print(f"  Finished in {t1-t0:.1f}s")
print(f"  det(D_even) = ...{str(det_even)[-40:]} ({len(str(abs(det_even)))} digits)")
print(f"  sign = {'+' if det_even > 0 else '-'}")
sys.stdout.flush()

# ── Ratio ──
print("\n" + "=" * 60)
print("EXACT RESULTS")
print("=" * 60)

if det_odd != 0 and det_even % det_odd == 0:
    ratio = det_even // det_odd
    print(f"Ratio = det(D_even)/det(D_odd) = {ratio}")
    print(f"Ratio is a clean integer: YES")
    
    # Power of 3
    r = abs(ratio)
    pow3 = 0
    while r % 3 == 0:
        r //= 3
        pow3 += 1
    print(f"Power of 3 in |ratio|: 3^{pow3}")
    print(f"Remaining factor: {r if ratio > 0 else -r}")
    print(f"Sign of ratio: {'+' if ratio > 0 else '-'}")
    
    # p_max exclusion
    print(f"\np_max = 11 exclusion check:")
    print(f"  det(D_even) % 11 = {det_even % 11}")
    print(f"  det(D_odd) % 11 = {det_odd % 11}")
    if det_even % 11 != 0 and det_odd % 11 != 0:
        print(f"  11 is ABSENT from both determinants — EXCLUSION CONFIRMED")
    elif det_even % 11 == 0 or det_odd % 11 == 0:
        print(f"  11 DIVIDES at least one determinant — EXCLUSION BROKEN")
    
    # Check prediction
    print(f"\nPrediction was: ratio = -3^5 * 23 = -5589")
    print(f"Actual ratio:   {ratio}")
    if ratio == -5589:
        print("PREDICTION CONFIRMED!")
    else:
        print(f"PREDICTION WRONG — actual pattern revealed")
        # Try to factor the ratio (small number, instant)
        print(f"Factorization of |ratio|: {factorint(abs(ratio))}")
else:
    print(f"det(D_odd) = {det_odd}")
    print(f"det(D_even) = {det_even}")
    if det_odd == 0:
        print("D_odd is SINGULAR — unexpected!")
    else:
        print(f"Ratio is NOT a clean integer")
        print(f"det(D_even) mod det(D_odd) = {det_even % det_odd}")

# Save exact values
print("\nSaving exact determinant values...")
with open("det_2310_exact.txt", "w") as f:
    f.write(f"det_even = {det_even}\n")
    f.write(f"det_odd = {det_odd}\n")
    if det_odd != 0 and det_even % det_odd == 0:
        f.write(f"ratio = {det_even // det_odd}\n")
print("Saved to det_2310_exact.txt")
