"""
v3_overnight.py — Overnight mod-p^k computation for m=30030.

Runs vectorized Gaussian elimination mod 3^k for k=1..8 on both
D_odd(30030) and D_even(30030) (2880x2880 matrices).

Validates approach on m=210 (24x24) and m=2310 (240x240) first.

KEY INSIGHT (proven before running):
  GF(3) nullity is 943 for BOTH blocks, but v_3(R) = -1.
  Since v_3(det_even) >= 943 and v_3(det_odd) >= 943,
  and v_3(det_even) - v_3(det_odd) = -1,
  we MUST have v_3(det_odd) >= 944.
  The "odd block GF(3) is exact" pattern BREAKS at m=30030.

  gap_even - gap_odd = -1, so gap_odd = gap_even + 1.
  Minimum: gap_even=0, gap_odd=1 => v_3(det_even)=943, v_3(det_odd)=944.

Output saved to v3_overnight_results.txt

Author: Antonio P. Matos / Fancyland LLC
Date: March 23, 2026
"""

import numpy as np
import math
import time
import sys
import os
from datetime import datetime

# ============================================================
# Output tee — write to both console and file
# ============================================================

RESULTS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "v3_overnight_results.txt")

class Tee:
    def __init__(self, filepath):
        self.terminal = sys.stdout
        self.log = open(filepath, 'w', encoding='utf-8')
    def write(self, msg):
        self.terminal.write(msg)
        self.log.write(msg)
        self.log.flush()
    def flush(self):
        self.terminal.flush()
        self.log.flush()

sys.stdout = Tee(RESULTS_FILE)

print(f"v3_overnight.py — started {datetime.now().isoformat()}")
print(f"Results saved to: {RESULTS_FILE}")
print()

# ============================================================
# Infrastructure
# ============================================================

def v_p(n, p):
    """p-adic valuation of integer n."""
    if n == 0:
        return float('inf')
    n = abs(int(n))
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def build_blocks_numpy(m):
    """Build D_even and D_odd as int64 numpy arrays."""
    coprimes = [r for r in range(1, m) if math.gcd(r, m) == 1]
    pairs = np.array([r for r in coprimes if r < m / 2], dtype=np.int64)
    n = len(pairs)
    ri = pairs.reshape(-1, 1)
    rj = pairs.reshape(1, -1)
    diff = (rj - ri) % m
    d_diff = np.minimum(diff, m - diff)
    summ = (ri + rj) % m
    d_sum = np.minimum(summ, m - summ)
    D_even = (d_diff + d_sum).astype(np.int64)
    D_odd = (d_diff - d_sum).astype(np.int64)
    return D_even, D_odd, n

def v3_gauss_mod_pk(M_int, p, k, label=""):
    """
    Vectorized Gaussian elimination mod p^k with p-adic pivot selection.

    Returns:
      total_v: sum of min(v_p(d_i), k) over all SNF diagonal entries d_i.
               This is a lower bound on v_p(det(M)) that improves with k.
      pivots_by_v: dict mapping valuation -> count of pivots at that valuation
      elapsed: wall-clock seconds
    """
    pk = p ** k
    n = M_int.shape[0]
    M = M_int.astype(np.int64) % pk
    M = ((M % pk) + pk) % pk

    total_v = 0
    row_current = 0
    pivots_by_v = {}
    t0 = time.time()

    for col in range(n):
        if col % 500 == 0 and col > 0:
            elapsed = time.time() - t0
            rate = col / elapsed
            eta = (n - col) / rate if rate > 0 else 0
            print(f"    [{label} k={k}] col {col}/{n}, "
                  f"v_so_far={total_v}, {elapsed:.0f}s, ~{eta:.0f}s left",
                  flush=True)

        # --- Pivot selection: find min p-adic valuation in column ---
        col_below = M[row_current:, col]
        nonzero_idx = np.nonzero(col_below)[0]

        if len(nonzero_idx) == 0:
            # Entire column is 0 mod p^k => valuation >= k
            total_v += k
            pivots_by_v[k] = pivots_by_v.get(k, 0) + 1
            continue

        best_v = k
        best_local = -1
        for idx in nonzero_idx:
            entry = int(col_below[idx])
            v = 0
            while entry % p == 0:
                entry //= p
                v += 1
            if v < best_v:
                best_v = v
                best_local = idx
                if v == 0:
                    break  # can't improve

        total_v += best_v
        pivots_by_v[best_v] = pivots_by_v.get(best_v, 0) + 1

        actual_pivot = int(best_local) + row_current
        if actual_pivot != row_current:
            M[[row_current, actual_pivot]] = M[[actual_pivot, row_current]].copy()

        # --- Elimination ---
        pivot_val = int(M[row_current, col])
        pv = p ** best_v
        pivot_unit = (pivot_val // pv) % (pk // pv) if pv < pk else 1

        if pivot_unit % p == 0:
            # Non-invertible unit part — can't eliminate cleanly
            row_current += 1
            continue

        remaining_pk = pk // pv
        if remaining_pk <= 1:
            row_current += 1
            continue

        try:
            pivot_inv = pow(int(pivot_unit), -1, int(remaining_pk))
        except (ValueError, ZeroDivisionError):
            row_current += 1
            continue

        # Vectorized row elimination
        below = M[row_current + 1:, col].copy()
        nonzero_mask = below != 0

        if np.any(nonzero_mask):
            nz_rows = np.where(nonzero_mask)[0]
            entries = below[nz_rows]
            entries_div = entries // pv
            factors = (entries_div * pivot_inv) % remaining_pk

            actual_rows = nz_rows + row_current + 1
            pivot_row = M[row_current, col + 1:].astype(np.int64)

            # Broadcasting: factors[:, None] * pivot_row[None, :]
            # Values bounded by (remaining_pk-1) * (pk-1), fits in int64
            updates = (factors[:, None] * pivot_row[None, :]) % pk
            M_sub = M[actual_rows][:, col + 1:]
            M[actual_rows, col + 1:] = (M_sub - updates + pk) % pk
            M[actual_rows, col] = 0

        row_current += 1

    elapsed = time.time() - t0
    return total_v, pivots_by_v, elapsed


# ============================================================
# VALIDATION on m=210 and m=2310
# ============================================================

print("=" * 70)
print("VALIDATION: mod-3^k vs exact v_3 at m=210 and m=2310")
print("=" * 70)

# --- m=210 (24x24) ---
print("\n--- m=210 (24x24) ---")
print("  Known: v_3(det_even)=11, v_3(det_odd)=7")
De210, Do210, n210 = build_blocks_numpy(210)

for block_name, block_mat, known_v3 in [("even", De210, 11), ("odd", Do210, 7)]:
    print(f"\n  {block_name} block (true v_3 = {known_v3}):")
    for k in range(1, 9):
        total, pivots, elapsed = v3_gauss_mod_pk(block_mat, 3, k,
                                                  label=f"210_{block_name}")
        converged = "  <-- EXACT" if total == known_v3 else ""
        print(f"    k={k} (mod {3**k:>5d}): total_v = {total:>4d}  "
              f"[{elapsed:.2f}s]{converged}")
        if total == known_v3:
            break

# --- m=2310 (240x240) ---
print("\n--- m=2310 (240x240) ---")
print("  Known: v_3(det_even)=77, v_3(det_odd)=75")
De2310, Do2310, n2310 = build_blocks_numpy(2310)

for block_name, block_mat, known_v3 in [("even", De2310, 77), ("odd", Do2310, 75)]:
    print(f"\n  {block_name} block (true v_3 = {known_v3}):")
    for k in range(1, 9):
        total, pivots, elapsed = v3_gauss_mod_pk(block_mat, 3, k,
                                                  label=f"2310_{block_name}")
        converged = "  <-- EXACT" if total == known_v3 else ""
        print(f"    k={k} (mod {3**k:>5d}): total_v = {total:>4d}  "
              f"pivots_by_v = {dict(sorted(pivots.items()))}  "
              f"[{elapsed:.2f}s]{converged}")
        if total == known_v3:
            break

# ============================================================
# MAIN COMPUTATION: m=30030 (2880x2880)
# ============================================================

print(f"\n{'=' * 70}")
print("MAIN COMPUTATION: m=30030 (2880x2880)")
print("=" * 70)

print("\nBuilding D_even(30030) and D_odd(30030)...")
t_build = time.time()
De30030, Do30030, n30030 = build_blocks_numpy(30030)
print(f"Built in {time.time() - t_build:.1f}s. Block size = {n30030}")
print(f"Known: v_3(R(30030)) = -1")
print(f"Known: GF(3) nullity(both) = 943")
print(f"Therefore: v_3(det_odd) = v_3(det_even) + 1")
print(f"And: v_3(det_odd) >= 944 (since GF(3) lower bound is 943 + gap_odd >= 1)")
print()

# Run mod 3^k for k=1..8 on both blocks
# k=1 (mod 3) is the GF(3) rank computation
# k=2 (mod 9) should start showing the asymmetry
# Higher k refines until convergence

results = {"odd": {}, "even": {}}

for block_name, block_mat in [("odd", Do30030), ("even", De30030)]:
    print(f"\n--- D_{block_name}(30030) ---")

    prev_total = 0
    for k in range(1, 9):
        print(f"\n  k={k} (mod 3^{k} = {3**k}):")
        total, pivots, elapsed = v3_gauss_mod_pk(block_mat, 3, k,
                                                  label=f"30030_{block_name}")
        results[block_name][k] = {
            "total_v": total,
            "pivots_by_v": dict(sorted(pivots.items())),
            "elapsed": elapsed
        }

        improvement = total - prev_total if k > 1 else total
        print(f"  total_v = {total}")
        print(f"  pivots_by_v = {dict(sorted(pivots.items()))}")
        print(f"  elapsed = {elapsed:.1f}s")
        if k > 1:
            print(f"  improvement over k={k-1}: +{improvement}")

        prev_total = total

        # Check if we can detect convergence:
        # If no pivots have valuation = k (the max), then k is large enough
        # and total_v is exact.
        max_v_pivots = pivots.get(k, 0)
        if max_v_pivots == 0:
            print(f"  *** CONVERGED: no pivots at max valuation k={k}. "
                  f"v_3(det_{block_name}(30030)) = {total} ***")
            break
        else:
            print(f"  {max_v_pivots} pivots at max valuation k={k} — "
                  f"not yet converged")


# ============================================================
# SUMMARY
# ============================================================

print(f"\n{'=' * 70}")
print("SUMMARY")
print("=" * 70)

print(f"\nmod-3^k totals (lower bounds on v_3):")
print(f"  {'k':>2s}  {'mod':>6s}  {'odd':>8s}  {'even':>8s}  {'diff':>8s}")
print(f"  {'-'*2}  {'-'*6}  {'-'*8}  {'-'*8}  {'-'*8}")

for k in range(1, 9):
    if k in results["odd"] and k in results["even"]:
        vo = results["odd"][k]["total_v"]
        ve = results["even"][k]["total_v"]
        diff = ve - vo
        print(f"  {k:>2d}  {3**k:>6d}  {vo:>8d}  {ve:>8d}  {diff:>8d}")

print(f"\n  Required: diff must converge to -1 (= v_3(R(30030)))")

# Check final values
max_k_odd = max(results["odd"].keys())
max_k_even = max(results["even"].keys())
final_odd = results["odd"][max_k_odd]["total_v"]
final_even = results["even"][max_k_even]["total_v"]

print(f"\n  Best estimates:")
print(f"    v_3(det_odd(30030))  >= {final_odd}")
print(f"    v_3(det_even(30030)) >= {final_even}")
print(f"    diff = {final_even - final_odd}")
if final_even - final_odd == -1:
    print(f"    *** CONSISTENT WITH v_3(R) = -1 ***")

# Cross-reference with known table
print(f"\n  Complete 3-adic table:")
print(f"  {'m':>7s}  {'v3(det_e)':>10s}  {'v3(det_o)':>10s}  {'diff':>6s}")
print(f"  {'-'*7}  {'-'*10}  {'-'*10}  {'-'*6}")
known = [(6, 0, 0), (30, 4, 1), (210, 11, 7), (2310, 77, 75)]
for m, ve, vo in known:
    print(f"  {m:>7d}  {ve:>10d}  {vo:>10d}  {ve-vo:>6d}")
print(f"  {30030:>7d}  {final_even:>10d}  {final_odd:>10d}  "
      f"{final_even - final_odd:>6d}")

# Gap analysis
print(f"\n  Gap analysis (v_3 - GF(3) nullity):")
print(f"  {'m':>7s}  {'null_e':>7s}  {'v3_e':>7s}  {'gap_e':>7s}  "
      f"{'null_o':>7s}  {'v3_o':>7s}  {'gap_o':>7s}")
gaps = [
    (210, 7, 11, 7, 7),
    (2310, 75, 77, 75, 75),
    (30030, 943, final_even, 943, final_odd),
]
for m, ne, ve, no, vo in gaps:
    print(f"  {m:>7d}  {ne:>7d}  {ve:>7d}  {ve-ne:>7d}  "
          f"{no:>7d}  {vo:>7d}  {vo-no:>7d}")

# Sequence analysis
print(f"\n  v_3(det_odd) full sequence: 0, 1, 7, 75, {final_odd}")
print(f"  v_3(det_even) full sequence: 0, 4, 11, 77, {final_even}")

# Factor the new values
for label, val in [("v_3(det_odd)", final_odd), ("v_3(det_even)", final_even)]:
    # Simple trial factorization
    n = val
    factors = []
    for p in range(2, int(n**0.5) + 2):
        while n % p == 0:
            factors.append(p)
            n //= p
    if n > 1:
        factors.append(n)
    print(f"  {label} = {val} = {' x '.join(map(str, factors))}")

print(f"\nCompleted {datetime.now().isoformat()}")
print(f"Total runtime: see timestamps above")
