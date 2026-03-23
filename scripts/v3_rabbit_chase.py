"""
v3_rabbit_chase.py — Chase the 3-adic valuation rabbit.

Computes the exact 3-adic valuations v₃(det_even) and v₃(det_odd) for
m = 6, 30, 210, 2310 from exact integer determinants (SymPy Bareiss),
then computes the GF(3) nullity of D_odd(30030) and D_even(30030)
to test the prediction v₃(det_odd(30030)) = 210 = m₃ (third primorial).

PART I:   Exact v₃ table from saved determinants (instant)
PART II:  GF(3) nullity of 2880×2880 blocks for m=30030
PART III: Higher-power lifting (mod 9, mod 27, mod 81, ...)

Author: Antonio P. Matos / Fancyland LLC
Date: March 23, 2026
"""

import numpy as np
import math
import time
import sys
import os
from fractions import Fraction

try:
    from sympy import Matrix, factorint, isprime
except ImportError:
    print("ERROR: sympy required. Install with: pip install sympy")
    sys.exit(1)

# ============================================================
# Infrastructure
# ============================================================

def v_p(n, p):
    """p-adic valuation of integer n. Returns (valuation, unit_part)."""
    if n == 0:
        return float('inf'), 0
    n = int(n)
    sign = 1 if n > 0 else -1
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v, sign * n

def build_blocks_sympy(m):
    """Build D_even and D_odd as SymPy matrices (exact integer determinants)."""
    coprimes = sorted([r for r in range(1, m) if math.gcd(r, m) == 1])
    pairs = [r for r in coprimes if r < m / 2]
    n = len(pairs)
    De = [[0] * n for _ in range(n)]
    Do = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            ri, rj = pairs[i], pairs[j]
            diff = (rj - ri) % m
            d_diff = min(diff, m - diff)
            summ = (ri + rj) % m
            d_sum = min(summ, m - summ)
            De[i][j] = d_diff + d_sum
            Do[i][j] = d_diff - d_sum
    return Matrix(De), Matrix(Do), n

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

def rank_mod_p(M_int, p):
    """
    Compute the rank of integer matrix M_int over GF(p) using
    Gaussian elimination. Returns (rank, num_zero_pivots).
    
    This is the same as det_mod_p_fast but we count zero pivots
    instead of accumulating the determinant product.
    """
    n = M_int.shape[0]
    M = M_int.astype(np.int64) % p
    M = ((M % p) + p) % p
    
    rank = 0
    pivot_col = 0
    pivot_row_map = []  # track which rows have pivots
    
    row_current = 0
    for col in range(n):
        # Find pivot in this column from row_current downward
        col_below = M[row_current:, col]
        nonzero_idx = np.nonzero(col_below)[0]
        
        if len(nonzero_idx) == 0:
            continue  # skip this column — zero pivot
        
        pivot = nonzero_idx[0] + row_current
        if pivot != row_current:
            M[[row_current, pivot]] = M[[pivot, row_current]].copy()
        
        pivot_val = int(M[row_current, col])
        inv_pivot = pow(pivot_val, p - 2, p)
        
        # Eliminate below
        below_col = M[row_current + 1:, col].astype(np.int64)
        nonzero_rows = np.nonzero(below_col)[0]
        
        if len(nonzero_rows) > 0:
            factors = np.empty(len(nonzero_rows), dtype=np.int64)
            for idx, row_idx in enumerate(nonzero_rows):
                factors[idx] = (int(below_col[row_idx]) * inv_pivot) % p
            
            pivot_row = M[row_current, col + 1:].astype(np.int64)
            actual_rows = nonzero_rows + row_current + 1
            
            if n - col - 1 > 0:
                submatrix = M[actual_rows][:, col + 1:]
                update = (factors[:, None] * pivot_row[None, :]) % p
                submatrix = (submatrix - update + p) % p
                M[actual_rows, col + 1:] = submatrix
        
        M[row_current + 1:, col] = 0
        rank += 1
        row_current += 1
        
        if row_current >= n:
            break
    
    return rank


def gauss_elim_mod_pk(D_int, p, k):
    """
    Gaussian elimination mod p^k with p-adic pivot selection.
    Returns the total p-adic valuation contribution from all pivots.
    
    For each column, find the entry with minimum p-adic valuation,
    use it as pivot, and track the total valuation.
    
    Uses Python arbitrary-precision integers (no numpy overflow).
    """
    pk = p ** k
    n = D_int.shape[0]
    
    # Convert to Python int matrix for exact mod arithmetic
    M = [[int(D_int[i, j]) % pk for j in range(n)] for i in range(n)]
    
    total_v = 0
    
    for col in range(n):
        if col % 200 == 0:
            print(f"    col {col}/{n}, v_so_far={total_v}", flush=True)
        
        # Find pivot: row with minimum p-adic valuation in column col
        best_v = k + 1  # sentinel: larger than any possible valuation mod p^k
        best_row = -1
        
        for row in range(col, n):
            entry = M[row][col] % pk
            if entry == 0:
                continue
            v = 0
            x = entry
            while x % p == 0 and v < k:
                x //= p
                v += 1
            if v < best_v:
                best_v = v
                best_row = row
                if v == 0:
                    break  # can't do better than 0
        
        if best_row == -1:
            # Entire column is 0 mod p^k — valuation >= k from this column
            total_v += k
            continue
        
        total_v += best_v
        
        # Swap rows
        if best_row != col:
            M[col], M[best_row] = M[best_row], M[col]
        
        # Pivot value, divided by p^best_v
        pivot_full = M[col][col]
        pivot_reduced = pivot_full // (p ** best_v)
        
        # We need pivot_reduced to be invertible mod p^(k - best_v)
        remaining_k = k - best_v
        if remaining_k <= 0:
            continue
        
        pk_remaining = p ** remaining_k
        
        # Check if pivot_reduced is invertible mod p
        if pivot_reduced % p == 0:
            # Not invertible — skip standard elimination, more complex handling needed
            continue
        
        try:
            pivot_inv = pow(pivot_reduced % pk_remaining, -1, pk_remaining)
        except ValueError:
            continue
        
        # Eliminate below
        for row in range(col + 1, n):
            entry = M[row][col]
            if entry % pk == 0:
                continue
            
            # entry = p^best_v * (entry // p^best_v) + remainder
            entry_div = entry // (p ** best_v)
            factor = (entry_div * pivot_inv) % pk_remaining
            
            for j in range(col + 1, n):
                M[row][j] = (M[row][j] - factor * M[col][j]) % pk
            M[row][col] = 0
    
    return total_v


# ============================================================
print("=" * 70)
print("PART I: EXACT 3-ADIC VALUATIONS FROM SAVED DETERMINANTS")
print("=" * 70)
# ============================================================

# m=6: compute directly
print("\n  m=6:")
De6, Do6, n6 = build_blocks_sympy(6)
det_e6 = int(De6.det())
det_o6 = int(Do6.det())
v3_e6, _ = v_p(det_e6, 3)
v3_o6, _ = v_p(det_o6, 3)
print(f"    det_even = {det_e6}, v_3 = {v3_e6}")
print(f"    det_odd  = {det_o6}, v_3 = {v3_o6}")
print(f"    diff = {v3_e6 - v3_o6}")

# m=30: compute directly
print("\n  m=30:")
De30, Do30, n30 = build_blocks_sympy(30)
det_e30 = int(De30.det())
det_o30 = int(Do30.det())
v3_e30, _ = v_p(det_e30, 3)
v3_o30, _ = v_p(det_o30, 3)
print(f"    det_even = {det_e30}, v_3 = {v3_e30}")
print(f"    det_odd  = {det_o30}, v_3 = {v3_o30}")
print(f"    diff = {v3_e30 - v3_o30}")

# m=210: compute directly (still fast — 24x24 Bareiss)
print("\n  m=210:")
t0 = time.time()
De210, Do210, n210 = build_blocks_sympy(210)
det_e210 = int(De210.det())
det_o210 = int(Do210.det())
v3_e210, unit_e210 = v_p(det_e210, 3)
v3_o210, unit_o210 = v_p(det_o210, 3)
print(f"    det_even = {det_e210}")
print(f"    det_odd  = {det_o210}")
print(f"    v_3(det_even) = {v3_e210}, v_3(det_odd) = {v3_o210}")
print(f"    diff = {v3_e210 - v3_o210}")
print(f"    [{time.time()-t0:.1f}s]")

# m=2310: load from saved file
print("\n  m=2310:")
det_file = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', '..', '..', 'backend', 'data', 'gbd_antigen_solver',
    'bvp_7_fusion_confinement', 'det_2310_exact.txt'
)
# Also try relative to workspace root
det_file_alt = r"C:\Users\apmat\AppData\Local\dev\prompt_studio\backend\data\gbd_antigen_solver\bvp_7_fusion_confinement\det_2310_exact.txt"

det_e2310 = None
det_o2310 = None

for fpath in [det_file, det_file_alt]:
    if os.path.exists(fpath):
        print(f"    Loading from {fpath}")
        with open(fpath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('det_even'):
                    # Parse: det_even = <number>
                    det_e2310 = int(line.split('=')[1].strip())
                elif line.startswith('det_odd'):
                    det_o2310 = int(line.split('=')[1].strip())
        break

if det_e2310 is None or det_o2310 is None:
    print("    Saved file not found — computing from scratch (this takes ~80s)...")
    t0 = time.time()
    De2310, Do2310, n2310 = build_blocks_sympy(2310)
    det_e2310 = int(De2310.det())
    det_o2310 = int(Do2310.det())
    print(f"    Computed in {time.time()-t0:.1f}s")

v3_e2310, unit_e2310 = v_p(det_e2310, 3)
v3_o2310, unit_o2310 = v_p(det_o2310, 3)
print(f"    v_3(det_even) = {v3_e2310}")
print(f"    v_3(det_odd)  = {v3_o2310}")
print(f"    diff = {v3_e2310 - v3_o2310}")
R2310 = Fraction(det_e2310, det_o2310)
print(f"    R(2310) = {R2310}")
v3_R2310, _ = v_p(abs(R2310.numerator), 3)
v3_R2310_den, _ = v_p(R2310.denominator, 3)
print(f"    v_3(R) = {v3_R2310 - v3_R2310_den} (check: should be 2)")

# Summary table
print(f"\n  {'='*60}")
print(f"  COMPLETE 3-ADIC VALUATION TABLE")
print(f"  {'='*60}")
print(f"  {'m':>7s} | {'v3(det_e)':>10s} | {'v3(det_o)':>10s} | {'diff':>6s} | {'v3(R)':>6s}")
print(f"  {'-'*7}-+-{'-'*10}-+-{'-'*10}-+-{'-'*6}-+-{'-'*6}")

table_data = [
    (6, v3_e6, v3_o6),
    (30, v3_e30, v3_o30),
    (210, v3_e210, v3_o210),
    (2310, v3_e2310, v3_o2310),
]

for m_val, ve, vo in table_data:
    d = ve - vo
    print(f"  {m_val:>7d} | {ve:>10d} | {vo:>10d} | {d:>6d} | {d:>6d}")

print(f"  {30030:>7d} | {'?':>10s} | {'?':>10s} | {-1:>6d} | {-1:>6d}")

# Check the Claude/Sonnet table: v3(det_odd) should be 0, 1, 7, 75
print(f"\n  v_3(det_odd) sequence: {v3_o6}, {v3_o30}, {v3_o210}, {v3_o2310}, ???")
print(f"  v_3(det_even) sequence: {v3_e6}, {v3_e30}, {v3_e210}, {v3_e2310}, ???")
print(f"  Differences: {v3_e6-v3_o6}, {v3_e30-v3_o30}, {v3_e210-v3_o210}, {v3_e2310-v3_o2310}, -1")

# ============================================================
print(f"\n{'='*70}")
print("PART II: GF(3) NULLITY OF D_odd(30030) AND D_even(30030)")
print("=" * 70)
# ============================================================

print("\n  Building 2880×2880 blocks for m=30030...")
t0 = time.time()
De30030, Do30030, n30030 = build_blocks_numpy(30030)
print(f"  Built in {time.time()-t0:.1f}s. Block size = {n30030}")

# GF(3) rank via Gaussian elimination mod 3
print("\n  Computing rank of D_odd(30030) mod 3...")
t0 = time.time()
rank_odd_3 = rank_mod_p(Do30030, 3)
nullity_odd_3 = n30030 - rank_odd_3
t_odd = time.time() - t0
print(f"  Rank(D_odd mod 3) = {rank_odd_3}")
print(f"  Nullity(D_odd mod 3) = {nullity_odd_3}")
print(f"  => v_3(det_odd(30030)) >= {nullity_odd_3}")
print(f"  [{t_odd:.1f}s]")

print("\n  Computing rank of D_even(30030) mod 3...")
t0 = time.time()
rank_even_3 = rank_mod_p(De30030, 3)
nullity_even_3 = n30030 - rank_even_3
t_even = time.time() - t0
print(f"  Rank(D_even mod 3) = {rank_even_3}")
print(f"  Nullity(D_even mod 3) = {nullity_even_3}")
print(f"  => v_3(det_even(30030)) >= {nullity_even_3}")
print(f"  [{t_even:.1f}s]")

print(f"\n  GF(3) nullity difference: {nullity_even_3} - {nullity_odd_3} = {nullity_even_3 - nullity_odd_3}")
print(f"  Known: v_3(R(30030)) = v_3(det_even) - v_3(det_odd) = -1")

# Test predictions
print(f"\n  PREDICTIONS:")
print(f"  v_3(det_odd(30030)) = 210 = m_3 (third primorial)?  nullity = {nullity_odd_3}, match = {nullity_odd_3 == 210}")
print(f"  v_3(det_even(30030)) = 209 = 11*19?  nullity = {nullity_even_3}, match = {nullity_even_3 == 209}")

# Also check: is nullity_odd related to any primorial?
for name, val in [("m_1=2", 2), ("m_2=6", 6), ("m_3=210", 210), 
                   ("phi(210)=48", 48), ("phi(2310)=480", 480),
                   ("phi(30030)/2=2880", 2880)]:
    if nullity_odd_3 == val:
        print(f"  *** nullity_odd = {val} = {name} ***")
for name, val in [("m_1=2", 2), ("m_2=6", 6), ("m_3=210", 210),
                   ("phi(210)=48", 48), ("phi(2310)=480", 480)]:
    if nullity_even_3 == val:
        print(f"  *** nullity_even = {val} = {name} ***")

# ============================================================
print(f"\n{'='*70}")
print("PART III: HIGHER-POWER LIFTING (mod 9, mod 27, mod 81)")
print("=" * 70)
# ============================================================

# First, validate the method on known cases (m=210, m=2310)
print("\n  Validation: checking GF(3) nullity against known valuations...")

# m=210
print("\n  m=210 (known: v_3(det_odd)=7, v_3(det_even)=11):")
De210np, Do210np, _ = build_blocks_numpy(210)
r_o210 = rank_mod_p(Do210np, 3)
r_e210 = rank_mod_p(De210np, 3)
null_o210 = Do210np.shape[0] - r_o210
null_e210 = De210np.shape[0] - r_e210
print(f"    GF(3) nullity(D_odd)  = {null_o210} (v_3 = {v3_o210})")
print(f"    GF(3) nullity(D_even) = {null_e210} (v_3 = {v3_e210})")
print(f"    nullity == v_3 for odd?  {null_o210 == v3_o210}")
print(f"    nullity == v_3 for even? {null_e210 == v3_e210}")

# m=2310 
print("\n  m=2310 (known: v_3(det_odd)={}, v_3(det_even)={}):".format(v3_o2310, v3_e2310))
De2310np, Do2310np, _ = build_blocks_numpy(2310)
t0 = time.time()
r_o2310 = rank_mod_p(Do2310np, 3)
r_e2310 = rank_mod_p(De2310np, 3)
null_o2310 = Do2310np.shape[0] - r_o2310
null_e2310 = De2310np.shape[0] - r_e2310
print(f"    GF(3) nullity(D_odd)  = {null_o2310} (v_3 = {v3_o2310})")
print(f"    GF(3) nullity(D_even) = {null_e2310} (v_3 = {v3_e2310})")
print(f"    nullity == v_3 for odd?  {null_o2310 == v3_o2310}")
print(f"    nullity == v_3 for even? {null_e2310 == v3_e2310}")
print(f"    [{time.time()-t0:.1f}s]")

# If GF(3) nullity != v_3, higher powers contribute.
# The gap = v_3 - nullity tells us how many additional factors of 3
# come from higher-order structure (Jordan blocks over GF(3)).

if null_o2310 != v3_o2310 or null_e2310 != v3_e2310:
    print("\n    GF(3) nullity does NOT equal v_3 — higher powers contribute!")
    print(f"    Gap for odd:  {v3_o2310 - null_o2310}")
    print(f"    Gap for even: {v3_e2310 - null_e2310}")
    print("    The GF(3) nullity is only a LOWER BOUND on v_3.")
    print("    Need mod 9, mod 27, ... lifting to get exact v_3 for m=30030.")
else:
    print("\n    GF(3) nullity EQUALS v_3 at m=210 and m=2310!")
    print("    If this pattern holds, GF(3) nullity directly gives v_3(det(30030)).")

# ============================================================
# PART IV: If GF(3) nullity is only a lower bound, try mod 9 lifting
# ============================================================

# For m=30030, try Gaussian elimination mod 9 on the odd block
# to get a tighter bound. Each entry is in [0,8].
# With 2880×2880 and Python arbitrary precision, this will be slower.

print(f"\n{'='*70}")
print("PART IV: MOD-9 AND MOD-27 ANALYSIS FOR m=30030")
print("=" * 70)

print("\n  Computing rank of D_odd(30030) mod 9 (p=3, k=2)...")
print("  (This uses p-adic pivot selection in pure Python — may be slow)")
t0 = time.time()

# For mod p^k, we use a simpler O(n^3) approach: 
# Gaussian elimination where at each step we track the minimum
# p-adic valuation of the pivot. This is slower but correct.

# Actually for mod 9 we can still use the vectorized approach 
# if we're careful about invertibility.
# In Z/9Z, an element is invertible iff gcd(element, 9) = 1.
# Elements {1,2,4,5,7,8} are invertible; {0,3,6} are not.

def rank_and_v_mod_pk_fast(M_int, p, k, max_cols=None):
    """
    Gaussian elimination over Z/p^k Z with p-adic pivot selection.
    Returns total p-adic valuation from pivots (lower bound on v_p(det)).
    
    Vectorized where possible, falls back to scalar for pivot selection.
    """
    pk = p ** k
    n = M_int.shape[0]
    if max_cols is None:
        max_cols = n
    
    # Use Python objects for exact arithmetic
    M = (M_int.astype(object) % pk + pk) % pk
    
    total_v = 0
    row_current = 0
    
    for col in range(min(n, max_cols)):
        if col % 500 == 0 and col > 0:
            elapsed = time.time() - t0
            rate = col / elapsed
            eta = (min(n, max_cols) - col) / rate if rate > 0 else 0
            print(f"    col {col}/{min(n,max_cols)}, v_so_far={total_v}, "
                  f"elapsed={elapsed:.0f}s, ETA={eta:.0f}s", flush=True)
        
        # Find pivot with minimum p-adic valuation
        best_v = k + 1
        best_row = -1
        
        for row in range(row_current, n):
            entry = int(M[row, col]) % pk
            if entry == 0:
                continue
            v = 0
            x = entry
            while x % p == 0 and v < k:
                x //= p
                v += 1
            if v < best_v:
                best_v = v
                best_row = row
                if v == 0:
                    break
        
        if best_row == -1:
            total_v += k
            continue
        
        total_v += best_v
        
        if best_row != row_current:
            M[[row_current, best_row]] = M[[best_row, row_current]].copy()
        
        # pivot_val / p^best_v must be invertible mod p
        pivot_val = int(M[row_current, col])
        pivot_reduced = pivot_val // (p ** best_v)
        
        if pivot_reduced % p == 0:
            # Can't eliminate cleanly — higher p-adic structure
            row_current += 1
            continue
        
        remaining_pk = p ** (k - best_v)
        if remaining_pk <= 1:
            row_current += 1
            continue
        
        try:
            pivot_inv = pow(pivot_reduced % remaining_pk, -1, remaining_pk)
        except (ValueError, ZeroDivisionError):
            row_current += 1
            continue
        
        # Eliminate below
        for row in range(row_current + 1, n):
            entry = int(M[row, col])
            if entry % pk == 0:
                continue
            entry_div = entry // (p ** best_v)
            factor = (entry_div * pivot_inv) % remaining_pk
            
            for j in range(col + 1, n):
                M[row, j] = (int(M[row, j]) - factor * int(M[row_current, j])) % pk
            M[row, col] = 0
        
        row_current += 1
    
    return total_v

# Only run mod-9 if it's going to be informative
# First check: did GF(3) nullity match known values?
if null_o2310 == v3_o2310:
    print("\n  GF(3) nullity was exact at m=2310 — skipping mod-9 lifting.")
    print(f"  Best estimate: v_3(det_odd(30030)) = {nullity_odd_3}")
    print(f"  Best estimate: v_3(det_even(30030)) = {nullity_even_3}")
else:
    print("\n  GF(3) nullity was NOT exact — running mod-9 on D_odd(30030)...")
    t0 = time.time()
    v_mod9 = rank_and_v_mod_pk_fast(Do30030, 3, 2)
    print(f"  v_3 from mod-9 elimination >= {v_mod9}")
    print(f"  [{time.time()-t0:.1f}s]")
    
    # And mod 27
    if v_mod9 < 300:
        print("\n  Running mod-27 on D_odd(30030)...")
        t0 = time.time()
        v_mod27 = rank_and_v_mod_pk_fast(Do30030, 3, 3)
        print(f"  v_3 from mod-27 elimination >= {v_mod27}")
        print(f"  [{time.time()-t0:.1f}s]")

# ============================================================
print(f"\n{'='*70}")
print("PART V: FULL p-ADIC VALUATION STRUCTURE OF R(30030)")
print("=" * 70)
# ============================================================

# R(30030) = -48317287/3. Numerator is prime.
# Therefore v_p(R) = 0 for all p except:
#   v_3(R) = -1 (from denominator)
# This is trivially read from the fraction.

R30030_num = 48317287
R30030_den = 3

print(f"\n  R(30030) = -{R30030_num}/{R30030_den}")
print(f"  {R30030_num} is prime: {isprime(R30030_num)}")
print(f"\n  p-adic valuations of R(30030):")
for p in [2, 3, 5, 7, 11, 13]:
    v_num, _ = v_p(R30030_num, p)
    v_den, _ = v_p(R30030_den, p)
    v_R = v_num - v_den
    print(f"    v_{p}(R(30030)) = {v_R}")

print(f"\n  All primes except 3 have v_p = 0 (because numerator is prime).")
print(f"  The ONLY structural asymmetry is v_3 = -1.")

# Cross-check: v_p structure of ALL ratios
print(f"\n  Complete v_p structure across all primorials:")
print(f"  {'m':>7s} | {'v_2':>4s} | {'v_3':>4s} | {'v_5':>4s} | {'v_7':>4s} | {'v_11':>4s} | {'v_13':>4s}")
print(f"  {'-'*7}-+-{'-'*4}-+-{'-'*4}-+-{'-'*4}-+-{'-'*4}-+-{'-'*4}-+-{'-'*4}")

R_values = {
    6: Fraction(-1),
    30: Fraction(-27),
    210: Fraction(-1053),
    2310: Fraction(-108927),
    30030: Fraction(-48317287, 3),
}

for m_val, R in R_values.items():
    row = f"  {m_val:>7d} |"
    for p in [2, 3, 5, 7, 11, 13]:
        if R == 0:
            v = 'inf'
        else:
            vn, _ = v_p(abs(R.numerator), p)
            vd, _ = v_p(R.denominator, p)
            v = vn - vd
        row += f" {v:>4} |"
    print(row)

# ============================================================
print(f"\n{'='*70}")
print("SUMMARY")
print("=" * 70)
# ============================================================

print(f"""
  CONFIRMED 3-ADIC VALUATION TABLE:
  
  m       v_3(det_e)  v_3(det_o)  diff = v_3(R)
  ------  ----------  ----------  -------------
  6       {v3_e6:>10d}  {v3_o6:>10d}  {v3_e6-v3_o6:>13d}
  30      {v3_e30:>10d}  {v3_o30:>10d}  {v3_e30-v3_o30:>13d}
  210     {v3_e210:>10d}  {v3_o210:>10d}  {v3_e210-v3_o210:>13d}
  2310    {v3_e2310:>10d}  {v3_o2310:>10d}  {v3_e2310-v3_o2310:>13d}
  30030   {nullity_even_3:>10d}* {nullity_odd_3:>10d}*          -1
  
  * = GF(3) nullity (lower bound; exact if pattern holds from m=210,2310)
  
  v_3(det_odd) sequence: {v3_o6}, {v3_o30}, {v3_o210}, {v3_o2310}, {nullity_odd_3}*
  v_3(det_even) sequence: {v3_e6}, {v3_e30}, {v3_e210}, {v3_e2310}, {nullity_even_3}*
  
  PREDICTION: v_3(det_odd(30030)) = 210 = m_3?  {nullity_odd_3 == 210}
  PREDICTION: v_3(det_even(30030)) = 209 = 11*19?  {nullity_even_3 == 209}
""")
