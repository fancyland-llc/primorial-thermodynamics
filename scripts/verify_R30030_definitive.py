"""
verify_R30030_definitive.py — Definitive verification of R(30030).

Strategy:
1. Expose the vacuous-pass bug in Step 2 of verify_R30030_v2.py
2. Fast vectorized numpy det mod p (100x faster than pure Python)
3. Recompute the 2 suspect residues (p=999999877, p=999999443)
4. Also compute at 3 NEW primes not in original set
5. Verify -48317287/3 against ALL results
"""

import numpy as np
import math
import time
from fractions import Fraction
from sympy import isprime

# ============================================================
print("=" * 70)
print("PART 1: EXPOSE THE VACUOUS-PASS BUG")
print("=" * 70)

# From Step 2 of verify_R30030_v2.py:
a_unreduced = -48317254144248150265649657
b_unreduced = 2999997960000205533

# The two suspect primes
p2 = 999999877
p8 = 999999443

# Show that b = 3 * p2 * p8
product = 3 * p2 * p8
print(f"b_unreduced = {b_unreduced}")
print(f"3 * {p2} * {p8} = {product}")
print(f"Match: {b_unreduced == product}")
print()

# Show that a = -48317287 * p2 * p8
k = p2 * p8
a_check = -48317287 * k
print(f"a_unreduced = {a_unreduced}")
print(f"-48317287 * {p2} * {p8} = {a_check}")
print(f"Match: {a_unreduced == a_check}")
print()

# Show the vacuous pass
print(f"b % {p2} = {b_unreduced % p2}  → verify_candidate skips to 'vacuous pass'")
print(f"a % {p2} = {a_unreduced % p2}  → a % p == 0, so it passes vacuously")
print(f"b % {p8} = {b_unreduced % p8}  → verify_candidate skips to 'vacuous pass'")
print(f"a % {p8} = {a_unreduced % p8}  → a % p == 0, so it passes vacuously")
print()
print(">>> The '10/10 pass' in Step 2 is WRONG.")
print(">>> It is really '8 verified passes + 2 vacuous (unchecked) passes'.")
print(">>> The actual residues at p=999999877 and p=999999443 were NEVER checked.")

# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: VERIFY PRIMALITY OF ALL 10 PRIMES")
print("=" * 70)

all_primes = [
    999999937, 999999893, 999999877, 999999613, 999999541,
    999999527, 999999503, 999999491, 999999443, 999999391,
]

for p in all_primes:
    is_p = isprime(p)
    suspect = " ← SUSPECT" if p in [999999877, 999999443] else ""
    print(f"  {p}: {'PRIME' if is_p else '*** NOT PRIME ***'}{suspect}")

# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: FAST VECTORIZED DETERMINANT MOD P")
print("=" * 70)

def build_blocks(m):
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


def det_mod_p_fast(M_int, p):
    """
    Vectorized Gaussian elimination mod p using numpy int64.
    ~100x faster than pure Python for large matrices.
    
    Safe for p < 3*10^9 (products < int64 max ≈ 9.2*10^18).
    """
    n = M_int.shape[0]
    # Convert to int64 mod p. Handle negative entries.
    M = M_int.astype(np.int64) % p
    M = ((M % p) + p) % p  # Ensure non-negative
    
    det_val = 1
    
    for col in range(n):
        # Find pivot: first nonzero in column below diagonal
        col_below = M[col:, col]
        nonzero_idx = np.nonzero(col_below)[0]
        
        if len(nonzero_idx) == 0:
            return 0  # Singular
        
        pivot = nonzero_idx[0] + col
        
        if pivot != col:
            M[[col, pivot]] = M[[pivot, col]].copy()
            det_val = (-det_val) % p
        
        pivot_val = int(M[col, col])
        det_val = (det_val * pivot_val) % p
        inv_pivot = pow(pivot_val, p - 2, p)
        
        if col == n - 1:
            break
        
        # Vectorized elimination of all rows below pivot
        # factors[i] = M[col+1+i, col] * inv_pivot mod p
        below_col = M[col + 1:, col].astype(np.int64)
        
        # Filter: only process rows with nonzero factor
        nonzero_rows = np.nonzero(below_col)[0]
        if len(nonzero_rows) == 0:
            continue
        
        # Compute factors for nonzero rows only
        # Use Python ints for the multiply to avoid int64 overflow edge cases
        factors = np.empty(len(nonzero_rows), dtype=np.int64)
        for idx, row_idx in enumerate(nonzero_rows):
            factors[idx] = (int(below_col[row_idx]) * inv_pivot) % p
        
        # pivot_row: M[col, col+1:]
        pivot_row = M[col, col + 1:].astype(np.int64)
        
        # For each nonzero row: M[row, col+1:] -= factor * pivot_row, mod p
        # Use chunked processing to avoid creating huge temporary arrays
        actual_rows = nonzero_rows + col + 1
        
        # Vectorized: outer product subtraction
        # result[i,j] = M[actual_rows[i], col+1+j] - factors[i] * pivot_row[j]
        # Do this in blocks to manage memory
        remaining = n - col - 1
        
        if remaining > 0:
            # Process all nonzero rows at once using outer product
            # factors: shape (k,), pivot_row: shape (remaining,)
            # outer: shape (k, remaining) — elements up to p^2 < int64 max
            
            # We need to be careful: factors[i] * pivot_row[j] could approach 10^18
            # and M[row,j] is < p, so the subtraction is in [-10^18, p]
            # Then mod p gives [0, p-1]
            
            submatrix = M[actual_rows][:, col + 1:]  # shape (k, remaining)
            
            # Compute factor * pivot_row for each factor
            # Use broadcasting: factors[:, None] * pivot_row[None, :] 
            update = (factors[:, None] * pivot_row[None, :]) % p
            submatrix = (submatrix - update + p) % p  # +p to ensure non-negative before mod
            
            M[actual_rows, col + 1:] = submatrix

        # Zero out the column below pivot
        M[col + 1:, col] = 0
    
    return det_val % p


# Validate on m=2310 first
print("Validating on m=2310 (known R = -108927)...")
t0 = time.time()
D_e_2310, D_o_2310, n_2310 = build_blocks(2310)
print(f"  Matrix size: {n_2310}x{n_2310}")

p_val = 999999937
t1 = time.time()
de_val = det_mod_p_fast(D_e_2310, p_val)
te = time.time() - t1
t2 = time.time()
do_val = det_mod_p_fast(D_o_2310, p_val)
to = time.time() - t2

inv_do = pow(do_val, p_val - 2, p_val)
R_val = (de_val * inv_do) % p_val
if R_val > p_val // 2:
    R_val -= p_val

print(f"  det_e mod p = {de_val}  [{te:.1f}s]")
print(f"  det_o mod p = {do_val}  [{to:.1f}s]")
print(f"  R mod p = {R_val}")
print(f"  Expected:  -108927")
print(f"  VALIDATION: {'PASS' if R_val == -108927 else 'FAIL'}")

# ============================================================
print(f"\n{'=' * 70}")
print("PART 4: BUILD m=30030 MATRICES AND RECOMPUTE")
print("=" * 70)

t0 = time.time()
D_e_30030, D_o_30030, n_30030 = build_blocks(30030)
print(f"m=30030: {n_30030}x{n_30030} blocks, built in {time.time()-t0:.1f}s")

# Candidate fraction
R_candidate = Fraction(-48317287, 3)
print(f"Candidate: R = {R_candidate}")
print(f"  Numerator {R_candidate.numerator} is prime: {isprime(abs(R_candidate.numerator))}")

# Primes to test: the 2 suspect ones + 3 new ones
test_primes = [
    (999999877, "SUSPECT (index 2)"),
    (999999443, "SUSPECT (index 8)"),
    (999999103, "NEW"),
    (999998981, "NEW"),
    (999998857, "NEW"),
]

results = []
for p, label in test_primes:
    print(f"\n  Computing R mod {p}  [{label}]...")
    
    t1 = time.time()
    de = det_mod_p_fast(D_e_30030, p)
    te = time.time() - t1
    
    t2 = time.time()
    do = det_mod_p_fast(D_o_30030, p)
    to = time.time() - t2
    
    print(f"    det_e = {de}, det_o = {do}  [{te+to:.1f}s]")
    
    if do == 0:
        print(f"    det_o = 0 → R undefined mod {p}")
        results.append((p, None, label))
        continue
    
    inv_do = pow(do, p - 2, p)
    r_actual = (de * inv_do) % p
    
    # What does -48317287/3 predict?
    inv3 = pow(3, p - 2, p)
    r_expected = (-48317287 * inv3) % p
    
    match = (r_actual == r_expected)
    print(f"    R mod p = {r_actual}")
    print(f"    -48317287/3 mod p = {r_expected}")
    print(f"    MATCH: {match}")
    
    results.append((p, r_actual, label, r_expected, match))

# ============================================================
print(f"\n{'=' * 70}")
print("PART 5: DEFINITIVE SUMMARY")
print("=" * 70)

# Original 8 consistent residues
original_consistent = [
    (317227550, 999999937),
    (650560833, 999999893),
    (317227442, 999999613),
    (317227418, 999999541),
    (650560589, 999999527),
    (650560573, 999999503),
    (650560565, 999999491),
    (317227368, 999999391),
]

# Original 2 suspect residues
original_suspect = [
    (633263010, 999999877),
    (630757884, 999999443),
]

print("\nOriginal 8 consistent residues vs -48317287/3:")
all_orig_pass = True
for r, p in original_consistent:
    inv3 = pow(3, p - 2, p)
    expected = (-48317287 * inv3) % p
    match = (expected == r)
    if not match:
        all_orig_pass = False
        print(f"  p={p}: FAIL (expected {expected}, got {r})")
print(f"  All 8 pass: {all_orig_pass}")

print(f"\nRecomputed suspect primes:")
for entry in results:
    p = entry[0]
    if entry[1] is None:
        print(f"  p={p}: det_o=0 (R undefined)")
    else:
        _, r_actual, label, r_expected, match = entry
        print(f"  p={p} [{label}]: {'MATCH' if match else 'MISMATCH'} (actual={r_actual}, expected={r_expected})")

print(f"\nNew primes:")
new_results = [e for e in results if "NEW" in e[2]]
for entry in new_results:
    if entry[1] is None:
        print(f"  p={entry[0]}: det_o=0")
    else:
        _, r_actual, label, r_expected, match = entry
        print(f"  p={entry[0]}: {'MATCH' if match else 'MISMATCH'}")

# Count total
total_pass = 8  # Original consistent
total_tested = 8
for entry in results:
    if entry[1] is not None:
        total_tested += 1
        if entry[4]:  # match
            total_pass += 1

print(f"\n{'=' * 70}")
print(f"TOTAL: {total_pass}/{total_tested} independent primes confirm R(30030) = -48317287/3")
print(f"{'=' * 70}")

if total_pass == total_tested:
    print("\n*** R(30030) = -48317287/3 DEFINITIVELY CONFIRMED ***")
    print(f"*** The 2 original suspect residues were COMPUTATIONAL ERRORS ***")
    
    # Full analysis
    print(f"\n  R(30030) = -48317287/3")
    print(f"  |numerator| = 48317287 is PRIME")
    print(f"  denominator = 3 = p_min")
    
    # Polynomial hierarchy
    neg2R = Fraction(-2) * R_candidate
    a4 = 217855
    residual = neg2R - a4 * 169 + 1
    b4_candidate = residual / 13
    print(f"\n  Polynomial hierarchy at k=4:")
    print(f"  -2R(30030) = {neg2R}")
    print(f"  b_4 = {b4_candidate}")
    print(f"  b_4 is integer: {b4_candidate.denominator == 1}")
    print(f"  b_4 denominator: {b4_candidate.denominator}")
    print(f"\n  *** HIERARCHY TERMINATES: b_4 = {b4_candidate} is NOT an integer ***")
    print(f"  *** The polynomial hierarchy holds for k=0,1,2,3 (m=6,30,210,2310) ***")
    print(f"  *** and BREAKS at k=4 (m=30030) with denominator 3 = p_min ***")
    
    # The sequence
    print(f"\n  Complete R-sequence:")
    data = [
        (6, Fraction(-1), "integer"),
        (30, Fraction(-27), "integer"),
        (210, Fraction(-1053), "integer"),
        (2310, Fraction(-108927), "integer"),
        (30030, R_candidate, f"rational, denom={R_candidate.denominator}"),
    ]
    for m, R, typ in data:
        print(f"    R({m:>5d}) = {str(R):>20s}  ({typ})")
    
    # Ratios
    print(f"\n  Consecutive ratios R(m_{k+1})/R(m_k):")
    vals = [d[1] for d in data]
    for i in range(1, len(vals)):
        ratio = vals[i] / vals[i-1]
        print(f"    R({data[i][0]})/R({data[i-1][0]}) = {ratio} = {float(ratio):.6f}")
    
elif total_pass < total_tested:
    mismatch_count = total_tested - total_pass
    print(f"\n  {mismatch_count} primes failed. R may not be -48317287/3.")
    print(f"  Need further investigation.")
