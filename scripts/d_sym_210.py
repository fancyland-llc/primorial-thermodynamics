"""
d_sym_210.py — Full computation of D_sym for m = 210 = 2*3*5*7
================================================================

This is the next primorial after m=30. We compute:
  1. The 48 coprime residues of Z/210Z
  2. The 48x48 symmetric distance matrix D_sym
  3. The palindromic decomposition under x <-> 210-x
  4. The 24x24 even (D_even) and odd (D_odd) blocks
  5. Cross-coupling (should be zero)
  6. Eigenvalues of both blocks
  7. Determinants, traces
  8. Whether chi_7 lives in the odd block (PREDICTED: yes, since 7 ≡ 3 mod 4)
  9. Whether chi_5 moves to even block (PREDICTED: yes, since 5 ≡ 1 mod 4)
  10. Gauss sum verification: |g(chi_7)| = sqrt(7)
  11. Determinant ratio — does it involve p_mid = 5?
  12. Prefactor: C = phi(210)/|g(chi_7)| = 48/sqrt(7)

Requires: numpy, sympy
"""

import numpy as np
from math import gcd
from sympy import factorint, isprime, sqrt as sym_sqrt, Rational, Matrix, simplify
from sympy import legendre_symbol, exp as sym_exp, pi as sym_pi, I as sym_I
from sympy import Abs as sym_Abs

m = 210

# ============================================================
# Stage 1: Coprime residues
# ============================================================
print("=" * 80)
print(f"STAGE 1: Coprime residues of Z/{m}Z")
print("=" * 80)

coprimes = sorted([x for x in range(1, m) if gcd(x, m) == 1])
phi_m = len(coprimes)
print(f"phi({m}) = {phi_m}")
print(f"Coprime residues ({phi_m} total): {coprimes[:10]}...{coprimes[-10:]}")

# Verify phi(210) = 48
assert phi_m == 48, f"Expected phi(210)=48, got {phi_m}"
print(f"  PASS: phi(210) = 48")

# ============================================================
# Stage 2: Symmetric distance matrix D_sym (48x48)
# ============================================================
print(f"\n{'=' * 80}")
print(f"STAGE 2: Symmetric distance matrix D_sym ({phi_m}x{phi_m})")
print("=" * 80)

def cyclic_distance(a, b, n):
    """Minimum distance on Z/nZ."""
    d = abs(a - b) % n
    return min(d, n - d)

D_sym = np.zeros((phi_m, phi_m), dtype=int)
for i in range(phi_m):
    for j in range(phi_m):
        D_sym[i, j] = cyclic_distance(coprimes[i], coprimes[j], m)

print(f"Matrix shape: {D_sym.shape}")
print(f"Symmetric: {np.allclose(D_sym, D_sym.T)}")
print(f"Diagonal all zero: {np.all(np.diag(D_sym) == 0)}")
print(f"Min off-diag: {D_sym[np.triu_indices(phi_m, 1)].min()}")
print(f"Max entry: {D_sym.max()}")

# ============================================================
# Stage 3: Palindromic decomposition x <-> 210 - x
# ============================================================
print(f"\n{'=' * 80}")
print(f"STAGE 3: Palindromic decomposition under x <-> {m} - x")
print("=" * 80)

# The palindromic involution maps x -> m - x
# For each coprime x, m-x is also coprime (since gcd(m-x, m) = gcd(x, m))
# The fixed points would be x = m/2, but m/2 = 105 and gcd(105, 210) = 105 ≠ 1
# So there are no fixed points — all 48 residues pair up into 24 orbits

palindrome_map = {}
for x in coprimes:
    partner = m - x
    assert partner in coprimes, f"{partner} not coprime to {m}!"
    palindrome_map[x] = partner

# Find the 24 pairs (x, m-x) with x < m-x
pairs = []
for x in coprimes:
    if x < m - x:
        pairs.append((x, m - x))

print(f"Number of palindromic pairs: {len(pairs)}")
assert len(pairs) == 24, f"Expected 24 pairs, got {len(pairs)}"
print(f"  PASS: 24 pairs (no fixed points)")
print(f"First 5 pairs: {pairs[:5]}")
print(f"Last 5 pairs:  {pairs[-5:]}")

# Build the transformation matrix P
# Order: first 24 entries are (x + (m-x))/sqrt(2), next 24 are (x - (m-x))/sqrt(2)
# i.e., even combinations first, odd combinations second

idx = {x: i for i, x in enumerate(coprimes)}

P = np.zeros((phi_m, phi_m))
for k, (x, mx) in enumerate(pairs):
    i = idx[x]
    j = idx[mx]
    # Even: (e_i + e_j) / sqrt(2)
    P[i, k] = 1 / np.sqrt(2)
    P[j, k] = 1 / np.sqrt(2)
    # Odd: (e_i - e_j) / sqrt(2)
    P[i, k + 24] = 1 / np.sqrt(2)
    P[j, k + 24] = -1 / np.sqrt(2)

# Verify P is orthogonal
assert np.allclose(P @ P.T, np.eye(phi_m), atol=1e-12), "P not orthogonal!"
print(f"  PASS: P is orthogonal")

# Transform D_sym
D_transformed = P.T @ D_sym @ P

D_even = D_transformed[:24, :24]
D_odd = D_transformed[24:, 24:]
D_cross_1 = D_transformed[:24, 24:]
D_cross_2 = D_transformed[24:, :24]

cross_norm = np.max(np.abs(D_cross_1))
print(f"D_even shape: {D_even.shape}")
print(f"D_odd shape:  {D_odd.shape}")
print(f"Cross-coupling max |entry|: {cross_norm:.2e}")
assert cross_norm < 1e-10, f"Cross-coupling too large: {cross_norm}"
print(f"  PASS: Exact block decomposition (coupling < 1e-10)")

# ============================================================
# Stage 4: Block properties
# ============================================================
print(f"\n{'=' * 80}")
print(f"STAGE 4: Block properties")
print("=" * 80)

trace_even = np.trace(D_even)
trace_odd = np.trace(D_odd)
det_even = np.linalg.det(D_even)
det_odd = np.linalg.det(D_odd)

print(f"trace(D_even) = {trace_even:.6f}")
print(f"trace(D_odd)  = {trace_odd:.6f}")
print(f"det(D_even)   = {det_even:.6e}")
print(f"det(D_odd)    = {det_odd:.6e}")

if abs(det_odd) > 1e-6:
    ratio = det_even / det_odd
    print(f"det(D_even)/det(D_odd) = {ratio:.6f}")
else:
    print(f"det(D_odd) is near zero — ratio undefined")

# Eigenvalues
eig_even = np.sort(np.linalg.eigvalsh(D_even))[::-1]
eig_odd = np.sort(np.linalg.eigvalsh(D_odd))[::-1]

print(f"\nEigenvalues of D_even (top 10): {eig_even[:10].round(4)}")
print(f"Eigenvalues of D_even (bot 5):  {eig_even[-5:].round(4)}")
print(f"\nEigenvalues of D_odd (top 10):  {eig_odd[:10].round(4)}")
print(f"Eigenvalues of D_odd (bot 5):   {eig_odd[-5:].round(4)}")

print(f"\nNumber of zero eigenvalues in D_even: {np.sum(np.abs(eig_even) < 1e-8)}")
print(f"Number of zero eigenvalues in D_odd:  {np.sum(np.abs(eig_odd) < 1e-8)}")

# ============================================================
# Stage 5: Legendre character separation
# ============================================================
print(f"\n{'=' * 80}")
print(f"STAGE 5: Legendre character separation")
print("=" * 80)

primes_of_m = [2, 3, 5, 7]

for p in primes_of_m[1:]:  # Skip 2 (Legendre symbol mod 2 is trivial)
    # chi_p(-1) determines even/odd
    chi_neg1 = legendre_symbol(-1, p)
    parity = "EVEN" if chi_neg1 == 1 else "ODD"
    mod4 = p % 4
    print(f"  chi_{p}(-1) = {chi_neg1}  ({p} ≡ {mod4} mod 4)  => palindromically {parity}")

# Compute chi_p on coprime residues
print(f"\n  Character values on coprime residues:")
for p in [3, 5, 7]:
    vals = [legendre_symbol(x % p, p) if x % p != 0 else 0 for x in coprimes]
    # Check: does chi_p(x) = chi_p(m-x)?  (even)
    # Or:    chi_p(x) = -chi_p(m-x)?  (odd)
    even_count = 0
    odd_count = 0
    zero_count = 0
    for k, (x, mx) in enumerate(pairs):
        vx = legendre_symbol(x % p, p) if x % p != 0 else 0
        vmx = legendre_symbol(mx % p, p) if mx % p != 0 else 0
        if vx == 0 or vmx == 0:
            zero_count += 1
        elif vx == vmx:
            even_count += 1
        elif vx == -vmx:
            odd_count += 1
    total_nonzero = even_count + odd_count
    print(f"  chi_{p}: even_pairs={even_count}, odd_pairs={odd_count}, zero_pairs={zero_count}")
    if even_count > odd_count:
        print(f"    => chi_{p} is PREDOMINANTLY EVEN")
    elif odd_count > even_count:
        print(f"    => chi_{p} is PREDOMINANTLY ODD")
    else:
        print(f"    => chi_{p} is MIXED")

# ============================================================
# Stage 6: Gauss sum verification
# ============================================================
print(f"\n{'=' * 80}")
print(f"STAGE 6: Gauss sums")
print("=" * 80)

for p in [3, 5, 7]:
    # g(chi_p) = sum_{a=1}^{p-1} (a/p) * exp(2*pi*i*a/p)
    g = 0.0 + 0.0j
    for a in range(1, p):
        ls = int(legendre_symbol(a, p))
        g += ls * np.exp(2j * np.pi * a / p)
    
    norm_sq = abs(g)**2
    print(f"  g(chi_{p}) = {g:.6f}")
    print(f"  |g(chi_{p})|^2 = {norm_sq:.6f} (expected: {p})")
    assert abs(norm_sq - p) < 1e-10, f"|g(chi_{p})|^2 = {norm_sq}, expected {p}"
    print(f"    PASS: |g(chi_{p})|^2 = {p}")

# Prefactor
C_210 = phi_m / np.sqrt(7)
C_30 = 8 / np.sqrt(5)
print(f"\n  C_210 = phi(210)/sqrt(7) = 48/sqrt(7) = {C_210:.6f}")
print(f"  C_30  = phi(30)/sqrt(5)  = 8/sqrt(5)  = {C_30:.6f}")
print(f"  C_210/C_30 = {C_210/C_30:.6f}")

# ============================================================
# Stage 7: Full spectral analysis
# ============================================================
print(f"\n{'=' * 80}")
print(f"STAGE 7: Full spectral analysis of D_sym(210)")
print("=" * 80)

eig_full = np.sort(np.linalg.eigvalsh(D_sym.astype(float)))[::-1]
print(f"Rank of D_sym: {np.linalg.matrix_rank(D_sym.astype(float))}")
print(f"Largest eigenvalue:  {eig_full[0]:.6f}")
print(f"Smallest eigenvalue: {eig_full[-1]:.6f}")
print(f"Spectral gap (lambda_1/lambda_0): {eig_full[1]/eig_full[0]:.6f}")

# How many distinct eigenvalues?
unique_eigs = []
for e in eig_full:
    if not any(abs(e - u) < 0.01 for u in unique_eigs):
        unique_eigs.append(e)
print(f"Number of approximately distinct eigenvalues: {len(unique_eigs)}")

# ============================================================
# Stage 8: Determinant ratio analysis
# ============================================================
print(f"\n{'=' * 80}")
print(f"STAGE 8: Determinant ratio analysis")
print("=" * 80)

print(f"det(D_even) = {det_even:.6e}")
print(f"det(D_odd)  = {det_odd:.6e}")

if abs(det_odd) > 1:
    ratio = det_even / det_odd
    print(f"det(D_even)/det(D_odd) = {ratio:.6f}")
    
    # For m=30: ratio = -27 = -3^3 = -p_mid^3
    # For m=210: p_mid = 5 (middle of {2,3,5,7}).
    # Is ratio = -5^k for some k?
    # Or: primes are {2,3,5,7}, middle pair is {3,5}.
    # Check various powers
    for base in [3, 5, 7, 15, 21, 35, 105]:
        for sgn in [1, -1]:
            for k in range(1, 8):
                target = sgn * base**k
                if abs(ratio - target) < abs(ratio) * 0.001:
                    print(f"  MATCH: ratio ≈ {'+' if sgn > 0 else '-'}{base}^{k} = {target}")
    
    # Also try ratio as fraction of small numbers
    for num in range(-1000, 1000):
        if num == 0:
            continue
        for den in range(1, 100):
            if abs(ratio - num/den) < abs(ratio) * 0.0001:
                print(f"  FRACTION: ratio ≈ {num}/{den} = {num/den:.6f}")
                break
elif abs(det_odd) < 1e-10:
    print(f"det(D_odd) ≈ 0 — D_odd is SINGULAR")
    print(f"This means the odd block has a zero eigenvalue")
    null_count = np.sum(np.abs(eig_odd) < 1e-6)
    print(f"Null space dimension of D_odd: {null_count}")

# ============================================================
# Stage 9: Summary
# ============================================================
print(f"\n{'=' * 80}")
print(f"SUMMARY: D_sym(210) key results")
print("=" * 80)

print(f"""
  m = {m}
  phi(m) = {phi_m}
  phi(m)/2 = {phi_m//2} (predicted number of exponents for m=210 regime)
  
  D_sym is {phi_m}x{phi_m} integer symmetric matrix
  Block decomposition under x <-> {m}-x:
    D_even: {D_even.shape}, trace = {trace_even:.2f}
    D_odd:  {D_odd.shape}, trace = {trace_odd:.2f}
    Cross-coupling: {cross_norm:.2e} (exact zero)
  
  Eigenvalues:
    D_even: {np.sum(eig_even > 0.01)} positive, {np.sum(np.abs(eig_even) < 0.01)} zero, {np.sum(eig_even < -0.01)} negative
    D_odd:  {np.sum(eig_odd > 0.01)} positive, {np.sum(np.abs(eig_odd) < 0.01)} zero, {np.sum(eig_odd < -0.01)} negative
  
  Determinants:
    det(D_even) = {det_even:.6e}
    det(D_odd)  = {det_odd:.6e}
  
  Legendre characters:
    chi_3(-1) = -1 (3 ≡ 3 mod 4) => ODD
    chi_5(-1) = +1 (5 ≡ 1 mod 4) => EVEN
    chi_7(-1) = -1 (7 ≡ 3 mod 4) => ODD   [KEY PREDICTION CONFIRMED]
  
  Gauss sums:
    |g(chi_3)|^2 = 3  ✓
    |g(chi_5)|^2 = 5  ✓
    |g(chi_7)|^2 = 7  ✓
  
  Predicted prefactor: C_210 = 48/sqrt(7) = {C_210:.6f}
""")

print("=== COMPUTATION COMPLETE ===")
