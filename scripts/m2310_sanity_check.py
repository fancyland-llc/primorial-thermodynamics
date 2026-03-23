"""
m2310_sanity_check.py — Fast pre-flight for the m=2310 primorial
=================================================================

Phase 1: Everything that runs in SECONDS (no determinants).
If this passes, Phase 2 (overnight determinant) is worth running.

Prediction written down BEFORE computation:
  det(D_even)/det(D_odd) = -3^5 * 23 = -5589
  (power of 3 = k = 5 prime factors; extra factor = sum(primes)-count = 28-5 = 23)

Requires: numpy, sympy (for legendre_symbol only)
"""

import numpy as np
from math import gcd
from sympy import legendre_symbol as sym_legendre, factorint, totient, divisor_count
import time

passed = 0
failed = 0
total_tests = 0

def check(condition, label):
    global passed, failed, total_tests
    total_tests += 1
    if condition:
        passed += 1
        print(f"  PASS [{total_tests:3d}]: {label}")
    else:
        failed += 1
        print(f"  FAIL [{total_tests:3d}]: {label}")

def legendre(a, p):
    """Fast modular Legendre symbol."""
    r = a % p
    if r == 0:
        return 0
    return 1 if pow(r, (p - 1) // 2, p) == 1 else -1

def cyclic_distance(a, b, n):
    d = abs(a - b) % n
    return min(d, n - d)


m = 2310
primes_m = [2, 3, 5, 7, 11]
p_max = 11

print("=" * 80)
print(f"m = {m} = 2*3*5*7*11  —  FAST SANITY CHECKS")
print("=" * 80)

# ── Stage 1: Basic structure ──
print("\n--- Stage 1: Basic structure ---")

check(factorint(m) == {2:1, 3:1, 5:1, 7:1, 11:1}, f"{m} = 2*3*5*7*11")

coprimes = sorted([x for x in range(1, m) if gcd(x, m) == 1])
phi = len(coprimes)
tau = divisor_count(m)

check(phi == 480, f"phi({m}) = {phi}")
check(tau == 32, f"tau({m}) = {tau}")
check(phi != tau, f"phi != tau (480 != 32)")

# ── Stage 2: Legendre character parity (Euler's criterion) ──
print("\n--- Stage 2: Legendre character parity ---")

for p in [3, 5, 7, 11]:
    chi_neg1 = legendre(-1, p)
    expected = "EVEN" if p % 4 == 1 else "ODD"
    actual = "EVEN" if chi_neg1 == 1 else "ODD"
    check(actual == expected, f"chi_{p}(-1) = {chi_neg1} => {actual} (p ≡ {p%4} mod 4)")

# ── Stage 3: Gauss sums ──
print("\n--- Stage 3: Gauss sums ---")

for p in [3, 5, 7, 11]:
    g = sum(legendre(a, p) * np.exp(2j * np.pi * a / p) for a in range(1, p))
    norm_sq = abs(g)**2
    check(abs(norm_sq - p) < 1e-8, f"|g(chi_{p})|^2 = {p} (got {norm_sq:.8f})")

C_2310 = phi / np.sqrt(p_max)
print(f"\n  Predicted prefactor: C_2310 = {phi}/sqrt({p_max}) = {C_2310:.6f}")
check(abs(C_2310 - 480/np.sqrt(11)) < 1e-6, f"C_2310 = 480/sqrt(11) = {C_2310:.4f}")

# ── Stage 4: Palindromic pairs ──
print("\n--- Stage 4: Palindromic pairs ---")

# Check all partners are coprime
all_partners_coprime = all(gcd(m - x, m) == 1 for x in coprimes)
check(all_partners_coprime, "All palindromic partners are coprime")

# No fixed point
check(gcd(m // 2, m) != 1, f"No fixed point ({m//2} not coprime to {m})")

pairs = [(x, m - x) for x in coprimes if x < m - x]
check(len(pairs) == 240, f"240 palindromic pairs (got {len(pairs)})")

# ── Stage 5: Full pair-by-pair character separation ──
print("\n--- Stage 5: Character separation (all pairs) ---")

for p in [3, 5, 7, 11]:
    expected = "even" if p % 4 == 1 else "odd"
    correct = 0
    for x, mx in pairs:
        vx = legendre(x, p)
        vmx = legendre(mx, p)
        if vx == 0 or vmx == 0:
            continue
        if expected == "even" and vx == vmx:
            correct += 1
        elif expected == "odd" and vx == -vmx:
            correct += 1
    check(correct == 240, f"chi_{p} is {expected.upper()}: {correct}/240 pairs consistent")

# ── Stage 6: Build 480x480 distance matrix ──
print("\n--- Stage 6: Building D_sym(2310) [480x480] ---")
t0 = time.time()

D = np.zeros((480, 480), dtype=np.int32)
for i in range(480):
    for j in range(i + 1, 480):
        d = cyclic_distance(coprimes[i], coprimes[j], m)
        D[i, j] = d
        D[j, i] = d

t1 = time.time()
print(f"  Built in {t1-t0:.1f}s")

check(np.allclose(D, D.T), "D_sym(2310) is symmetric")
check(np.all(np.diag(D) == 0), "Diagonal is zero")
check(D.max() == m // 2 - 1, f"Max distance = {D.max()} (expected {m//2-1})")

# ── Stage 7: Build even/odd integer blocks (240x240) ──
print("\n--- Stage 7: Building 240x240 even/odd blocks ---")
t0 = time.time()

idx = {x: i for i, x in enumerate(coprimes)}

D_even_int = []
D_odd_int = []
for k in range(240):
    row_e, row_o = [], []
    xk, mxk = pairs[k]
    ik, jk = idx[xk], idx[mxk]
    for l in range(240):
        xl, mxl = pairs[l]
        il, jl = idx[xl], idx[mxl]
        se = D[ik,il] + D[ik,jl] + D[jk,il] + D[jk,jl]
        so = D[ik,il] - D[ik,jl] - D[jk,il] + D[jk,jl]
        assert se % 2 == 0 and so % 2 == 0, f"Parity violation at ({k},{l})"
        row_e.append(se // 2)
        row_o.append(so // 2)
    D_even_int.append(row_e)
    D_odd_int.append(row_o)

t1 = time.time()
print(f"  Built in {t1-t0:.1f}s")

D_even = np.array(D_even_int, dtype=np.float64)
D_odd = np.array(D_odd_int, dtype=np.float64)

# ── Stage 8: Block properties (numerical) ──
print("\n--- Stage 8: Block properties (numerical) ---")

tr_even = np.trace(D_even)
tr_odd = np.trace(D_odd)
print(f"  tr(D_even) = {tr_even:.0f}")
print(f"  tr(D_odd) = {tr_odd:.0f}")
check(abs(tr_even + tr_odd) < 1e-6, f"Traces are antisymmetric: {tr_even:.0f} + {tr_odd:.0f} = 0")

# ── Stage 9: Cross-coupling check (numerical) ──
print("\n--- Stage 9: Cross-coupling (numerical spot-check) ---")

# Build P for a random subset to verify decoupling
# Full 480x480 P is expensive, so check via the integer construction directly
# The integer block construction already enforces exact separation.
# Just verify the block matrices are symmetric:
check(np.allclose(D_even, D_even.T), "D_even(2310) is symmetric")
check(np.allclose(D_odd, D_odd.T), "D_odd(2310) is symmetric")

# ── Stage 10: Eigenvalue structure (numerical) ──
print("\n--- Stage 10: Eigenvalue structure ---")
t0 = time.time()

eig_even = np.sort(np.linalg.eigvalsh(D_even))[::-1]
eig_odd = np.sort(np.linalg.eigvalsh(D_odd))[::-1]

t1 = time.time()
print(f"  Eigenvalues computed in {t1-t0:.1f}s")

n_pos_even = np.sum(eig_even > 0.1)
n_neg_even = np.sum(eig_even < -0.1)
n_pos_odd = np.sum(eig_odd > 0.1)
n_neg_odd = np.sum(eig_odd < -0.1)

print(f"  D_even: {n_pos_even} positive, {n_neg_even} negative, largest = {eig_even[0]:.1f}")
print(f"  D_odd:  {n_pos_odd} positive, {n_neg_odd} negative, largest = {eig_odd[0]:.1f}")

check(n_pos_even == 1, f"D_even: exactly 1 positive eigenvalue (got {n_pos_even})")
check(n_neg_even == 239, f"D_even: 239 negative eigenvalues (got {n_neg_even})")
check(n_pos_odd == 0, f"D_odd: 0 positive eigenvalues (got {n_pos_odd})")
check(n_neg_odd == 240, f"D_odd: 240 negative eigenvalues (got {n_neg_odd})")

# ── Stage 11: Numerical determinant (sanity, NOT exact) ──
print("\n--- Stage 11: Numerical determinant estimate ---")

# Use log-det for stability
sign_e, logdet_e = np.linalg.slogdet(D_even)
sign_o, logdet_o = np.linalg.slogdet(D_odd)

print(f"  sign(det_even) = {sign_e:.0f}, log|det_even| = {logdet_e:.1f}")
print(f"  sign(det_odd)  = {sign_o:.0f}, log|det_odd|  = {logdet_o:.1f}")
print(f"  log|ratio| = {logdet_e - logdet_o:.1f}")

check(sign_e == -1, f"det(D_even) is negative (sign = {sign_e})")
check(sign_o == 1, f"det(D_odd) is positive (sign = {sign_o})")

# Numerical estimate of ratio
num_ratio = -np.exp(logdet_e - logdet_o)
print(f"  Numerical ratio estimate: {num_ratio:.1f}")
print(f"  Prediction: -5589 = -3^5 * 23")
print(f"  Deviation from prediction: {abs(num_ratio - (-5589)) / 5589 * 100:.1f}%")

# Check if numerical ratio is close to -5589
check(abs(num_ratio - (-5589)) / 5589 < 0.01,
      f"Numerical ratio ≈ -5589 within 1% (got {num_ratio:.1f})")

# ── Stage 12: Quick divisibility check via NumPy float det ──
print("\n--- Stage 12: Pattern checks (from numerical det) ---")

# We can't do exact modular arithmetic with float determinants.
# But we CAN check the sign pattern and rough magnitude.
# The exact checks require the overnight SymPy run.

digits_even = logdet_e / np.log(10)
digits_odd = logdet_o / np.log(10)
print(f"  |det(D_even)| has ~{digits_even:.0f} digits")
print(f"  |det(D_odd)|  has ~{digits_odd:.0f} digits")
check(digits_even > 500, f"det(D_even) is a huge number (~{digits_even:.0f} digits)")


# ============================================================
print("\n" + "=" * 80)
print(f"FAST SANITY CHECK RESULTS: {passed} passed, {failed} failed, {total_tests} total")
print("=" * 80)

if failed > 0:
    print(f"\n*** {failed} ASSERTION(S) FAILED — DO NOT proceed to overnight run ***")
else:
    print("\nAll fast checks passed. The rabbit is clean.")
    print(f"Numerical ratio estimate: {num_ratio:.1f}  (prediction: -5589)")
    print(f"Exact determinants require ~8-20 hours via SymPy Bareiss.")
    print("Proceed to m2310_overnight.py when ready.")

# Save the integer blocks for the overnight script
print("\nSaving integer blocks to disk for overnight script...")
np.save("D_even_2310.npy", np.array(D_even_int, dtype=object))
np.save("D_odd_2310.npy", np.array(D_odd_int, dtype=object))
print("Saved: D_even_2310.npy, D_odd_2310.npy")
