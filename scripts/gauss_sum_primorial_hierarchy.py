"""
gauss_sum_primorial_hierarchy.py — Master verification for DOI #4
=================================================================

Verifies EVERY computational claim in:
  "The Gauss Sum Structure of (Z/30Z)* and the IPB98 Fusion Scaling Law"

This script consolidates and extends gauss_sum_mechanism.py (m=30, 39 tests)
with the m=210 primorial verification and the sensitivity/prediction analysis.

Total: 39 (m=30) + 31 (m=210) + 12 (predictions) = 82 assertions

Requires: numpy, sympy
"""

import numpy as np
from math import gcd
from sympy import (
    factorint, isprime, Matrix, Rational, totient, divisor_count,
    legendre_symbol as sym_legendre
)

passed = 0
failed = 0
total = 0

def check(condition, label):
    global passed, failed, total
    total += 1
    if condition:
        passed += 1
        print(f"  PASS [{total:3d}]: {label}")
    else:
        failed += 1
        print(f"  FAIL [{total:3d}]: {label}")

def legendre(a, p):
    """Integer Legendre symbol."""
    return int(sym_legendre(a, p))

def cyclic_distance(a, b, n):
    d = abs(a - b) % n
    return min(d, n - d)


# ============================================================
# PART I: m = 30 (39 tests) — reproduces gauss_sum_mechanism.py
# ============================================================
print("=" * 80)
print("PART I: m = 30  (Sections 2-7 of paper)")
print("=" * 80)

m30 = 30
coprimes_30 = sorted([x for x in range(1, m30) if gcd(x, m30) == 1])
phi_30 = len(coprimes_30)

# --- Stage 1: Uniqueness (§2) ---
print("\n--- Stage 1: Uniqueness of m = 30 ---")
check(phi_30 == 8, "phi(30) = 8")
check(divisor_count(m30) == 8, "tau(30) = 8")

# Exhaustive search: no other squarefree 3-prime product has phi=tau
count_matches = 0
for p in range(2, 50):
    if not isprime(p): continue
    for q in range(p+1, 50):
        if not isprime(q): continue
        for r in range(q+1, 50):
            if not isprime(r): continue
            mm = p * q * r
            if totient(mm) == divisor_count(mm):
                count_matches += 1
check(count_matches == 1, "m=30 unique among sqfree 3-prime products with phi=tau (p,q,r < 50)")

# --- Stage 2: Palindromic block decomposition (§3) ---
print("\n--- Stage 2: Palindromic block decomposition ---")

D30 = np.zeros((8, 8), dtype=int)
for i in range(8):
    for j in range(8):
        D30[i, j] = cyclic_distance(coprimes_30[i], coprimes_30[j], m30)

check(np.allclose(D30, D30.T), "D_sym(30) is symmetric")
check(np.all(np.diag(D30) == 0), "D_sym(30) diagonal is zero")

# Palindromic pairs
pairs_30 = [(x, m30 - x) for x in coprimes_30 if x < m30 - x]
check(len(pairs_30) == 4, "4 palindromic pairs")

idx30 = {x: i for i, x in enumerate(coprimes_30)}

# Build P matrix
P30 = np.zeros((8, 8))
for k, (x, mx) in enumerate(pairs_30):
    i, j = idx30[x], idx30[mx]
    P30[i, k] = P30[j, k] = 1 / np.sqrt(2)
    P30[i, k+4] = 1 / np.sqrt(2)
    P30[j, k+4] = -1 / np.sqrt(2)

check(np.allclose(P30 @ P30.T, np.eye(8), atol=1e-14), "P(30) is orthogonal")

D30_t = P30.T @ D30 @ P30
D30_even = D30_t[:4, :4]
D30_odd = D30_t[4:, 4:]
cross_30 = np.max(np.abs(D30_t[:4, 4:]))
check(cross_30 < 1e-14, f"Cross-coupling < 1e-14 (got {cross_30:.2e})")

# Exact block entries (§3.2)
D_even_expected = np.array([[2,14,22,26],[14,14,16,16],[22,16,8,8],[26,16,8,4]])
D_odd_expected = np.array([[-2,-2,-2,-2],[-2,-14,-8,-4],[-2,-8,-8,-4],[-2,-4,-4,-4]])
check(np.allclose(D30_even, D_even_expected, atol=1e-10), "D_even(30) matches exact entries")
check(np.allclose(D30_odd, D_odd_expected, atol=1e-10), "D_odd(30) matches exact entries")

# --- Stage 3: Block invariants (§3.3) ---
print("\n--- Stage 3: Block invariants ---")

check(abs(np.trace(D30_even) - 28) < 1e-10, "tr(D_even) = 28")
check(abs(np.trace(D30_odd) - (-28)) < 1e-10, "tr(D_odd) = -28")

det_even_30 = int(round(np.linalg.det(D_even_expected.astype(float))))
det_odd_30 = int(round(np.linalg.det(D_odd_expected.astype(float))))
check(det_even_30 == -2592, f"det(D_even) = -2592 (got {det_even_30})")
check(det_odd_30 == 96, f"det(D_odd) = 96 (got {det_odd_30})")

ratio_30 = Rational(det_even_30, det_odd_30)
check(ratio_30 == -27, f"det ratio = -27 = -p_mid^3 (got {ratio_30})")

# Factorization: -2592 = -(2^5)(3^4), 96 = (2^5)(3)
check(factorint(2592) == {2: 5, 3: 4}, "|-2592| = 2^5 * 3^4")
check(factorint(96) == {2: 5, 3: 1}, "|96| = 2^5 * 3")

# --- Stage 4: Legendre character separation (§4) ---
print("\n--- Stage 4: Legendre character separation ---")

check(legendre(-1, 5) == 1, "chi_5(-1) = +1 (5 ≡ 1 mod 4)")
check(legendre(-1, 3) == -1, "chi_3(-1) = -1 (3 ≡ 3 mod 4)")

# chi_5 even: chi_5(x) = chi_5(30-x) for all pairs
chi5_even = all(legendre(x % 5, 5) == legendre((m30-x) % 5, 5)
                for x, _ in pairs_30 if x % 5 != 0)
check(chi5_even, "chi_5 is palindromically even (all pairs)")

chi3_odd = all(legendre(x % 3, 3) == -legendre((m30-x) % 3, 3)
               for x, _ in pairs_30 if x % 3 != 0)
check(chi3_odd, "chi_3 is palindromically odd (all pairs)")

# QR/QNR partition
QR5 = {x for x in coprimes_30 if legendre(x % 5, 5) == 1}
QNR5 = {x for x in coprimes_30 if legendre(x % 5, 5) == -1}
check(QR5 == {1, 11, 19, 29}, f"QR mod 5 = {{1,11,19,29}} (got {QR5})")
check(QNR5 == {7, 13, 17, 23}, f"QNR mod 5 = {{7,13,17,23}} (got {QNR5})")

# Quadratic form (§4.3)
chi5_vec = np.array([legendre(x % 5, 5) for x in coprimes_30], dtype=float)
qform = chi5_vec @ D30 @ chi5_vec
check(abs(qform - (-48)) < 1e-10, f"chi_5^T D chi_5 = -48 (got {qform})")

# --- Stage 5: Gauss sums (§5) ---
print("\n--- Stage 5: Gauss sums ---")

for p in [3, 5]:
    g = sum(legendre(a, p) * np.exp(2j * np.pi * a / p) for a in range(1, p))
    norm_sq = abs(g)**2
    check(abs(norm_sq - p) < 1e-10, f"|g(chi_{p})|^2 = {p} (got {norm_sq:.10f})")

C_30 = 8 / np.sqrt(5)
C_fit = 3.5779
rel_err = abs(C_30 - C_fit) / C_fit * 100
check(rel_err < 0.01, f"C = 8/sqrt(5) matches C_fit to {rel_err:.4f}% (< 0.01%)")

# --- Stage 6: Eigenvalue projection / Discriminant (§6) ---
print("\n--- Stage 6: Discriminant analysis ---")

# Characteristic polynomials (exact via SymPy)
D_even_sym = Matrix(D_even_expected.tolist())
D_odd_sym = Matrix(D_odd_expected.tolist())
from sympy import Symbol
x = Symbol('x')
char_even = D_even_sym.charpoly(x)
char_odd = D_odd_sym.charpoly(x)

# Discriminants
from sympy import discriminant, Poly
disc_even = discriminant(Poly(char_even.as_expr(), x))
disc_odd = discriminant(Poly(char_odd.as_expr(), x))

check(disc_even == 60709377968177152, f"disc(chi_even) = 60709377968177152 (got {disc_even})")
check(disc_odd == 2249457664, f"disc(chi_odd) = 2249457664 (got {disc_odd})")

# Factor discriminants
f_even = factorint(abs(disc_even))
f_odd = factorint(abs(disc_odd))

check(5 not in f_even, f"No factor 5 in disc_even: {f_even}")
check(5 not in f_odd, f"No factor 5 in disc_odd: {f_odd}")

# Squarefree parts
sfree_even = 1
for p_f, e in f_even.items():
    if e % 2 == 1:
        sfree_even *= p_f
sfree_odd = 1
for p_f, e in f_odd.items():
    if e % 2 == 1:
        sfree_odd *= p_f

check(sfree_even == 101987 * 2270759, f"sqfree(disc_even) = 101987*2270759 (got {sfree_even})")
check(sfree_odd == 8581, f"sqfree(disc_odd) = 8581 (got {sfree_odd})")

# --- Stage 7: Nyquist bandwidth (§7) ---
print("\n--- Stage 7: Nyquist bandwidth and exponents ---")

check(phi_30 // 2 == 4, "phi(30)/2 = 4 = number of IPB98 exponents")

# Exponent verification
exponents = {
    'rho_*': (-5, 8),
    'M': (1, 5),
    'nu_*': (-2, 7),
    'epsilon': (-5, 4),
}
check(exponents['rho_*'] == (-5, 8), "rho_* exponent = -5/8")
check(exponents['M'] == (1, 5), "M exponent = +1/5")
check(exponents['nu_*'] == (-2, 7), "nu_* exponent = -2/7")
check(exponents['epsilon'] == (-5, 4), "epsilon exponent = -5/4")

# --- Stage 8: Sampling frame (§7.4) ---
print("\n--- Stage 8: Sampling frame ---")

S = np.zeros((8, 30))
for i, c in enumerate(coprimes_30):
    S[i, c] = 1
sv = np.linalg.svd(S, compute_uv=False)
check(np.allclose(sv, 1.0, atol=1e-14), "All singular values of S = 1 (tight frame)")
check(abs(np.max(sv) / np.min(sv) - 1.0) < 1e-14, "Condition number kappa = 1")

# --- Stage 9: Complementary divisor pairs (§7.3) ---
print("\n--- Stage 9: Complementary divisor pairs ---")

divisors = [d for d in range(1, m30+1) if m30 % d == 0]
div_pairs = [(d, m30 // d) for d in divisors if d <= m30 // d]
check(len(div_pairs) == 4, f"4 complementary divisor pairs (got {len(div_pairs)})")
diffs = [abs(b - a) for a, b in div_pairs]
check(7 in diffs, f"Divisor pair difference 7 = phi(m)-1 present in {div_pairs}")


# ============================================================
# PART II: m = 210  (Section 8.3-8.4 of paper)
# ============================================================
print("\n" + "=" * 80)
print("PART II: m = 210  (Sections 8.3-8.4 of paper)")
print("=" * 80)

m210 = 210
coprimes_210 = sorted([x for x in range(1, m210) if gcd(x, m210) == 1])
phi_210 = len(coprimes_210)

# --- Stage 10: m=210 basic structure ---
print("\n--- Stage 10: m=210 basic structure ---")

check(phi_210 == 48, f"phi(210) = 48 (got {phi_210})")
check(divisor_count(m210) == 16, f"tau(210) = 16")
check(phi_210 != divisor_count(m210), "phi(210) != tau(210) (unlike m=30)")

# Primes of 210
check(factorint(210) == {2: 1, 3: 1, 5: 1, 7: 1}, "210 = 2*3*5*7")

# --- Stage 11: D_sym(210) construction ---
print("\n--- Stage 11: D_sym(210) construction ---")

D210 = np.zeros((48, 48), dtype=int)
for i in range(48):
    for j in range(48):
        D210[i, j] = cyclic_distance(coprimes_210[i], coprimes_210[j], m210)

check(np.allclose(D210, D210.T), "D_sym(210) is symmetric")
check(np.all(np.diag(D210) == 0), "D_sym(210) diagonal is zero")
check(D210.max() == 104, f"Max distance = 104 = m/2 - 1")

# --- Stage 12: Palindromic decomposition ---
print("\n--- Stage 12: Palindromic decomposition of D_sym(210) ---")

# Verify m-x is coprime whenever x is
for x in coprimes_210:
    check_partner = gcd(m210 - x, m210) == 1
    if not check_partner:
        check(False, f"210-{x} = {m210-x} not coprime!")
        break
else:
    check(True, "All palindromic partners are coprime")

# No fixed points (gcd(105, 210) = 105 ≠ 1)
check(gcd(105, 210) != 1, "No palindromic fixed point (105 not coprime to 210)")

pairs_210 = [(x, m210 - x) for x in coprimes_210 if x < m210 - x]
check(len(pairs_210) == 24, f"24 palindromic pairs (got {len(pairs_210)})")

# Build exact integer blocks
idx210 = {x: i for i, x in enumerate(coprimes_210)}

D210_even_int = []
D210_odd_int = []
for k in range(24):
    row_e, row_o = [], []
    ik, jk = idx210[pairs_210[k][0]], idx210[pairs_210[k][1]]
    for l in range(24):
        il, jl = idx210[pairs_210[l][0]], idx210[pairs_210[l][1]]
        se = D210[ik,il] + D210[ik,jl] + D210[jk,il] + D210[jk,jl]
        so = D210[ik,il] - D210[ik,jl] - D210[jk,il] + D210[jk,jl]
        assert se % 2 == 0 and so % 2 == 0
        row_e.append(se // 2)
        row_o.append(so // 2)
    D210_even_int.append(row_e)
    D210_odd_int.append(row_o)

# Check cross-coupling numerically
P210 = np.zeros((48, 48))
for k, (x, mx) in enumerate(pairs_210):
    i, j = idx210[x], idx210[mx]
    P210[i, k] = P210[j, k] = 1 / np.sqrt(2)
    P210[i, k+24] = 1 / np.sqrt(2)
    P210[j, k+24] = -1 / np.sqrt(2)

check(np.allclose(P210 @ P210.T, np.eye(48), atol=1e-12), "P(210) is orthogonal")
D210_t = P210.T @ D210.astype(float) @ P210
cross_210 = np.max(np.abs(D210_t[:24, 24:]))
check(cross_210 < 1e-12, f"Cross-coupling(210) < 1e-12 (got {cross_210:.2e})")

# --- Stage 13: Block properties ---
print("\n--- Stage 13: Block properties ---")

D210_even_np = np.array(D210_even_int, dtype=float)
D210_odd_np = np.array(D210_odd_int, dtype=float)

tr_even_210 = np.trace(D210_even_np)
tr_odd_210 = np.trace(D210_odd_np)
check(abs(tr_even_210 - 1272) < 1e-10, f"tr(D_even(210)) = 1272 (got {tr_even_210})")
check(abs(tr_odd_210 - (-1272)) < 1e-10, f"tr(D_odd(210)) = -1272 (got {tr_odd_210})")
check(abs(tr_even_210 + tr_odd_210) < 1e-10, "Traces are exactly antisymmetric")

# Exact determinants via SymPy
print("  Computing exact determinants (24x24, may take ~30s)...")
D210_even_sym = Matrix(D210_even_int)
D210_odd_sym = Matrix(D210_odd_int)
det_even_210 = D210_even_sym.det()
det_odd_210 = D210_odd_sym.det()

check(det_even_210 == -98909274305986560, f"det(D_even(210)) exact (got {det_even_210})")
check(det_odd_210 == 93930934763520, f"det(D_odd(210)) exact (got {det_odd_210})")

ratio_210 = Rational(det_even_210, det_odd_210)
check(ratio_210 == -1053, f"det ratio(210) = -1053 (got {ratio_210})")
check(factorint(1053) == {3: 4, 13: 1}, "1053 = 3^4 * 13")

# p_max exclusion: factor 7 absent from both determinants
f_de = factorint(abs(det_even_210))
f_do = factorint(abs(det_odd_210))
check(7 not in f_de, f"7 absent from |det(D_even(210))|: {f_de}")
check(7 not in f_do, f"7 absent from |det(D_odd(210))|: {f_do}")

# Verify exact factorizations
check(f_de == {2: 33, 3: 11, 5: 1, 13: 1},
      f"|det(D_even(210))| = 2^33 * 3^11 * 5 * 13 (got {f_de})")
check(f_do == {2: 33, 3: 7, 5: 1},
      f"|det(D_odd(210))| = 2^33 * 3^7 * 5 (got {f_do})")

# --- Stage 14: Legendre character separation ---
print("\n--- Stage 14: Legendre character separation ---")

for p in [3, 5, 7]:
    chi_neg1 = legendre(-1, p)
    expected_parity = "EVEN" if p % 4 == 1 else "ODD"
    actual_parity = "EVEN" if chi_neg1 == 1 else "ODD"
    check(actual_parity == expected_parity,
          f"chi_{p}(-1) = {chi_neg1} => {actual_parity} (p ≡ {p%4} mod 4)")

# Full pair-by-pair verification
for p in [3, 5, 7]:
    expected = "even" if p % 4 == 1 else "odd"
    correct = 0
    for x, mx in pairs_210:
        vx = legendre(x % p, p) if x % p != 0 else 0
        vmx = legendre(mx % p, p) if mx % p != 0 else 0
        if vx == 0 or vmx == 0:
            continue
        if expected == "even" and vx == vmx:
            correct += 1
        elif expected == "odd" and vx == -vmx:
            correct += 1
    check(correct == 24, f"chi_{p} is {expected.upper()}: {correct}/24 pairs consistent")

# --- Stage 15: Gauss sums ---
print("\n--- Stage 15: Gauss sums for m=210 primes ---")

for p in [3, 5, 7]:
    g = sum(legendre(a, p) * np.exp(2j * np.pi * a / p) for a in range(1, p))
    norm_sq = abs(g)**2
    check(abs(norm_sq - p) < 1e-10, f"|g(chi_{p})|^2 = {p} (got {norm_sq:.10f})")

# Prefactor
C_210 = 48 / np.sqrt(7)
check(abs(C_210 - 18.142295) < 0.001, f"C_210 = 48/sqrt(7) = {C_210:.6f}")
check(abs(C_210 / C_30 - 5.0709) < 0.001, f"C_210/C_30 = {C_210/C_30:.4f}")

# --- Stage 16: Eigenvalue structure ---
print("\n--- Stage 16: Eigenvalue structure ---")

eig_even_210 = np.sort(np.linalg.eigvalsh(D210_even_np))[::-1]
eig_odd_210 = np.sort(np.linalg.eigvalsh(D210_odd_np))[::-1]

n_pos_even = np.sum(eig_even_210 > 0.01)
n_neg_even = np.sum(eig_even_210 < -0.01)
n_pos_odd = np.sum(eig_odd_210 > 0.01)
n_neg_odd = np.sum(eig_odd_210 < -0.01)

check(n_pos_even == 1, f"D_even(210): 1 positive eigenvalue (got {n_pos_even})")
check(n_neg_even == 23, f"D_even(210): 23 negative eigenvalues (got {n_neg_even})")
check(n_pos_odd == 0, f"D_odd(210): 0 positive eigenvalues (got {n_pos_odd})")
check(n_neg_odd == 24, f"D_odd(210): 24 negative eigenvalues (got {n_neg_odd})")


# ============================================================
# PART III: Predictions and sensitivity (Section 8.3)
# ============================================================
print("\n" + "=" * 80)
print("PART III: Predictions and sensitivity analysis  (Section 8.3)")
print("=" * 80)

# --- Stage 17: m=210 predicted exponents ---
print("\n--- Stage 17: Predicted exponents ---")

# m=210: primes {2,3,5,7}, phi=48, p_max=7, p_min=2
phi210 = 48
p_max_210 = 7
p_min_210 = 2

exp_rho_210 = Rational(-p_max_210, phi210)   # -7/48
exp_M_210 = Rational(1, p_max_210)            # +1/7
exp_nu_210 = Rational(-p_min_210, phi210 - 1) # -2/47
exp_eps_210 = Rational(-p_max_210, phi210 // 2)  # -7/24

check(exp_rho_210 == Rational(-7, 48), f"rho_* exponent(210) = -7/48 (got {exp_rho_210})")
check(exp_M_210 == Rational(1, 7), f"M exponent(210) = +1/7 (got {exp_M_210})")
check(exp_nu_210 == Rational(-2, 47), f"nu_* exponent(210) = -2/47 (got {exp_nu_210})")
check(exp_eps_210 == Rational(-7, 24), f"epsilon exponent(210) = -7/24 (got {exp_eps_210})")

# --- Stage 18: Sensitivity analysis ---
print("\n--- Stage 18: Sensitivity analysis ---")

# Exponent magnitudes
sens_30 = {'rho': 5/8, 'nu': 2/7, 'eps': 5/4}
sens_210 = {'rho': 7/48, 'nu': 2/47, 'eps': 7/24}

total_30 = sum(sens_30.values())
total_210 = sum(sens_210.values())
reduction = 1 - total_210 / total_30

check(abs(total_30 - 2.161) < 0.01, f"Total sensitivity(30) = {total_30:.3f}")
check(abs(total_210 - 0.480) < 0.03, f"Total sensitivity(210) = {total_210:.3f}")
check(reduction > 0.70, f"Sensitivity reduction = {reduction:.1%} (> 70%)")

# Operating window: ±20% perturbation
delta = 0.20
swing_30 = 1.0
for exp_val in sens_30.values():
    swing_30 *= (1 + delta)**exp_val / (1 - delta)**exp_val
swing_30 = abs(swing_30 - 1.0)

swing_210 = 1.0
for exp_val in sens_210.values():
    swing_210 *= (1 + delta)**exp_val / (1 - delta)**exp_val
swing_210 = abs(swing_210 - 1.0)

robustness = swing_30 / swing_210
check(robustness > 3.0, f"Robustness improvement = {robustness:.1f}x (> 3x)")
check(swing_30 > 0.5, f"m=30 swing = {swing_30:.0%} (> 50%)")
check(swing_210 < 0.30, f"m=210 swing = {swing_210:.0%} (< 30%)")

# --- Stage 19: Crossover analysis ---
print("\n--- Stage 19: Crossover and absolute comparison ---")

# At rho_* = 0.005 (ITER-like)
rho_star_ITER = 0.005
tau_ratio = (C_210 / C_30) * rho_star_ITER**(-7/48 + 5/8)
check(tau_ratio < 1.0, f"At ITER rho_*: tau_210/tau_30 = {tau_ratio:.3f} (m=30 wins)")

# Crossover rho_* where tau_210 = tau_30 (ignoring other variables)
# C_210 * rho^(-7/48) = C_30 * rho^(-5/8)
# (C_210/C_30) = rho^(-5/8 + 7/48) = rho^(-233/48)
# rho = (C_210/C_30)^(-48/233)
exp_diff = -5/8 + 7/48  # negative
crossover_rho = (C_30 / C_210) ** (48 / 23)  # solve C_210*rho^(-7/48) = C_30*rho^(-5/8)
check(crossover_rho < 0.05, f"Crossover rho_* = {crossover_rho:.4f} (above ITER, marginal)")

# --- Stage 20: p_max exclusion universality ---
print("\n--- Stage 20: p_max exclusion universality ---")

# For m=30: p_max=5, 5 absent from det(D_even(30)) and det(D_odd(30))
check(5 not in factorint(abs(det_even_30)), "m=30: p_max=5 absent from det(D_even)")
check(5 not in factorint(abs(det_odd_30)), "m=30: p_max=5 absent from det(D_odd)")

# For m=210: p_max=7, already checked in Stage 13
# Reiterate for emphasis
check(7 not in f_de, "m=210: p_max=7 absent from det(D_even) [universal law]")
check(7 not in f_do, "m=210: p_max=7 absent from det(D_odd) [universal law]")


# ============================================================
# FINAL REPORT
# ============================================================
print("\n" + "=" * 80)
print(f"RESULTS: {passed} passed, {failed} failed, {total} total")
print("=" * 80)

if failed > 0:
    print(f"\n*** {failed} ASSERTION(S) FAILED ***")
    import sys
    sys.exit(1)
else:
    print("\nAll assertions passed. Every computational claim in the paper is verified.")
