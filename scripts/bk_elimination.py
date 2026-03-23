"""
bk_elimination.py — Proof by inference: systematic elimination of b_k hypotheses.
Compute what b_k ISN'T. When the topology collapses, the rabbit is caught.

Known data:
  a = [1, 3, 55, 2107, 217855]  (a_k = leading coeff of level-k polynomial)
  b = [-2, -4, -84, -3372, ?]   (b_k = linear coeff of level-k polynomial)
  p = [3, 5, 7, 11, 13]         (primes added at each level)

Proved: a_k = a_{k-1}*p_k^2 + b_{k-1}*p_k
Proved: a_k = 1 - 2*R(m_k) where R(m_k) is the primorial ratio
Proved: -2*R(S_k ∪ {q}) = a_k*q^2 + b_k*q - 1 (universal constant -1)
"""

import sympy as sp
from sympy import Matrix, isprime, factorint
from fractions import Fraction
import math
import time

# ============================================================
# KNOWN DATA
# ============================================================
a = [1, 3, 55, 2107, 217855]
b = [-2, -4, -84, -3372]  # b_4 unknown
p = [3, 5, 7, 11, 13]     # primes at each level
R_primorial = [0, -1, -1, -27, -1053, -108927]  # R(m_0)=R(1)?, R(m_1)=R(6)=-1, ...

# Actually: R(m_k) = -(a_k - 1)/2
# k=0: a_0=1, R=0  (trivial, m=2, phi=1)
# k=1: a_1=3, R=-1  (m=6)
# k=2: a_2=55, R=-27  (m=30) 
# k=3: a_3=2107, R=-1053  (m=210)
# k=4: a_4=217855(?), but a_4 = a_3*13^2 + b_3*13 = 2107*169 + (-3372)*13

a4_check = 2107*169 + (-3372)*13
print(f"a_4 = {a4_check}")
print(f"R(m_5=2310) from a_4: {-(a4_check - 1)//2} (should be -108927)")
print(f"  Match: {-(a4_check-1)/2 == -108927}")

c = [a[k]//p[k] if k > 0 else None for k in range(5)]
print(f"\nc_k = a_k / p_k: {c}")

passes = 0
fails = 0
total = 0

def test(name, condition):
    global passes, fails, total
    total += 1
    if condition:
        passes += 1
        print(f"  [KILL] {name}")
    else:
        print(f"  [SURVIVES] {name}")
        fails += 1

# ============================================================
print("\n" + "="*70)
print("ROUND 1: b_k = f(a_k) for some simple f?")
print("="*70)

# Test b_k = x*a_k + y (affine)
# k=0: -2 = x + y
# k=1: -4 = 3x + y  => 2x = -2, x=-1, y=-1
x_lin, y_lin = -1, -1
print(f"  Affine fit from k=0,1: b = {x_lin}*a + {y_lin}")
for k in range(4):
    pred = x_lin * a[k] + y_lin
    print(f"    k={k}: pred={pred}, actual={b[k]}, {'✓' if pred==b[k] else '✗'}")
test("b_k = -a_k - 1", all(x_lin*a[k]+y_lin == b[k] for k in range(4)))

# Test b_k = x*a_k^2 + y*a_k + z (quadratic in a)
# 3 unknowns, use k=0,1,2
A_q = Matrix([[1, 1, 1], [9, 3, 1], [3025, 55, 1]])
b_q = Matrix([-2, -4, -84])
sol_q = A_q.solve(b_q)
print(f"\n  Quadratic fit b = {sol_q[0]}*a^2 + {sol_q[1]}*a + {sol_q[2]}")
pred_q3 = sol_q[0]*2107**2 + sol_q[1]*2107 + sol_q[2]
print(f"    Predicts b_3 = {pred_q3}, actual = {b[3]}")
test("b_k = quadratic in a_k", pred_q3 == b[3])

# ============================================================
print("\n" + "="*70)
print("ROUND 2: b_k = constant Markov on (a_{k-1}, b_{k-1})?")
print("="*70)
# b_k = c*a_{k-1}*p_k^2 + d*b_{k-1}*p_k
# Same structure as a-recurrence but different constants
# k=1: -4 = 9c - 6d
# k=2: -84 = 75c - 20d
A_m = Matrix([[9, -6], [75, -20]])
b_m = Matrix([-4, -84])
sol_m = A_m.solve(b_m)
print(f"  b_k = {sol_m[0]}*a_{{k-1}}*p_k^2 + {sol_m[1]}*b_{{k-1}}*p_k")
pred_b3 = sol_m[0]*a[2]*p[2]**2 + sol_m[1]*b[2]*p[2]
print(f"  Predicts b_3 = {pred_b3}, actual = {b[3]}")
test("b_k = c*a_{k-1}*p^2 + d*b_{k-1}*p (constant c,d)", pred_b3 == b[3])

# Also try: b_k = c*a_{k-1} + d*b_{k-1} (no prime factors)
A_m2 = Matrix([[1, -2], [3, -4]])
b_m2 = Matrix([-4, -84])
sol_m2 = A_m2.solve(b_m2)
print(f"\n  b_k = {sol_m2[0]}*a_{{k-1}} + {sol_m2[1]}*b_{{k-1}}")
pred_b3_v2 = sol_m2[0]*a[2] + sol_m2[1]*b[2]
print(f"  Predicts b_3 = {pred_b3_v2}, actual = {b[3]}")
test("b_k = c*a_{k-1} + d*b_{k-1} (no primes)", pred_b3_v2 == b[3])

# Try: b_k = c*a_k + d*b_{k-1}
# k=1: -4 = 3c + (-2)d
# k=2: -84 = 55c + (-4)d
A_m3 = Matrix([[3, -2], [55, -4]])
b_m3 = Matrix([-4, -84])
sol_m3 = A_m3.solve(b_m3)
print(f"\n  b_k = {sol_m3[0]}*a_k + {sol_m3[1]}*b_{{k-1}}")
pred_b3_v3 = sol_m3[0]*a[3] + sol_m3[1]*b[2]
print(f"  Predicts b_3 = {pred_b3_v3}, actual = {b[3]}")
test("b_k = c*a_k + d*b_{k-1} (mixed)", pred_b3_v3 == b[3])

# ============================================================
print("\n" + "="*70)
print("ROUND 3: b_k = f(p_k, a_{k-1}, b_{k-1}) linear with p-dependent coeffs?")
print("="*70)
# b_k = alpha(p_k)*a_{k-1} + beta(p_k)*b_{k-1}
# Try alpha = p^2*A, beta = p*B (same shape as a-recurrence)
# Already killed above. Try alpha = -p^2, beta = p
for k in range(1, 4):
    pred = -p[k-1]**2 * a[k-1] + p[k-1] * b[k-1]
    print(f"  k={k}: -p^2*a_prev + p*b_prev = {pred}, actual = {b[k]}")
test("b_k = -p^2*a_{k-1} + p*b_{k-1}", 
     all(-p[k-1]**2*a[k-1] + p[k-1]*b[k-1] == b[k] for k in range(1,4)))

# Try b_k = -p*a_{k-1} + (stuff)
# Try ALL f(p)*a_{k-1} + g(p)*b_{k-1} where f,g are low-degree polys in p
print("\n  Exhaustive search: b_k = (Ap^2+Bp+C)*a_{k-1} + (Dp^2+Ep+F)*b_{k-1}")
# 6 unknowns, 3 equations (k=1,2,3). Underdetermined. Try with fewer params.
# b_k = (Ap+B)*a_{k-1} + (Cp+D)*b_{k-1}  -- 4 unknowns, 3 equations
# k=1(p=3): (3A+B)*1 + (3C+D)*(-2) = -4 => 3A+B - 6C - 2D = -4
# k=2(p=5): (5A+B)*3 + (5C+D)*(-4) = -84 => 15A+3B - 20C - 4D = -84
# k=3(p=7): (7A+B)*55 + (7C+D)*(-84) = -3372 => 385A+55B - 588C - 84D = -3372
# Underdetermined (3 eq, 4 unkn). Try 1-parameter family.
A_ex = Matrix([
    [3, 1, -6, -2],
    [15, 3, -20, -4],
    [385, 55, -588, -84]
])
b_ex = Matrix([-4, -84, -3372])
# Solve in rationals
print("  System (Ap+B)*a_{k-1} + (Cp+D)*b_{k-1} = b_k:")
try:
    sol_ex = A_ex.solve(b_ex)  # particular solution
    print(f"    Solution: A={sol_ex[0]}, B={sol_ex[1]}, C={sol_ex[2]}, D={sol_ex[3]}")
    # Check if solution is rational with small denominators
    all_int = all(s.is_integer for s in sol_ex)
    print(f"    All integer? {all_int}")
    if not all_int:
        print(f"    Coefficients: {[float(s) for s in sol_ex]}")
except Exception as e:
    # Might need nullspace approach for underdetermined
    print(f"    Solving via augmented nullspace...")
    aug = A_ex.col_insert(4, b_ex)
    rref, pivots = aug.rref()
    print(f"    RREF pivots: {pivots}")
    # If rank < 4, there's a family of solutions
    for i in range(3):
        print(f"    Row {i}: {rref.row(i).tolist()}")

# ============================================================
print("\n" + "="*70)
print("ROUND 4: b_k related to R(p_{k+1})?")
print("="*70)
# R(p) = -((p-1)^2 - 2)/2
def R_prime(p):
    return -((p-1)**2 - 2) // 2

for k in range(4):
    p_next = p[k]  # prime at level k
    R_val = R_prime(p_next)
    print(f"  k={k}: b[k]={b[k]}, R(p_{k})=R({p_next})={R_val}, -2R={-2*R_val}")
test("b_k = -2*R(p_k)", all(b[k] == -2*R_prime(p[k]) for k in range(4)))
test("b_k = R(p_k)", all(b[k] == R_prime(p[k]) for k in range(4)))

# ============================================================
print("\n" + "="*70)
print("ROUND 5: c_k = a_k/p_k factorization and recurrence")
print("="*70)
# c_k = a_{k-1}*p_k + b_{k-1}, and a_k = p_k * c_k
c_seq = []
for k in range(1, 5):
    ck = a[k-1]*p[k-1] + b[k-1]
    c_seq.append(ck)
    print(f"  c_{k} = a_{k-1}*p_{k} + b_{k-1} = {a[k-1]}*{p[k-1]} + {b[k-1]} = {ck}")
    print(f"         a_{k} = p_{k}*c_{k} = {p[k-1]}*{ck} = {p[k-1]*ck} (actual: {a[k]})")
    assert p[k-1]*ck == a[k], f"Factorization failed at k={k}"
    if isprime(ck):
        print(f"         c_{k} = {ck} is PRIME")
    else:
        print(f"         c_{k} = {ck} = {dict(factorint(ck))}")

print(f"\n  c sequence: {c_seq}")
print(f"  All prime? {all(isprime(c) for c in c_seq)}")

# Does c_k satisfy c_k = c_{k-1}*p_k + something?
print("\n  c_k recurrence test:")
for k in range(1, len(c_seq)):
    diff = c_seq[k] - c_seq[k-1]*p[k]
    ratio = Fraction(c_seq[k], c_seq[k-1])
    print(f"    c_{k+1}/c_{k} = {c_seq[k]}/{c_seq[k-1]} = {float(ratio):.4f}")
    print(f"    c_{k+1} - c_{k}*p_{k+1} = {diff}")

# ============================================================
print("\n" + "="*70)
print("ROUND 6: -b_k/2 sequence analysis")
print("="*70)
neg_b_half = [-bk//2 for bk in b]
print(f"  -b_k/2: {neg_b_half}")
print(f"  Factorizations:")
for k, v in enumerate(neg_b_half):
    if v == 0:
        print(f"    k={k}: 0")
    elif v == 1:
        print(f"    k={k}: 1")
    else:
        print(f"    k={k}: {v} = {dict(factorint(v))}")

# Ratios
print(f"  Ratios: ", end="")
for k in range(1, len(neg_b_half)):
    print(f"{neg_b_half[k]}/{neg_b_half[k-1]} = {Fraction(neg_b_half[k], max(1,neg_b_half[k-1]))}", end="  ")
print()

# ============================================================
print("\n" + "="*70)
print("ROUND 7: d_k = a_k - b_k sequence")
print("="*70)
d_seq = [a[k] - b[k] for k in range(4)]
print(f"  d_k = a_k - b_k: {d_seq}")
for dk in d_seq:
    if isprime(dk):
        print(f"    {dk} is PRIME")
    else:
        print(f"    {dk} = {dict(factorint(dk))}")

# a_k + b_k
s_seq = [a[k] + b[k] for k in range(4)]
print(f"\n  s_k = a_k + b_k: {s_seq}")
for sk in s_seq:
    if sk == 0:
        print(f"    0")
    else:
        sign = "-" if sk < 0 else ""
        val = abs(sk)
        if val == 1:
            print(f"    {sk}")
        elif isprime(val):
            print(f"    {sk} ({val} is PRIME)")
        else:
            print(f"    {sk} = {sign}{dict(factorint(val))}")

# (a_k + b_k + 1)
print(f"\n  a_k + b_k + 1: {[a[k]+b[k]+1 for k in range(4)]}")

# ============================================================
print("\n" + "="*70)
print("ROUND 8: OEIS-style analysis of -b_k/2 = {1, 2, 42, 1686}")
print("="*70)
seq = [1, 2, 42, 1686]
# Test: Catalan numbers? C_0=1, C_1=1, C_2=2, C_3=5, C_4=14, C_5=42
print(f"  Catalan(5) = 42 ✓ (matches k=2)")
print(f"  Catalan(?) for 1686: ", end="")
# C_n = binomial(2n,n)/(n+1)
for n in range(20):
    cn = sp.binomial(2*n, n) // (n+1)
    if cn == 1686:
        print(f"C_{n} = 1686 ✓")
        break
    if cn > 1686:
        print(f"NOT a Catalan number (C_{n-1}={sp.binomial(2*(n-1),n-1)//n}, C_{n}={cn})")
        break

# Double factorials, subfactorials, etc.
print(f"\n  Testing known sequences:")
# 1, 2, 42, 1686
# Is it (2k)! / k! / something?
for k in range(4):
    val = seq[k]
    print(f"    k={k}: {val}")
    # Check (2k)! / something
    for d in range(1, 20):
        if math.factorial(d) == val:
            print(f"      = {d}!")
            break

# Product formula test
print(f"\n  Product tests:")
print(f"  Product of (p_i - 1) for first k primes:")
prod_pm1 = 1
for k in range(4):
    prod_pm1 *= (p[k] - 1)
    print(f"    k={k}: product = {prod_pm1}, -b/2 = {seq[k]}, ratio = {Fraction(seq[k], prod_pm1)}")

print(f"\n  Product of (p_i + 1) for first k primes:")
prod_pp1 = 1
for k in range(4):
    prod_pp1 *= (p[k] + 1)
    print(f"    k={k}: product = {prod_pp1}, -b/2 = {seq[k]}, ratio = {Fraction(seq[k], prod_pp1)}")

# ============================================================
print("\n" + "="*70)
print("ROUND 9: COMPUTE TRACES — tr(D_even) for primorials")
print("="*70)

def compute_D_sym_and_traces(m):
    """Compute D_sym matrix and even/odd block traces for modulus m."""
    t0 = time.time()
    
    # Coprime residues
    coprimes = sorted([r for r in range(1, m) if math.gcd(r, m) == 1])
    phi_m = len(coprimes)
    
    # Distance matrix
    D = {}
    for i, a_val in enumerate(coprimes):
        for j, b_val in enumerate(coprimes):
            D[(i,j)] = min((b_val - a_val) % m, (a_val - b_val) % m)
    
    # Palindromic pairs: r < m/2 with gcd(r,m)=1
    pairs = []
    for r in coprimes:
        if r < m/2:
            pairs.append(r)
    n = len(pairs)
    
    if n == 0:
        return {'m': m, 'phi': phi_m, 'block_size': 0, 'tr_even': 0, 'tr_odd': 0, 'time': time.time()-t0}
    
    # Build coprime index
    cop_idx = {r: i for i, r in enumerate(coprimes)}
    
    # Even and odd block diagonals (we only need diag for trace)
    tr_even = sp.Integer(0)
    tr_odd = sp.Integer(0)
    
    for idx, r_i in enumerate(pairs):
        i1 = cop_idx[r_i]
        i2 = cop_idx[m - r_i]
        
        # D_even(idx, idx) = (D[i1,i1] + D[i1,i2] + D[i2,i1] + D[i2,i2]) / 2
        d_even_diag = sp.Rational(D[(i1,i1)] + D[(i1,i2)] + D[(i2,i1)] + D[(i2,i2)], 2)
        # D_odd(idx, idx) = (D[i1,i1] - D[i1,i2] - D[i2,i1] + D[i2,i2]) / 2
        d_odd_diag = sp.Rational(D[(i1,i1)] - D[(i1,i2)] - D[(i2,i1)] + D[(i2,i2)], 2)
        
        tr_even += d_even_diag
        tr_odd += d_odd_diag
    
    elapsed = time.time() - t0
    return {
        'm': m, 'phi': phi_m, 'block_size': n,
        'tr_even': tr_even, 'tr_odd': tr_odd,
        'time': elapsed
    }

# Compute traces for primorials
primorial_list = [6, 30, 210, 2310]
T_seq = []  # T_k = tr(D_even) / (phi/2)

for m_val in primorial_list:
    result = compute_D_sym_and_traces(m_val)
    n = result['block_size']
    tr_e = result['tr_even']
    tr_o = result['tr_odd']
    phi = result['phi']
    T = tr_e / (phi // 2) if phi > 2 else tr_e
    T_seq.append(T)
    print(f"  m={m_val}: phi={phi}, block={n}x{n}, tr(D_even)={tr_e}, tr(D_odd)={tr_o}")
    print(f"    tr_even + tr_odd = {tr_e + tr_o} (should be 0 since diag of D_sym is 0)")
    print(f"    T_k = tr(D_even)/(phi/2) = {tr_e}/{phi//2} = {T}")
    print(f"    Time: {result['time']:.2f}s")

print(f"\n  T sequence (tr_even / (phi/2)): {T_seq}")

# ============================================================
print("\n" + "="*70)
print("ROUND 10: Does {a, b, T} form a closed system?")
print("="*70)
print(f"  a = {a[:4]}")
print(f"  b = {b}")
print(f"  T = {T_seq}")

# Test: b_k = f(T_k, a_k, p_k)?
print(f"\n  Test b_k = x*T_k + y:")
if len(T_seq) >= 2:
    # Use k=0,1 to fit x,y; check k=2,3
    A_bt = Matrix([[T_seq[0], 1], [T_seq[1], 1]])
    b_bt = Matrix([b[0], b[1]])
    try:
        sol_bt = A_bt.solve(b_bt)
        print(f"    b = {sol_bt[0]}*T + {sol_bt[1]}")
        for k in range(min(4, len(T_seq))):
            pred = sol_bt[0]*T_seq[k] + sol_bt[1]
            print(f"    k={k}: pred={pred}, actual={b[k]}, {'✓' if pred==b[k] else '✗'}")
    except:
        print(f"    Linear fit failed")

# Test: T_k = x*a_k + y*b_k + z?
print(f"\n  Test T_k = x*a_k + y*b_k + z:")
if len(T_seq) >= 3:
    A_tab = Matrix([
        [a[0], b[0], 1],
        [a[1], b[1], 1],
        [a[2], b[2], 1]
    ])
    b_tab = Matrix([T_seq[0], T_seq[1], T_seq[2]])
    try:
        sol_tab = A_tab.solve(b_tab)
        print(f"    T = {sol_tab[0]}*a + {sol_tab[1]}*b + {sol_tab[2]}")
        if len(T_seq) >= 4:
            pred_T3 = sol_tab[0]*a[3] + sol_tab[1]*b[3] + sol_tab[2]
            print(f"    Predicts T_3 = {pred_T3}, actual = {T_seq[3]}, {'✓' if pred_T3==T_seq[3] else '✗'}")
    except:
        print(f"    Linear fit over (a,b) failed")

# Test: T_k follows its own recurrence T_k = T_{k-1}*p_k^2 + E*p_k?
print(f"\n  Test T_k recurrence (same as a-recurrence):")
for k in range(1, len(T_seq)):
    if T_seq[k-1] != 0:
        E = sp.Rational(T_seq[k] - T_seq[k-1]*p[k-1]**2, p[k-1])
        print(f"    k={k}: T_{k} = T_{k-1}*p^2 + E*p => E = {E} ({'integer' if E.is_integer else 'NOT integer'})")

# ============================================================
print("\n" + "="*70)
print("ROUND 11: Compute FULL D_even and D_odd for small primorials")
print("="*70)
print("Computing full blocks to find ALL invariants...")

def compute_full_blocks(m):
    """Compute full D_even, D_odd matrices."""
    coprimes = sorted([r for r in range(1, m) if math.gcd(r, m) == 1])
    phi_m = len(coprimes)
    
    D = {}
    for i, a_val in enumerate(coprimes):
        for j, b_val in enumerate(coprimes):
            D[(i,j)] = min((b_val - a_val) % m, (a_val - b_val) % m)
    
    pairs = [r for r in coprimes if r < m/2]
    n = len(pairs)
    cop_idx = {r: i for i, r in enumerate(coprimes)}
    
    D_even = sp.zeros(n, n)
    D_odd = sp.zeros(n, n)
    
    for i_idx, r_i in enumerate(pairs):
        i1 = cop_idx[r_i]
        i2 = cop_idx[m - r_i]
        for j_idx, r_j in enumerate(pairs):
            j1 = cop_idx[r_j]
            j2 = cop_idx[m - r_j]
            D_even[i_idx, j_idx] = sp.Rational(D[(i1,j1)] + D[(i1,j2)] + D[(i2,j1)] + D[(i2,j2)], 2)
            D_odd[i_idx, j_idx] = sp.Rational(D[(i1,j1)] - D[(i1,j2)] - D[(i2,j1)] + D[(i2,j2)], 2)
    
    return D_even, D_odd, n

invariants = {}
for m_val in [6, 30, 210]:
    t0 = time.time()
    D_e, D_o, n = compute_full_blocks(m_val)
    
    det_e = D_e.det()
    det_o = D_o.det()
    tr_e = D_e.trace()
    tr_o = D_o.trace()
    
    # Sum of all elements
    sum_e = sum(D_e[i,j] for i in range(n) for j in range(n))
    sum_o = sum(D_o[i,j] for i in range(n) for j in range(n))
    
    # Sum of squared elements (Frobenius squared)
    frob_e = sum(D_e[i,j]**2 for i in range(n) for j in range(n))
    frob_o = sum(D_o[i,j]**2 for i in range(n) for j in range(n))
    
    # Trace of D^2
    D_e_sq = D_e * D_e
    D_o_sq = D_o * D_o
    tr_e2 = D_e_sq.trace()
    tr_o2 = D_o_sq.trace()
    
    el = time.time() - t0
    
    R_val = sp.Rational(det_e, det_o)
    
    invariants[m_val] = {
        'n': n, 'det_e': det_e, 'det_o': det_o, 'R': R_val,
        'tr_e': tr_e, 'tr_o': tr_o,
        'sum_e': sum_e, 'sum_o': sum_o,
        'frob_e': frob_e, 'frob_o': frob_o,
        'tr_e2': tr_e2, 'tr_o2': tr_o2,
    }
    
    print(f"\n  m={m_val} ({n}x{n}) [{el:.2f}s]:")
    print(f"    det(even)={det_e}, det(odd)={det_o}, R={R_val}")
    print(f"    tr(even)={tr_e}, tr(odd)={tr_o}")
    print(f"    sum(even)={sum_e}, sum(odd)={sum_o}")
    print(f"    ||even||_F^2={frob_e}, ||odd||_F^2={frob_o}")
    print(f"    tr(even^2)={tr_e2}, tr(odd^2)={tr_o2}")

# ============================================================
print("\n" + "="*70)
print("ROUND 12: Invariant ratios — hunting for b_k's source")
print("="*70)

print("  Per-level invariant extraction:")
for m_val in [6, 30, 210]:
    inv = invariants[m_val]
    n = inv['n']
    phi = 2*n  # palindromic, so phi/2 = n
    
    # Normalized invariants
    print(f"\n  m={m_val} (n={n}):")
    print(f"    R = det_e/det_o = {inv['R']}")
    print(f"    tr_e/n = {sp.Rational(inv['tr_e'], n)}")
    print(f"    sum_e/n^2 = {sp.Rational(inv['sum_e'], n**2)}")
    print(f"    tr(e^2)/n = {sp.Rational(inv['tr_e2'], n)}")
    print(f"    frob_e/n^2 = {sp.Rational(inv['frob_e'], n**2)}")
    
    # Ratio of odd to even invariants
    if inv['tr_o'] != 0:
        print(f"    tr_e/tr_o = {sp.Rational(inv['tr_e'], inv['tr_o'])}")
    if inv['frob_o'] != 0:
        print(f"    frob_e/frob_o = {sp.Rational(inv['frob_e'], inv['frob_o'])}")
    if inv['tr_o2'] != 0:
        print(f"    tr(e^2)/tr(o^2) = {sp.Rational(inv['tr_e2'], inv['tr_o2'])}")

# ============================================================
print("\n" + "="*70)
print("ROUND 13: The b_k from the polynomial at SECOND data point")
print("="*70)
# b_k is determined by evaluating the level-k polynomial at q = first available prime > p_k
# We have: -2*R(S_k ∪ {q}) = a_k*q^2 + b_k*q - 1
# So: b_k = (-2*R(S_k ∪ {q}) - a_k*q^2 + 1) / q
# 
# We have computed R at many values. Let's verify b_k from known R values:
# Level 0 (S={3}): polynomial -2R(3*q) = 1*q^2 + (-2)*q - 1. Verified.
# Level 1 (S={3,5}): a_1=3, b_1=-4. Check: -2R(15*q) = 3q^2 - 4q - 1
# Level 2 (S={3,5,7}): a_2=55, b_2=-84. Check: -2R(105*q) = 55q^2 - 84q - 1

# Can we get b_k from the EXPRESSION at q = p_k (the prime that defined THIS level)?
# -2R(m_k) = a_{k-1}*p_k^2 + b_{k-1}*p_k - 1
# This gives b_{k-1} from R(m_k) and a_{k-1}. Already known.

# What about the RATIO of the polynomial at two points?
# At q and q': (-2R(q) + 1) / (-2R(q') + 1) = (a*q^2 + b*q) / (a*q'^2 + b*q')
#                                               = q*(a*q + b) / (q'*(a*q' + b))
# This means b = [(q*result_q' * q') - (q'*result_q * q)] ... getting complicated

# The KEY structural relation:
# From two data points at level k: 
#   -2R(m_k * q1) = a_k * q1^2 + b_k * q1 - 1  ... (A)
#   -2R(m_k * q2) = a_k * q2^2 + b_k * q2 - 1  ... (B)
# (A) - (B): -2(R(q1) - R(q2)) = a_k*(q1^2 - q2^2) + b_k*(q1 - q2)
#                                = (q1 - q2)*(a_k*(q1+q2) + b_k)
# So: -2(R(q1) - R(q2)) / (q1 - q2) = a_k*(q1 + q2) + b_k
# Thus: b_k = -2(R(q1) - R(q2))/(q1 - q2) - a_k*(q1 + q2)

# This is just the standard 2-point formula. But it shows b_k is the 
# "discrete derivative correction" of the ratio function R at level k.

# What if b_k comes from the SECOND DERIVATIVE? 
# The polynomial is exactly quadratic, so second derivative is constant = 2*a_k.
# First derivative at q: 2*a_k*q + b_k.
# The derivative at q=0 is exactly b_k.
# So b_k = d/dq[-2R(S_k ∪ {q})]|_{q=0}

print("  b_k is the derivative at q=0 of the level-k polynomial.")
print("  This means b_k measures the LINEAR SENSITIVITY of -2R to adding a 'zero-size' prime.")
print()

# Physical meaning: as q → 0, the polynomial → -1 (constant term).
# b_k = slope at q=0. It measures how fast R departs from the trivial limit.
# And a_k = curvature (second deriv / 2).

# ============================================================
print("\n" + "="*70)
print("ROUND 14: Cross-level identities — can b_k be written in terms of sub-level data?")
print("="*70)

# The level-k polynomial is -2R(S_k ∪ {q}) = a_k*q^2 + b_k*q - 1
# At q = p_{k+1}, this gives -2R(m_{k+1}) = a_k*p_{k+1}^2 + b_k*p_{k+1} - 1
# And a_{k+1} = -2R(m_{k+1}) + 1 = a_k*p_{k+1}^2 + b_k*p_{k+1}
# So: b_k = (a_{k+1} - a_k*p_{k+1}^2) / p_{k+1}

# This is circular (b_k ↔ a_{k+1}). But what if we express b_k differently?

# From a_k = p_k * c_k where c_k = a_{k-1}*p_k + b_{k-1}:
# a_{k+1} = a_k * p_{k+1}^2 + b_k * p_{k+1} = p_{k+1} * (a_k * p_{k+1} + b_k)
# So: c_{k+1} = a_k * p_{k+1} + b_k
# Thus: b_k = c_{k+1} - a_k * p_{k+1}

print("  Key identity: b_k = c_{k+1} - a_k * p_{k+1}")
print(f"  c sequence: {c_seq}")
for k in range(min(3, len(c_seq)-1)):
    bk_from_c = c_seq[k+1] - a[k+1-1]*p[k+1-1]  # careful with indexing
    # Actually: c_{k+1} is for level k+1, which uses p_{k+1}
    # c_{k+1} = a_k * p_{k+1} + b_k
    bk_check = c_seq[k] - a[k]*p[k]  
    # Wait, let me be careful. c_seq[0] = c_1, c_seq[1] = c_2, etc.
    # c_{k+1} = a_k * p_{k+1} + b_k  (for the NEXT level transition)
    pass

# Let me redo with careful indexing:
# c_1 = a_0*p_1 + b_0 = 1*3 + (-2) = 1  => c_seq[0] = 1
# c_2 = a_1*p_2 + b_1 = 3*5 + (-4) = 11 => c_seq[1] = 11
# c_3 = a_2*p_3 + b_2 = 55*7 + (-84) = 301 => c_seq[2] = 301
# c_4 = a_3*p_4 + b_3 = 2107*11 + (-3372) = 19805 => c_seq[3] = 19805

# And a_k = p_k * c_k:
# a_1 = p_1 * c_1 = 3*1 = 3
# a_2 = p_2 * c_2 = 5*11 = 55
# a_3 = p_3 * c_3 = 7*301 = 2107
# a_4 = p_4 * c_4 = 11*19805 = 217855

# So b_k = c_{k+1} - a_k * p_{k+1}:
# b_0 = c_1 - a_0*p_1 = 1 - 1*3 = -2 ✓
# b_1 = c_2 - a_1*p_2 = 11 - 3*5 = -4 ✓
# b_2 = c_3 - a_2*p_3 = 301 - 55*7 = 301-385 = -84 ✓
# b_3 = c_4 - a_3*p_4 = 19805 - 2107*11 = 19805-23177 = -3372 ✓

print("\n  Verification: b_k = c_{k+1} - a_k * p_{k+1}")
for k in range(4):
    bk_check = c_seq[k] - a[k]*p[k]
    print(f"    b_{k} = c_{k+1} - a_{k}*p_{k+1} = {c_seq[k]} - {a[k]}*{p[k]} = {bk_check} (actual: {b[k]}) {'✓' if bk_check==b[k] else '✗'}")

# So the hierarchy is really about c_k!
# c_1 = 1, c_2 = 11, c_3 = 301, c_4 = 19805
# c_{k+1} = a_k * p_{k+1} + b_k = p_k * c_k * p_{k+1} + (c_k - a_{k-1}*p_k)
# Wait: b_{k-1} = c_k - a_{k-1}*p_k, so:
# c_{k+1} = a_k * p_{k+1} + b_k
# And b_k = c_{k+1} - a_k * p_{k+1}  -- this is circular.
# 
# But: a_k = p_k * c_k. So:
# c_{k+1} = p_k * c_k * p_{k+1} + b_k
#
# And b_k also satisfies: from the polynomial evaluated at some q:
# b_k = (-2R(S_k ∪ {q}) - a_k*q^2 + 1) / q

# ============================================================
print("\n" + "="*70)
print("ROUND 15: c_k recurrence — IS THERE ONE?")
print("="*70)
# c_1=1, c_2=11, c_3=301, c_4=19805
# Does c_{k+1} = f(c_k, p_k, p_{k+1})?
# c_2 = 11 = f(1, 3, 5)
# c_3 = 301 = f(11, 5, 7)
# c_4 = 19805 = f(301, 7, 11)

# Try: c_{k+1} = c_k * p_k * p_{k+1} + something
# c_2 = 1*3*5 + X = 15 + X => X = -4 = b_1
# c_3 = 11*5*7 + X = 385 + X => X = -84 = b_2
# c_4 = 301*7*11 + X = 23177 + X => X = -3372 = b_3
# So: c_{k+1} = c_k * p_k * p_{k+1} + b_k
print("  c_{k+1} = c_k * p_k * p_{k+1} + b_k:")
for k in range(3):
    pred = c_seq[k] * p[k] * p[k+1] + b[k+1]
    print(f"    c_{k+2} = {c_seq[k]}*{p[k]}*{p[k+1]} + {b[k+1]} = {c_seq[k]*p[k]*p[k+1]} + {b[k+1]} = {pred} (actual: {c_seq[k+1]}) {'✓' if pred==c_seq[k+1] else '✗'}")

# That's just rewriting a_{k+1} = p_{k+1} * c_{k+1} and a_{k+1} = a_k*p_{k+1}^2 + b_k*p_{k+1}
# c_{k+1} = a_k*p_{k+1} + b_k = p_k*c_k*p_{k+1} + b_k. Still circular.

# WHAT IF WE SUBSTITUTE? b_k = c_{k+1} - p_k*c_k*p_{k+1}
# Then: c_{k+2} = p_{k+1}*c_{k+1}*p_{k+2} + b_{k+1}
#               = p_{k+1}*c_{k+1}*p_{k+2} + c_{k+2} - p_{k+1}*c_{k+1}*p_{k+2}
# That's trivially c_{k+2} = c_{k+2}. The b_k is exactly the non-product part.

# ============================================================
print("\n" + "="*70)
print("ROUND 16: The polynomial -2R(S_k ∪ {q}) evaluated at MANY points")
print("="*70)
# For S_0 = {3}: -2R(3q) = q^2 - 2q - 1  for prime q > 3
# Check: q=5: 25-10-1=14. R(15)=-7. -2*(-7)=14 ✓
# q=7: 49-14-1=34. R(21)=-17. -2*(-17)=34 ✓
# q=2: 4-4-1=-1. R(6)=-(-1+1)/2... hmm, q=2 is special

# For S_1 = {3,5}: -2R(15q) = 3q^2 - 4q - 1 for prime q > 5
# q=7: 3*49-28-1=147-28-1=118. R(105)=-59... wait, R(105)=R(210)=-1053? No.
# Actually R(m) is det_even/det_odd, not the ratio for the specific m.
# Let me recompute from the polynomial for consistency.

# The polynomial gives -2R(S ∪ {q}) where R is the det ratio for m = 2*prod(S)*q
# Wait, no. S is the set of ODD primes. m = 2*prod(S ∪ {q}).

# Level 0: S_0 = {} (empty), polynomial: -2R({q}) = ? This would be m = 2*q
# R(2q) = R(q) by the factor-of-2 invariance. And R(q) = -((q-1)^2-2)/2.
# So -2R({q}) = (q-1)^2 - 2 = q^2 - 2q - 1.
# THIS IS the level-0 polynomial! a_0=1, b_0=-2. ✓

# Level 1: S_1 = {3}, polynomial: -2R({3,q}) = 3q^2 - 4q - 1 for prime q > 3
# R({3,q}) is the det ratio for m = 2*3*q = 6q.
# q=5: -2R(30) = 3*25-20-1 = 54. R(30) = -27. -2*(-27) = 54 ✓
# q=7: -2R(42) = 3*49-28-1 = 118. R(42) = -59.
# q=11: -2R(66) = 3*121-44-1 = 318. R(66) = -159.

# I should compute R(42) and R(66) to verify level-1 polynomial at MORE points
print("  Verifying level-1 polynomial at additional points...")

def compute_ratio(m):
    """Compute det(D_even)/det(D_odd) for modulus m."""
    coprimes = sorted([r for r in range(1, m) if math.gcd(r, m) == 1])
    phi_m = len(coprimes)
    
    D = {}
    for i, a_val in enumerate(coprimes):
        for j, b_val in enumerate(coprimes):
            D[(i,j)] = min((b_val - a_val) % m, (a_val - b_val) % m)
    
    pairs = [r for r in coprimes if r < m/2]
    n = len(pairs)
    if n == 0:
        return None
    
    cop_idx = {r: i for i, r in enumerate(coprimes)}
    
    D_even = sp.zeros(n, n)
    D_odd = sp.zeros(n, n)
    
    for i_idx, r_i in enumerate(pairs):
        i1 = cop_idx[r_i]
        i2 = cop_idx[m - r_i]
        for j_idx, r_j in enumerate(pairs):
            j1 = cop_idx[r_j]
            j2 = cop_idx[m - r_j]
            D_even[i_idx, j_idx] = sp.Rational(D[(i1,j1)] + D[(i1,j2)] + D[(i2,j1)] + D[(i2,j2)], 2)
            D_odd[i_idx, j_idx] = sp.Rational(D[(i1,j1)] - D[(i1,j2)] - D[(i2,j1)] + D[(i2,j2)], 2)
    
    det_e = D_even.det()
    det_o = D_odd.det()
    if det_o == 0:
        return None
    return sp.Rational(det_e, det_o)

# Quick verifications for small moduli
print("\n  Level 1 (S={3}): -2R(6q) = 3q^2 - 4q - 1")
for q in [5, 7, 11, 13]:
    m_val = 6*q
    t0 = time.time()
    R_val = compute_ratio(m_val)
    elapsed = time.time() - t0
    poly_val = 3*q**2 - 4*q - 1
    match = (-2*R_val == poly_val) if R_val is not None else False
    print(f"    q={q}: m={m_val}, R={R_val}, -2R={-2*R_val if R_val else 'N/A'}, poly={poly_val}, {'✓' if match else '✗'} [{elapsed:.2f}s]")

# ============================================================
print("\n" + "="*70)
print("ROUND 17: Sum/product identities for b_k")
print("="*70)

# Instead of b_k alone, look at the POLYNOMIAL value at special points
# P_k(q) = a_k*q^2 + b_k*q - 1
# P_k(1) = a_k + b_k - 1
# P_k(-1) = a_k - b_k - 1
# P_k(p_k) = a_k*p_k^2 + b_k*p_k - 1 = a_{k+1} - 1 = -2R(m_{k+1})

P_at_1 = [a[k] + b[k] - 1 for k in range(4)]
P_at_m1 = [a[k] - b[k] - 1 for k in range(4)]
P_at_2 = [a[k]*4 + b[k]*2 - 1 for k in range(4)]

print(f"  P_k(1) = a_k + b_k - 1: {P_at_1}")
print(f"  P_k(-1) = a_k - b_k - 1: {P_at_m1}")
print(f"  P_k(2): {P_at_2}")

for name, seq in [("P_k(1)", P_at_1), ("P_k(-1)", P_at_m1), ("P_k(2)", P_at_2)]:
    print(f"\n  {name} = {seq}")
    for v in seq:
        if abs(v) <= 1:
            print(f"    {v}")
        elif isprime(abs(v)):
            print(f"    {v} (PRIME)")
        else:
            print(f"    {v} = {'-' if v < 0 else ''}{dict(factorint(abs(v)))}")
    if all(v != 0 for v in seq[1:]):
        print(f"    Ratios: {[Fraction(seq[k], seq[k-1]) for k in range(1, len(seq)) if seq[k-1] != 0]}")

# ============================================================
print("\n" + "="*70)
print("ROUND 18: The discriminant of the polynomial")
print("="*70)
# P_k(q) = a_k*q^2 + b_k*q - 1
# discriminant = b_k^2 + 4*a_k (since c = -1)
disc = [b[k]**2 + 4*a[k] for k in range(4)]
print(f"  disc_k = b_k^2 + 4*a_k: {disc}")
for k, d in enumerate(disc):
    sqrt_d = sp.sqrt(d)
    is_sq = sqrt_d.is_integer
    print(f"    k={k}: disc={d} = {dict(factorint(d))}, sqrt={float(sqrt_d):.4f}, perfect_square={is_sq}")

# Roots of the polynomial
print(f"\n  Roots of a_k*q^2 + b_k*q - 1 = 0:")
for k in range(4):
    r1 = (-b[k] + sp.sqrt(b[k]**2 + 4*a[k])) / (2*a[k])
    r2 = (-b[k] - sp.sqrt(b[k]**2 + 4*a[k])) / (2*a[k])
    print(f"    k={k}: roots = {float(r1):.6f}, {float(r2):.6f}")
    print(f"           product = {float(r1*r2):.6f} (should be -1/a_k = {-1/a[k]:.6f})")
    print(f"           sum = {float(r1+r2):.6f} (should be -b_k/a_k = {-b[k]/a[k]:.6f})")

# ============================================================
print("\n" + "="*70)
print("ROUND 19: MASTER CHECK — b_k as polynomial in p_{k+1}")
print("="*70)
# We know a_k and b_k determine the polynomial at level k.
# If we evaluate at q=p_{k+1}, we get a_{k+1}.
# What if b_k itself satisfies: b_k = a_{k-1}*p_k*X + b_{k-1}*Y for clean X,Y functions of p?

# More carefully: 
# a_k = a_{k-1}*p_k^2 + b_{k-1}*p_k  (PROVED)
# b_k = ???

# What if b_k = a_{k-1}*f(p_k) + b_{k-1}*g(p_k) for polynomial f,g?
# k=1: b_1 = a_0*f(3) + b_0*g(3) = f(3) - 2*g(3) = -4
# k=2: b_2 = a_1*f(5) + b_1*g(5) = 3*f(5) - 4*g(5) = -84
# k=3: b_3 = a_2*f(7) + b_2*g(7) = 55*f(7) - 84*g(7) = -3372

# Try f(p) = A*p^2 + B*p + C, g(p) = D*p + E (5 unknowns, 3 equations -- underdetermined)
# Try f(p) = A*p + B, g(p) = C*p + D (4 unknowns, 3 equations)
# k=1: (3A+B) - 2(3C+D) = -4 => 3A+B-6C-2D = -4
# k=2: 3(5A+B) - 4(5C+D) = -84 => 15A+3B-20C-4D = -84
# k=3: 55(7A+B) - 84(7C+D) = -3372 => 385A+55B-588C-84D = -3372

# 3 equations in 4 unknowns. Let's parametrize by D:
# From eq1: 3A+B-6C = -4+2D ... (I)
# From eq2: 15A+3B-20C = -84+4D ... (II)
# From eq3: 385A+55B-588C = -3372+84D ... (III)

# (II) - 3*(I): 15A+3B-20C - 9A-3B+18C = -84+4D+12-6D
#               6A-2C = -72-2D => 3A-C = -36-D ... (IV)
# (III) - 55*(I)/1: wait, let me use sympy

A_var, B_var, C_var, D_var = sp.symbols('A B C D')
eq1 = sp.Eq(3*A_var + B_var - 6*C_var - 2*D_var, -4)
eq2 = sp.Eq(15*A_var + 3*B_var - 20*C_var - 4*D_var, -84)
eq3 = sp.Eq(385*A_var + 55*B_var - 588*C_var - 84*D_var, -3372)

sol_fg = sp.solve([eq1, eq2, eq3], [A_var, B_var, C_var, D_var])
print(f"  b_k = (Ap+B)*a_{{k-1}} + (Cp+D)*b_{{k-1}}:")
print(f"  Solution: {sol_fg}")

# If solution has a free parameter, check if any integer value of D_var gives clean coefficients
if D_var in sol_fg:
    # It's underdetermined; sol_fg gives A,B,C in terms of D
    print(f"\n  Free parameter D. Testing integer values:")
    for D_test in range(-10, 11):
        A_val = sol_fg[A_var].subs(D_var, D_test)
        B_val = sol_fg[B_var].subs(D_var, D_test)
        C_val = sol_fg[C_var].subs(D_var, D_test)
        if A_val.is_integer and B_val.is_integer and C_val.is_integer:
            print(f"    D={D_test}: A={A_val}, B={B_val}, C={C_val} => f(p)={A_val}p+{B_val}, g(p)={C_val}p+{D_test}")
            # PREDICT b_4 = (Ap_5+B)*a_3 + (Cp_5+D)*b_3 where p_5=13
            b4_pred = (A_val*13 + B_val)*a[3] + (C_val*13 + D_test)*b[3]
            print(f"      => b_4 = ({A_val}*13+{B_val})*{a[3]} + ({C_val}*13+{D_test})*{b[3]} = {b4_pred}")

# Also try f(p) = A*p^2 + B, g(p) = C (3 unknowns, 3 equations — determined!)
print(f"\n  b_k = (Ap^2+B)*a_{{k-1}} + C*b_{{k-1}}:")
A_v2, B_v2, C_v2 = sp.symbols('A2 B2 C2')
eq1b = sp.Eq((9*A_v2 + B_v2)*1 + C_v2*(-2), -4)
eq2b = sp.Eq((25*A_v2 + B_v2)*3 + C_v2*(-4), -84)
eq3b = sp.Eq((49*A_v2 + B_v2)*55 + C_v2*(-84), -3372)
sol_fg2 = sp.solve([eq1b, eq2b, eq3b], [A_v2, B_v2, C_v2])
print(f"  Solution: {sol_fg2}")
if sol_fg2:
    A_s, B_s, C_s = sol_fg2[A_v2], sol_fg2[B_v2], sol_fg2[C_v2]
    print(f"  A={A_s}, B={B_s}, C={C_s}")
    # Predict b_4
    b4_pred2 = (169*A_s + B_s)*a[3] + C_s*b[3]
    print(f"  Predicts b_4 = ({169*A_s+B_s})*{a[3]} + {C_s}*{b[3]} = {b4_pred2}")

# ============================================================
print("\n" + "="*70)
print("ROUND 20: NUCLEAR OPTION — brute force over ALL recurrence forms")
print("="*70)
# Try: b_k = sum of terms using {a_{k-1}, b_{k-1}, a_k, p_k, p_k^2, 1}
# Each term is a product of one of {a_{k-1}, b_{k-1}, 1} with one of {p_k^2, p_k, 1}
# That gives 9 basis functions. 3 equations. Solve exactly.

# Basis: a_{k-1}*p^2, a_{k-1}*p, a_{k-1}, b_{k-1}*p^2, b_{k-1}*p, b_{k-1}, p^2, p, 1
# 9 unknowns, 3 equations. Too many free params.

# Restrict to 3 basis vectors. Try ALL combinations of 3 from 9.
from itertools import combinations

basis_names = ['a*p^2', 'a*p', 'a', 'b*p^2', 'b*p', 'b', 'p^2', 'p', '1']

def basis_val(name, k):
    """Evaluate basis function at level k (using k-1 data to predict b_k)."""
    ak = a[k-1]
    bk = b[k-1]
    pk = p[k-1]  # prime added at level k
    if name == 'a*p^2': return ak * pk**2
    elif name == 'a*p': return ak * pk
    elif name == 'a': return ak
    elif name == 'b*p^2': return bk * pk**2
    elif name == 'b*p': return bk * pk
    elif name == 'b': return bk
    elif name == 'p^2': return pk**2
    elif name == 'p': return pk
    elif name == '1': return 1

print("  Testing all 3-term recurrences b_k = c1*X + c2*Y + c3*Z:")
solutions_found = 0
for combo in combinations(range(9), 3):
    names = [basis_names[i] for i in combo]
    # Build 3x3 system
    rows = []
    for k in range(1, 4):
        row = [basis_val(names[j], k) for j in range(3)]
        rows.append(row)
    A_sys = sp.Matrix(rows)
    b_sys = sp.Matrix([b[1], b[2], b[3]])
    
    det_sys = A_sys.det()
    if det_sys == 0:
        continue
    
    sol_sys = A_sys.solve(b_sys)
    # Check if all coefficients are "clean" (integer or simple fraction)
    all_clean = all(s.is_rational and s.q <= 12 for s in sol_sys)
    all_integer = all(s.is_integer for s in sol_sys)
    
    if all_integer:
        # Verify at the training data:
        passed = True
        for k in range(1, 4):
            pred = sum(sol_sys[j]*basis_val(names[j], k) for j in range(3))
            if pred != b[k]:
                passed = False
                break
        if passed:
            solutions_found += 1
            c1, c2, c3 = sol_sys
            b4_pred = sum(sol_sys[j]*basis_val(names[j], 4) for j in range(3))
            print(f"\n  *** INTEGER SOLUTION: b_k = {c1}*({names[0]}) + {c2}*({names[1]}) + {c3}*({names[2]})")
            print(f"      Predicts b_4 = {b4_pred}")
            
    elif all_clean and not all_integer:
        # Also report clean rational solutions
        passed = True
        for k in range(1, 4):
            pred = sum(sol_sys[j]*basis_val(names[j], k) for j in range(3))
            if pred != b[k]:
                passed = False
                break
        if passed:
            c1, c2, c3 = sol_sys
            b4_pred = sum(sol_sys[j]*basis_val(names[j], 4) for j in range(3))
            if solutions_found < 20:  # limit output
                print(f"\n  RATIONAL SOLUTION: b_k = ({c1})*({names[0]}) + ({c2})*({names[1]}) + ({c3})*({names[2]})")
                print(f"      Predicts b_4 = {b4_pred}")
            solutions_found += 1

if solutions_found == 0:
    print("  NO 3-term recurrence with integer or small-rational coefficients found!")
else:
    print(f"\n  Total solutions found: {solutions_found}")

# ============================================================
print("\n" + "="*70)
print("SCORECARD")  
print("="*70)
print(f"  Total hypotheses tested: {total}")
print(f"  Killed: {passes}")
print(f"  Survived: {fails}")
print(f"\n  {'TOPOLOGY COLLAPSED' if fails == 0 else f'{fails} hypotheses survive'}")
