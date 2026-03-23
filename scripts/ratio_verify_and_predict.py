"""
ratio_verify_and_predict.py — Verify quadratic exactness + predictions
======================================================================

Level 2 polynomial: -2R({3,5,q}) = 55q² - 84q - 1
Predicts: R(255) = R(3·5·17) = -7233

Level 3 polynomial: -2R({3,5,7,q}) = 2107q² - 3372q - 1
Predicts: R(3·5·7·17) = -275799

Also verify the a-recurrence: a_k = eval of previous poly at p_k, plus 1.
"""

from math import gcd
from sympy import Matrix, factorint
from fractions import Fraction
import time

PASS = 0
FAIL = 0

def report(name, condition, detail=""):
    global PASS, FAIL
    status = "PASS" if condition else "FAIL"
    if condition: PASS += 1
    else: FAIL += 1
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

def compute_ratio_int(m):
    """Compute exact integer ratio det(D_even)/det(D_odd)."""
    coprimes = sorted([x for x in range(1, m) if gcd(x, m) == 1])
    phi = len(coprimes)
    if phi < 2 or phi % 2 != 0:
        return None
    pairs = [(r, m - r) for r in coprimes if r < m - r]
    half = len(pairs)
    idx = {x: i for i, x in enumerate(coprimes)}
    
    D = [[0]*phi for _ in range(phi)]
    for i in range(phi):
        for j in range(i+1, phi):
            d = (coprimes[j] - coprimes[i]) % m
            D[i][j] = min(d, m - d)
            D[j][i] = D[i][j]
    
    D_even, D_odd = [], []
    for k in range(half):
        row_e, row_o = [], []
        xk, mxk = pairs[k]
        ik, jk = idx[xk], idx[mxk]
        for l in range(half):
            xl, mxl = pairs[l]
            il, jl = idx[xl], idx[mxl]
            se = D[ik][il] + D[ik][jl] + D[jk][il] + D[jk][jl]
            so = D[ik][il] - D[ik][jl] - D[jk][il] + D[jk][jl]
            row_e.append(se // 2)
            row_o.append(so // 2)
        D_even.append(row_e)
        D_odd.append(row_o)
    
    det_e = Matrix(D_even).det()
    det_o = Matrix(D_odd).det()
    
    if det_o == 0:
        return None
    if det_e % det_o != 0:
        return Fraction(int(det_e), int(det_o))
    return int(det_e // det_o)


# ═══════════════════════════════════════════════════════════════════════
#  STAGE 1: Verify Level 2 polynomial at 4th point (q=17)
# ═══════════════════════════════════════════════════════════════════════

print("=" * 80)
print("  STAGE 1: Verify Level 2 quadratic at 4th point")
print("  -2R({3,5,q}) = 55q² - 84q - 1")
print("=" * 80)

# Predict R(3·5·17 = 255)
q = 17
pred = -(55*q**2 - 84*q - 1) // 2
print(f"\n  Prediction: R(255) = {pred}")

t0 = time.time()
actual = compute_ratio_int(255)
dt = time.time() - t0
print(f"  Computed:   R(255) = {actual} ({dt:.1f}s)")
report(f"Level 2 poly predicts R(255) = {pred}", actual == pred)

# 5th point: q=19
q = 19
pred = -(55*q**2 - 84*q - 1) // 2
print(f"\n  Prediction: R(285) = {pred}")
t0 = time.time()
actual = compute_ratio_int(285)
dt = time.time() - t0
print(f"  Computed:   R(285) = {actual} ({dt:.1f}s)")
report(f"Level 2 poly predicts R(285) = {pred}", actual == pred)

# 6th point: q=23
q = 23
pred = -(55*q**2 - 84*q - 1) // 2
print(f"\n  Prediction: R(345) = {pred}")
t0 = time.time()
actual = compute_ratio_int(345)
dt = time.time() - t0
print(f"  Computed:   R(345) = {actual} ({dt:.1f}s)")
report(f"Level 2 poly predicts R(345) = {pred}", actual == pred)

# ═══════════════════════════════════════════════════════════════════════
#  STAGE 2: Verify Level 3 polynomial at 3rd point (q=17)
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  STAGE 2: Verify Level 3 quadratic at 3rd point")
print("  -2R({3,5,7,q}) = 2107q² - 3372q - 1")
print("=" * 80)

# q=17: m = 3·5·7·17 = 1785, φ=768, blocks=384×384
q = 17
m = 3*5*7*17
pred_level3 = -(2107*q**2 - 3372*q - 1) // 2
print(f"\n  Prediction: R({m}) = {pred_level3}")
print(f"  φ({m}) = 768, blocks = 384×384 — computing...")
t0 = time.time()
actual_level3 = compute_ratio_int(m)
dt = time.time() - t0
print(f"  Computed:   R({m}) = {actual_level3} ({dt:.1f}s)")
report(f"Level 3 poly predicts R({m}) = {pred_level3}", actual_level3 == pred_level3)

if actual_level3 is not None:
    print(f"  |R| = {abs(actual_level3)}")
    if isinstance(actual_level3, int):
        print(f"  Factorization: {dict(factorint(abs(actual_level3)))}")
    # p_max check
    det_check_msg = f"  p_max=17: {'EXCLUDED' if abs(actual_level3) % 17 != 0 else 'PRESENT'}"
    print(det_check_msg)

# ═══════════════════════════════════════════════════════════════════════
#  STAGE 3: Complete hierarchy summary
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  STAGE 3: THE COMPLETE POLYNOMIAL HIERARCHY")
print("=" * 80)

hierarchy = [
    {"level": 0, "S": "∅", "a": 1, "b": -2, "primes": []},
    {"level": 1, "S": "{3}", "a": 3, "b": -4, "primes": [3]},
    {"level": 2, "S": "{3,5}", "a": 55, "b": -84, "primes": [3, 5]},
    {"level": 3, "S": "{3,5,7}", "a": 2107, "b": -3372, "primes": [3, 5, 7]},
]

# a_4 = a_3·11² + b_3·11 = 2107·121 + (-3372)·11
a4 = 2107*121 + (-3372)*11
hierarchy.append({
    "level": 4, "S": "{3,5,7,11}", "a": a4, "b": "?",
    "primes": [3, 5, 7, 11]
})

print(f"\n  Level  S             a_k        b_k        polynomial -2R(S∪{{q}})")
print(f"  " + "-" * 75)

for h in hierarchy:
    b_str = str(h['b']) if h['b'] != '?' else '?'
    if h['b'] != '?':
        poly_str = f"{h['a']}q² + ({h['b']})q - 1"
    else:
        poly_str = f"{h['a']}q² + ?q - 1"
    print(f"  {h['level']:>5}  {h['S']:<14} {h['a']:>8}   {b_str:>8}   {poly_str}")

# ═══════════════════════════════════════════════════════════════════════
#  STAGE 4: The a-recurrence chain
# ═══════════════════════════════════════════════════════════════════════

print(f"\n  a-RECURRENCE: a_{{k+1}} = a_k · p_{{k+1}}² + b_k · p_{{k+1}}")
print(f"  " + "-" * 50)
primes_chain = [3, 5, 7, 11]
for i in range(4):
    a_prev = hierarchy[i]['a']
    b_prev = hierarchy[i]['b']
    p_next = primes_chain[i]
    a_next = a_prev * p_next**2 + b_prev * p_next
    a_actual = hierarchy[i+1]['a']
    ok = "✓" if a_next == a_actual else "✗"
    print(f"  a_{i+1} = {a_prev}·{p_next}² + ({b_prev})·{p_next} = {a_next} {ok}")

# ═══════════════════════════════════════════════════════════════════════
#  STAGE 5: Primorial ratios from the chain
# ═══════════════════════════════════════════════════════════════════════

print(f"\n" + "=" * 80)
print(f"  STAGE 5: PRIMORIAL RATIOS FROM THE HIERARCHY")
print(f"=" * 80)

# R_k = -(a_k - 1)/2
primorials = [
    (6, "{3}", 3),
    (30, "{3,5}", 55),
    (210, "{3,5,7}", 2107),
    (2310, "{3,5,7,11}", a4),
]

print(f"\n  primorial#  S               a_k      R_k = -(a_k-1)/2    Verified")
print(f"  " + "-" * 70)
known_ratios = {6: -1, 30: -27, 210: -1053, 2310: -108927}

for m, s, a in primorials:
    R = -(a - 1) // 2
    ok = "✓" if R == known_ratios[m] else "✗"
    print(f"  {m:>10}  {s:<14}  {a:>8}   {R:>12}   {ok}")

# ═══════════════════════════════════════════════════════════════════════
#  STAGE 6: b coefficient analysis
# ═══════════════════════════════════════════════════════════════════════

print(f"\n" + "=" * 80)
print(f"  STAGE 6: b COEFFICIENT ANALYSIS")
print(f"=" * 80)

b_vals = [-2, -4, -84, -3372]
a_vals = [1, 3, 55, 2107]
p_vals = [3, 5, 7, 11]  # primes at each level transition

print(f"\n  b sequence: {b_vals}")
print(f"  a sequence: {a_vals}")

# Define c_k = a_k + b_k (polynomial evaluated at q=1, plus 1)
c_vals = [a + b for a, b in zip(a_vals, b_vals)]
print(f"  c_k = a_k + b_k: {c_vals}")

# Define d_k = a_k - b_k (related to poly evaluated at q=-1)
d_vals = [a - b for a, b in zip(a_vals, b_vals)]
print(f"  d_k = a_k - b_k: {d_vals}")

# Key ratio: b_k / a_k
print(f"  b_k/a_k: {[Fraction(b,a) for b,a in zip(b_vals, a_vals)]}")
print(f"  b_k/a_k (float): {[b/a for b,a in zip(b_vals, a_vals)]}")

# Does d satisfy a recurrence? d = a - b
# d₀=3, d₁=7, d₂=139, d₃=5479
# d_{k+1} = d_k · p_{k+1}² + ?_k · p_{k+1} ?

print(f"\n  Testing: does d_k follow a recurrence like a_k?")
for i in range(3):
    # d_{i+1} = d_i · p_{i+1}² + e_i · p_{i+1}?
    p = p_vals[i]
    residual = d_vals[i+1] - d_vals[i] * p**2
    e = Fraction(residual, p)
    print(f"    d_{i+1} = {d_vals[i]}·{p}² + e·{p} → e = {e}")

# Check b_k factorizations
print(f"\n  b_k factorizations:")
for i, b in enumerate(b_vals):
    if abs(b) > 1:
        print(f"    b_{i} = {b} = {dict(factorint(abs(b)))}")
    else:
        print(f"    b_{i} = {b}")

# Check: is -2·b_k a perfect form?
print(f"\n  -2·b_k: {[-2*b for b in b_vals]}")
# 4, 8, 168, 6744
# 4 = 4, 8 = 8, 168 = 8·21, 6744 = 8·843 = 8·3·281
g_vals = [-2*b for b in b_vals]
print(f"  Factor out 4: {[g//4 for g in g_vals]}")
# 1, 2, 42, 1686

seq = [g//4 for g in g_vals]
print(f"  Sequence {seq}")
print(f"  Ratios: {seq[1]/seq[0]:.1f}, {seq[2]/seq[1]:.1f}, {seq[3]/seq[2]:.1f}")

# ═══════════════════════════════════════════════════════════════════════
#  SCORECARD
# ═══════════════════════════════════════════════════════════════════════

print(f"\n" + "=" * 80)
print(f"  SCORECARD: {PASS} PASS / {FAIL} FAIL out of {PASS+FAIL}")
print(f"=" * 80)
