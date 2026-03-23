"""
ratio_pattern_hunter.py — Hunt the pattern in det(D_even)/det(D_odd)
====================================================================

Known primorial ratios:
  m=6    (2·3):       ratio = ???  (to be computed — 1×1 blocks!)
  m=30   (2·3·5):     ratio = -27 = -3³
  m=210  (2·3·5·7):   ratio = -1053 = -3⁴·13
  m=2310 (2·3·5·7·11): ratio = -108927 = -3²·7²·13·19

Strategy:
  1. Compute ratios for ALL squarefree m with manageable φ(m)
  2. Factor each ratio
  3. Look for patterns specific to primorials vs non-primorials
  4. Test algebraic formulas
"""

from math import gcd
from sympy import Matrix, factorint, isprime, totient, primefactors
from sympy import sqrt as sym_sqrt
from itertools import combinations
import time

# ═══════════════════════════════════════════════════════════════════════
#  CORE: Build D_sym and palindromic blocks for any m
# ═══════════════════════════════════════════════════════════════════════

def compute_ratio(m, verbose=False):
    """Compute det(D_even)/det(D_odd) for the symmetric distance matrix on (Z/mZ)*."""
    coprimes = sorted([x for x in range(1, m) if gcd(x, m) == 1])
    phi = len(coprimes)
    
    if phi < 2 or phi % 2 != 0:
        return None  # Need even φ for palindromic decomposition
    
    # Build pairs (r, m-r) with r < m/2
    pairs = [(r, m - r) for r in coprimes if r < m - r]
    half = len(pairs)
    
    if half == 0:
        return None
    
    idx = {x: i for i, x in enumerate(coprimes)}
    
    # Build D_sym
    D = [[0] * phi for _ in range(phi)]
    for i in range(phi):
        for j in range(i + 1, phi):
            d = (coprimes[j] - coprimes[i]) % m
            d = min(d, m - d)
            D[i][j] = d
            D[j][i] = d
    
    # Build palindromic blocks (exact integers)
    D_even = []
    D_odd = []
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
    
    # Exact determinants
    Me = Matrix(D_even)
    Mo = Matrix(D_odd)
    det_e = Me.det()
    det_o = Mo.det()
    
    if det_o == 0:
        return {
            'm': m, 'phi': phi, 'half': half,
            'det_even': int(det_e), 'det_odd': 0,
            'ratio': None, 'factors': None,
            'tr_even': sum(D_even[i][i] for i in range(half)),
            'tr_odd': sum(D_odd[i][i] for i in range(half)),
        }
    
    # Check if ratio is integer
    is_int = (det_e % det_o == 0)
    ratio = det_e // det_o if is_int else float(det_e) / float(det_o)
    
    result = {
        'm': m, 'phi': phi, 'half': half,
        'det_even': int(det_e), 'det_odd': int(det_o),
        'ratio': int(ratio) if is_int else ratio,
        'is_integer': is_int,
        'factors': dict(factorint(abs(int(ratio)))) if is_int and ratio != 0 else None,
        'sign': '+' if ratio > 0 else '-',
        'tr_even': sum(D_even[i][i] for i in range(half)),
        'tr_odd': sum(D_odd[i][i] for i in range(half)),
        'primes_of_m': sorted(primefactors(m)),
        'is_primorial': False,  # set below
    }
    
    if verbose:
        print(f"  m={m}: φ={phi}, blocks={half}×{half}, "
              f"det_e={det_e}, det_o={det_o}, ratio={ratio}")
    
    return result


# ═══════════════════════════════════════════════════════════════════════
#  Generate squarefree numbers with φ(m) ≤ threshold
# ═══════════════════════════════════════════════════════════════════════

def squarefree_numbers(max_m):
    """Generate squarefree numbers up to max_m."""
    result = []
    for n in range(2, max_m + 1):
        factors = factorint(n)
        if all(v == 1 for v in factors.values()):
            result.append(n)
    return result


def is_primorial(m):
    """Check if m is a primorial (product of first k primes)."""
    primes = sorted(primefactors(m))
    if len(primes) == 0:
        return False
    # Check: primes should be {2, 3, 5, 7, ...} consecutive
    expected = []
    p = 2
    for _ in range(len(primes)):
        expected.append(p)
        p = next_prime(p)
    return primes == expected


def next_prime(n):
    """Next prime after n."""
    n += 1
    while not isprime(n):
        n += 1
    return n


# ═══════════════════════════════════════════════════════════════════════
#  MAIN: Compute everything
# ═══════════════════════════════════════════════════════════════════════

print("=" * 80)
print("  RATIO PATTERN HUNTER")
print("  det(D_even)/det(D_odd) for symmetric distance matrix on (Z/mZ)*")
print("=" * 80)

# Phase 1: Compute ratios for all squarefree m with manageable blocks
MAX_BLOCK = 24  # max block size (24×24 for m=210)
MAX_M = 250     # search up to here

print(f"\nPhase 1: Computing ratios for squarefree m ≤ {MAX_M} with block size ≤ {MAX_BLOCK}...")
print("-" * 80)

candidates = squarefree_numbers(MAX_M)
results = []

for m in candidates:
    phi = int(totient(m))
    if phi < 2 or phi % 2 != 0:
        continue
    half = phi // 2
    if half > MAX_BLOCK:
        continue
    
    t0 = time.time()
    r = compute_ratio(m)
    dt = time.time() - t0
    
    if r is None:
        continue
    
    r['is_primorial'] = is_primorial(m)
    r['time'] = dt
    results.append(r)

print(f"\nComputed {len(results)} ratios.")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 2: Display ALL results sorted by m
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 2: ALL RATIOS")
print("=" * 80)
print(f"{'m':>6} {'φ':>4} {'blk':>4} {'ratio':>12} {'int?':>5} {'sign':>5} "
      f"{'factorization':>30} {'primes(m)':>20} {'prim?':>5}")
print("-" * 100)

for r in sorted(results, key=lambda x: x['m']):
    fac_str = ""
    if r['factors']:
        parts = []
        for p, e in sorted(r['factors'].items()):
            parts.append(f"{p}^{e}" if e > 1 else str(p))
        fac_str = " · ".join(parts)
    
    prim_str = "YES" if r['is_primorial'] else ""
    int_str = "Y" if r.get('is_integer', True) else "N"
    primes_str = "·".join(str(p) for p in r['primes_of_m'])
    
    ratio_val = r['ratio'] if r['ratio'] is not None else "DIV0"
    
    print(f"{r['m']:>6} {r['phi']:>4} {r['half']:>4} {ratio_val:>12} {int_str:>5} "
          f"{r['sign']:>5} {fac_str:>30} {primes_str:>20} {prim_str:>5}")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 3: Primorials only — the money shot
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 3: PRIMORIAL RATIOS (the money shot)")
print("=" * 80)

primorial_results = [r for r in results if r['is_primorial']]

# Add m=2310 from saved data
primorial_results.append({
    'm': 2310, 'phi': 480, 'half': 240,
    'det_even': -15850621001836021333221174768983163854584839879704914491212378895303358368015006442035718257372076738197504825297786565812992507478244589568000000000000000,
    'det_odd': 145515996968942698625879485976692315537789894881020449394662286625936254262166464164401096673662881913552239805537530325933813540061184000000000000000,
    'ratio': -108927,
    'is_integer': True,
    'factors': {3: 2, 7: 2, 13: 1, 19: 1},
    'sign': '-',
    'tr_even': 138480,
    'tr_odd': -138480,
    'primes_of_m': [2, 3, 5, 7, 11],
    'is_primorial': True,
})

primorial_results.sort(key=lambda x: x['m'])

print(f"\n{'k':>3} {'m':>6} {'φ':>4} {'blk':>4} {'ratio':>12} "
      f"{'factorization':>30} {'tr_even':>10} {'p_max':>5}")
print("-" * 85)

for i, r in enumerate(primorial_results):
    k = len(r['primes_of_m'])
    p_max = max(r['primes_of_m'])
    fac_str = ""
    if r['factors']:
        parts = []
        for p, e in sorted(r['factors'].items()):
            parts.append(f"{p}^{e}" if e > 1 else str(p))
        fac_str = " · ".join(parts)
    
    print(f"{k:>3} {r['m']:>6} {r['phi']:>4} {r['half']:>4} {r['ratio']:>12} "
          f"{fac_str:>30} {r['tr_even']:>10} {p_max:>5}")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 4: Pattern analysis on primorial ratios
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 4: PATTERN ANALYSIS")
print("=" * 80)

abs_ratios = [abs(r['ratio']) for r in primorial_results]
ms = [r['m'] for r in primorial_results]
phis = [r['phi'] for r in primorial_results]
pmaxes = [max(r['primes_of_m']) for r in primorial_results]
traces = [r['tr_even'] for r in primorial_results]

print("\n--- Basic sequences ---")
print(f"m values:        {ms}")
print(f"|ratio| values:  {abs_ratios}")
print(f"φ(m) values:     {phis}")
print(f"p_max values:    {pmaxes}")
print(f"tr(D_even):      {traces}")

# Test: ratio / trace
print("\n--- Ratio / Trace ---")
for r in primorial_results:
    if r['tr_even'] != 0:
        rat = r['ratio'] / r['tr_even']
        print(f"  m={r['m']:>5}: ratio/tr = {r['ratio']} / {r['tr_even']} = {rat:.6f}")

# Test: ratio / φ
print("\n--- Ratio / φ(m) ---")
for r in primorial_results:
    rat = r['ratio'] / r['phi']
    print(f"  m={r['m']:>5}: ratio/φ = {r['ratio']} / {r['phi']} = {rat:.4f}")

# Test: ratio / (φ/2)
print("\n--- Ratio / (φ/2) ---")
for r in primorial_results:
    rat = r['ratio'] / (r['phi'] // 2)
    print(f"  m={r['m']:>5}: ratio/(φ/2) = {r['ratio']} / {r['phi']//2} = {rat:.4f}")

# Test: ratio / m
print("\n--- Ratio / m ---")
for r in primorial_results:
    rat = r['ratio'] / r['m']
    print(f"  m={r['m']:>5}: ratio/m = {r['ratio']} / {r['m']} = {rat:.6f}")

# Test: ratio * p_max / φ²
print("\n--- ratio * p_max / φ² ---")
for r in primorial_results:
    p = max(r['primes_of_m'])
    rat = r['ratio'] * p / (r['phi'] ** 2)
    print(f"  m={r['m']:>5}: ratio·p_max/φ² = {rat:.6f}")

# Test: ratio / tr²
print("\n--- ratio / tr² ---")
for r in primorial_results:
    if r['tr_even'] != 0:
        rat = r['ratio'] / (r['tr_even'] ** 2)
        print(f"  m={r['m']:>5}: ratio/tr² = {rat:.10f}")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 5: Deep algebraic tests
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 5: ALGEBRAIC FORMULA TESTS")
print("=" * 80)

# For each primorial m = p₁·p₂·...·pₖ, compute various products over primes
print("\n--- Products over prime factors of m ---")
for r in primorial_results:
    primes = r['primes_of_m']
    m = r['m']
    phi = r['phi']
    ratio = r['ratio']
    p_max = max(primes)
    
    # Various products
    prod_p = 1
    for p in primes: prod_p *= p
    prod_p_minus_1 = 1
    for p in primes: prod_p_minus_1 *= (p - 1)
    prod_p_plus_1 = 1
    for p in primes: prod_p_plus_1 *= (p + 1)
    prod_p_minus_2 = 1
    for p in primes[1:]: prod_p_minus_2 *= (p - 2)  # skip p=2
    prod_p_sq_minus_1 = 1
    for p in primes: prod_p_sq_minus_1 *= (p**2 - 1)
    
    print(f"\n  m={m} (primes={primes}):")
    print(f"    ratio = {ratio}")
    print(f"    ∏p = {prod_p} = m")
    print(f"    ∏(p-1) = {prod_p_minus_1} = φ(m)")
    print(f"    ∏(p+1) = {prod_p_plus_1}")
    print(f"    ∏(p-2) [p>2] = {prod_p_minus_2}")
    print(f"    ∏(p²-1) = {prod_p_sq_minus_1}")
    print(f"    ratio / ∏(p-2) = {ratio / prod_p_minus_2 if prod_p_minus_2 != 0 else 'N/A'}")
    print(f"    ratio / ∏(p+1) = {ratio / prod_p_plus_1:.6f}")
    print(f"    ratio / ∏(p²-1) = {ratio / prod_p_sq_minus_1:.6f}")
    
    # Test: is |ratio| related to a product involving only primes NOT in the factorization of m?
    ratio_primes = set(r['factors'].keys()) if r['factors'] else set()
    m_primes = set(primes)
    foreign_primes = ratio_primes - m_primes
    shared_primes = ratio_primes & m_primes
    print(f"    Primes in |ratio|: {sorted(ratio_primes)}")
    print(f"    Primes in m:       {sorted(m_primes)}")
    print(f"    Foreign primes:    {sorted(foreign_primes)}")
    print(f"    Shared primes:     {sorted(shared_primes)}")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 6: Non-primorials with same prime sets
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 6: SAME PRIMES, DIFFERENT PRODUCTS (does ratio depend on product or prime set?)")
print("=" * 80)

# Group results by their prime factor sets
from collections import defaultdict
by_primes = defaultdict(list)
for r in results:
    key = tuple(r['primes_of_m'])
    by_primes[key].append(r)

# Only show groups with multiple members or with primorial
for key in sorted(by_primes.keys()):
    group = by_primes[key]
    if len(group) > 1 or any(r['is_primorial'] for r in group):
        print(f"\n  Prime set = {{{','.join(str(p) for p in key)}}}:")
        for r in sorted(group, key=lambda x: x['m']):
            prim = " ★PRIMORIAL" if r['is_primorial'] else ""
            fac_str = ""
            if r['factors']:
                parts = [f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(r['factors'].items())]
                fac_str = " · ".join(parts)
            print(f"    m={r['m']:>5} φ={r['phi']:>3} ratio={r['ratio']:>10} = "
                  f"{r['sign']}{fac_str}{prim}")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 7: Integer ratio universality check
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 7: IS THE RATIO ALWAYS AN INTEGER?")
print("=" * 80)

non_integer_count = 0
always_negative = True
for r in results:
    if not r.get('is_integer', True):
        print(f"  NON-INTEGER: m={r['m']}, ratio={r['ratio']:.6f}")
        non_integer_count += 1
    if r['ratio'] is not None and r['ratio'] > 0:
        always_negative = False
        print(f"  POSITIVE: m={r['m']}, ratio={r['ratio']}")

if non_integer_count == 0:
    print(f"  ALL {len(results)} ratios are exact integers!")
else:
    print(f"  {non_integer_count} non-integer ratios found.")

if always_negative:
    print(f"  ALL ratios are NEGATIVE!")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 8: Consecutive ratio analysis 
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 8: CONSECUTIVE PRIMORIAL RATIO ANALYSIS")
print("=" * 80)

print("\nk  m       ratio           R(k)/R(k-1)        log₃|ratio|")
print("-" * 70)
prev_ratio = None
for r in primorial_results:
    k = len(r['primes_of_m'])
    ratio = r['ratio']
    import math
    log3 = math.log(abs(ratio)) / math.log(3) if ratio != 0 else 0
    
    if prev_ratio is not None:
        step = ratio / prev_ratio
        print(f"{k}  {r['m']:>5}  {ratio:>12}     {step:>12.4f}           {log3:.4f}")
    else:
        print(f"{k}  {r['m']:>5}  {ratio:>12}     {'---':>12}           {log3:.4f}")
    prev_ratio = ratio

# ═══════════════════════════════════════════════════════════════════════
#  Phase 9: Test specific conjectural formulas for primorial ratios
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 9: CONJECTURAL FORMULA TESTS")
print("=" * 80)

# Test formula: ratio = -∏_{i<j} (p_i + p_j) / something?
print("\n--- Sum-of-pairs products ---")
for r in primorial_results:
    primes = r['primes_of_m']
    ratio = r['ratio']
    
    # ∏_{i<j} (p_i + p_j)
    prod_sum = 1
    for i in range(len(primes)):
        for j in range(i+1, len(primes)):
            prod_sum *= (primes[i] + primes[j])
    
    # ∏_{i<j} (p_j - p_i)
    prod_diff = 1
    for i in range(len(primes)):
        for j in range(i+1, len(primes)):
            prod_diff *= (primes[j] - primes[i])
    
    # ∏_{i<j} (p_i * p_j - 1)
    prod_pq_minus1 = 1
    for i in range(len(primes)):
        for j in range(i+1, len(primes)):
            prod_pq_minus1 *= (primes[i] * primes[j] - 1)
    
    print(f"  m={r['m']:>5}: ratio={ratio:>10}  "
          f"∏(pi+pj)={prod_sum:>10}  ∏(pj-pi)={prod_diff:>10}  "
          f"∏(pipj-1)={prod_pq_minus1:>10}")
    print(f"           ratio/∏(pi+pj)={ratio/prod_sum:.6f}  "
          f"ratio/∏(pj-pi)={ratio/prod_diff:.6f}  "
          f"ratio/∏(pipj-1)={ratio/prod_pq_minus1:.6f}")

# Test: does ratio relate to det(Vandermonde) on the primes?
print("\n--- Vandermonde-related ---")
for r in primorial_results:
    primes = r['primes_of_m']
    ratio = r['ratio']
    
    # Vandermonde determinant = ∏_{i<j} (p_j - p_i)
    vdm = 1
    for i in range(len(primes)):
        for j in range(i+1, len(primes)):
            vdm *= (primes[j] - primes[i])
    
    print(f"  m={r['m']:>5}: ratio={ratio:>10}  Vandermonde={vdm:>10}  ratio/V={ratio/vdm:.6f}")

# Test: relationship between ratio and trace
print("\n--- Trace relationships ---")
for r in primorial_results:
    ratio = r['ratio']
    tr = r['tr_even']
    phi = r['phi']
    half = r['half']
    
    # tr = sum of min(2r, m-2r) over coprime r < m/2
    # For primorial m, tr = (φ/2)(m/2 - 1) ... no, let's compute
    # tr(D_even) = sum of diagonal = sum_k D(r_k, m-r_k) = sum_k min(2r_k, m-2r_k)
    
    if tr != 0:
        print(f"  m={r['m']:>5}: ratio={ratio:>10}  tr={tr:>10}  "
              f"ratio·φ/tr={ratio*phi/tr:.6f}  "
              f"|ratio|/tr={abs(ratio)/tr:.6f}  "
              f"ratio/tr³·φ²={ratio*phi**2/tr**3:.10f}")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 10: The A/B decomposition
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 10: A/B DECOMPOSITION INVARIANTS")
print("  D_even = A + B, D_odd = A - B")
print("  A(i,j) = circular_dist(r_i - r_j, m)")
print("  B(i,j) = circular_dist(r_i + r_j, m)")
print("=" * 80)

for r_data in primorial_results:
    m = r_data['m']
    if m > 210:  # skip 2310 for speed
        continue
    
    coprimes = sorted([x for x in range(1, m) if gcd(x, m) == 1])
    pairs = [(x, m - x) for x in coprimes if x < m - x]
    half = len(pairs)
    rs = [p[0] for p in pairs]  # the smaller residues
    
    # Build A and B matrices
    A = [[0]*half for _ in range(half)]
    B = [[0]*half for _ in range(half)]
    for i in range(half):
        for j in range(half):
            d_diff = (rs[j] - rs[i]) % m
            A[i][j] = min(d_diff, m - d_diff)
            s = (rs[i] + rs[j]) % m
            B[i][j] = min(s, m - s)
    
    Ma = Matrix(A)
    Mb = Matrix(B)
    det_a = Ma.det()
    det_b = Mb.det()
    tr_a = sum(A[i][i] for i in range(half))
    tr_b = sum(B[i][i] for i in range(half))
    
    ratio = r_data['ratio']
    
    print(f"\n  m={m} (blocks {half}×{half}):")
    print(f"    det(A)={det_a}  det(B)={det_b}")
    print(f"    tr(A)={tr_a}  tr(B)={tr_b}")
    print(f"    det(A+B)/det(A-B) = {ratio}")
    if det_a != 0:
        print(f"    det(A+B)/det(A) = {(det_a + det_b) if half == 1 else '...'}")
    if det_b != 0:
        factors_a = dict(factorint(abs(det_a))) if det_a != 0 else {}
        factors_b = dict(factorint(abs(det_b))) if det_b != 0 else {}
        print(f"    |det(A)| factors: {factors_a}")
        print(f"    |det(B)| factors: {factors_b}")

# ═══════════════════════════════════════════════════════════════════════
#  Phase 11: p_max exclusion check for ALL m (not just primorials)
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  Phase 11: p_max EXCLUSION — UNIVERSAL OR PRIMORIAL-SPECIFIC?")
print("=" * 80)

for r in sorted(results, key=lambda x: x['m']):
    if r['ratio'] is None or r['ratio'] == 0 or r['factors'] is None:
        continue
    p_max = max(r['primes_of_m'])
    ratio_primes = set(r['factors'].keys())
    excluded = p_max not in ratio_primes
    status = "EXCLUDED" if excluded else "PRESENT"
    prim = "★" if r['is_primorial'] else " "
    print(f"  {prim} m={r['m']:>5} p_max={p_max:>3} ratio={r['ratio']:>10} "
        f"→ p_max {status} from |ratio|")

# ═══════════════════════════════════════════════════════════════════════
#  FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("  FINAL SUMMARY")
print("=" * 80)

print(f"\nTotal m values tested: {len(results)}")
print(f"Non-integer ratios: {non_integer_count}")
print(f"Always negative: {always_negative}")
print(f"\nPrimorial sequence: {{", end="")
print(", ".join(str(r['ratio']) for r in primorial_results), end="")
print("}")
print(f"Absolute values:   {{", end="")
print(", ".join(str(abs(r['ratio'])) for r in primorial_results), end="")
print("}")
print(f"\nFactorizations:")
for r in primorial_results:
    fac_str = ""
    if r['factors']:
        parts = [f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(r['factors'].items())]
        fac_str = " · ".join(parts)
    print(f"  m={r['m']:>5}: |ratio| = {abs(r['ratio'])} = {fac_str}")

print("\n🐇 The hunt continues...\n")
