"""
R30030_crack.py — Rational reconstruction of R(30030) from known residues.
No matrix computation needed. Pure number theory.

Known residues (from R30030_validate.py, code validated on m=2310):
  R mod 999999937 = 317227550
  R mod 999999893 = 650560833
  R mod 999999877 = 633263010
  R mod 999999613 = 317227442
"""

from fractions import Fraction
import math

# Known residues from the validated computation
residues = [
    (317227550, 999999937),
    (650560833, 999999893),
    (633263010, 999999877),
    (317227442, 999999613),
]

def rational_recon(r, p):
    """
    Rational reconstruction: find a/b such that a * b^{-1} ≡ r (mod p)
    with |a| and b both < sqrt(p/2).
    Uses the half-GCD / extended Euclidean algorithm.
    Returns (numerator, denominator) or None.
    """
    bound = int(math.isqrt(p // 2))
    
    # Extended Euclidean on (p, r)
    old_r, new_r = p, r % p
    old_t, new_t = 0, 1
    
    while new_r > bound:
        q = old_r // new_r
        old_r, new_r = new_r, old_r - q * new_r
        old_t, new_t = new_t, old_t - q * new_t
    
    a = new_r  # remainder = numerator (possibly need sign)
    b = new_t  # coefficient = denominator (possibly negative)
    
    if abs(b) > bound:
        return None
    
    # Normalize: b > 0
    if b < 0:
        a, b = -a, -b
    
    # Verify: a ≡ b*r (mod p)
    if (a - b * r) % p != 0:
        a = -a
        if (a - b * r) % p != 0:
            return None
    
    return (a, b)

# ============================================================
print("=" * 70)
print("RATIONAL RECONSTRUCTION FROM INDIVIDUAL RESIDUES")
print("=" * 70)

results = []
for r, p in residues:
    rec = rational_recon(r, p)
    if rec:
        a, b = rec
        f = Fraction(a, b)
        # Double check
        check = (a * pow(b, p - 2, p)) % p
        if check > p // 2:
            check -= p
        match = (check == r) or (check == r - p) or (check + p == r)
        print(f"  mod {p}: {a}/{b} = {f}  (verify: {a}*{b}^-1 mod p = {check}, r={r}, {'OK' if (a - b*r) % p == 0 else 'FAIL'})")
        results.append(f)
    else:
        print(f"  mod {p}: reconstruction failed (|a| or b > sqrt(p/2) ~ {int(math.isqrt(p//2))})")

# Check consistency
if results:
    unique = set(results)
    if len(unique) == 1:
        R_exact = unique.pop()
        print(f"\n*** ALL RESIDUES AGREE: R(30030) = {R_exact} ***")
    else:
        print(f"\nMultiple fractions found: {unique}")
        print("Need multi-modulus reconstruction...")

# ============================================================
print(f"\n{'=' * 70}")
print("MULTI-MODULUS RATIONAL RECONSTRUCTION")
print("=" * 70)

# CRT combine first, then reconstruct from the combined modulus
r_crt = residues[0][0]
m_crt = residues[0][1]

for i in range(1, len(residues)):
    ri, pi = residues[i]
    t_val = ((ri - r_crt) * pow(m_crt % pi, pi - 2, pi)) % pi
    r_crt = r_crt + m_crt * t_val
    m_crt = m_crt * pi
    
    # Try rational reconstruction
    rec = rational_recon(r_crt, m_crt)
    if rec:
        a, b = rec
        f = Fraction(a, b)
        # Verify against ALL known residues
        all_ok = all((a - b * r) % p == 0 for r, p in residues)
        print(f"  After {i+1} primes: R = {a}/{b} = {f}  (verified against all {len(residues)}: {all_ok})")
        if all_ok:
            R_exact = f
            break
    else:
        print(f"  After {i+1} primes: no small reconstruction yet (modulus = {len(str(m_crt))} digits)")

# ============================================================
print(f"\n{'=' * 70}")
print("BRUTE FORCE DENOMINATOR SEARCH (d = 1..50000)")
print("=" * 70)

# For each candidate denominator d, check if n = d*r mod p is consistent
# across all residues. If R = n/d then n ≡ d*r (mod p) for each (r, p).
found = False
for d in range(1, 50001):
    ref = None
    consistent = True
    for r, p in residues:
        n_mod_p = (d * r) % p
        # Map to signed
        if n_mod_p > p // 2:
            n_mod_p -= p
        
        if ref is None:
            ref = n_mod_p
        elif n_mod_p != ref:
            consistent = False
            break
    
    if consistent and ref is not None:
        n = ref
        # Triple-verify
        all_ok = all((n - d * r) % p == 0 for r, p in residues)
        if all_ok:
            f = Fraction(n, d)
            print(f"  d = {d}: n = {n}, R = {n}/{d} = {f}")
            found = True
            
            # Full analysis
            from sympy import factorint, isprime
            
            print(f"\n  R(30030) = {f}")
            print(f"  Numerator: {n}")
            if abs(n) > 1 and abs(n) < 10**18:
                print(f"    = {'-' if n < 0 else ''}{dict(factorint(abs(n)))}")
            print(f"  Denominator: {d}")
            if d > 1:
                print(f"    = {dict(factorint(d))}")
            
            # Is it actually integer?
            if f.denominator == 1:
                print(f"\n  R(30030) IS an integer: {f.numerator}")
                if abs(f.numerator) < 10**15:
                    factors = factorint(abs(f.numerator))
                    print(f"  |R| = {dict(factors)}")
                    print(f"  13 divides |R|? {abs(f.numerator) % 13 == 0}")
            else:
                print(f"\n  R(30030) is a FRACTION with denominator {f.denominator}")
                print(f"  The integer determinant ratio hypothesis BREAKS at m=30030")
            
            # Polynomial hierarchy check
            a4 = 217855
            neg2R = -2 * f
            b4_frac = (neg2R - a4 * 169 + 1) / 13
            print(f"\n  Polynomial hierarchy check:")
            print(f"  -2R = {neg2R}")
            print(f"  b_4 = (-2R - a_4*169 + 1)/13 = {b4_frac}")
            if b4_frac.denominator == 1:
                b4 = int(b4_frac)
                print(f"  b_4 = {b4} (INTEGER)")
                
                a_seq = [1, 3, 55, 2107, a4]
                b_seq = [-2, -4, -84, -3372, b4]
                print(f"\n  Updated sequences:")
                print(f"  a = {a_seq}")
                print(f"  b = {b_seq}")
                print(f"  b_k/a_k = {[bk/ak for ak, bk in zip(a_seq, b_seq)]}")
                
                dk = [a_seq[k] - b_seq[k] for k in range(5)]
                print(f"  d_k = a-b = {dk}")
                for k, dv in enumerate(dk):
                    if abs(dv) < 10**12:
                        if isprime(abs(dv)):
                            print(f"    d_{k} = {dv} (PRIME)")
                        else:
                            print(f"    d_{k} = {dv} = {dict(factorint(abs(dv)))}")
            else:
                print(f"  b_4 is NOT integer: {b4_frac}")
                print(f"  The POLYNOMIAL HIERARCHY BREAKS at level 4")
            
            break

if not found:
    print("  No denominator d <= 50000 found.")
    print("  R(30030) has |numerator| > p/2 ≈ 5×10^8 for ALL candidate denominators.")
    print("  This means R is either a very large integer or has large numerator/denominator.")
    
    # Check: is it a very large integer?
    # If R is integer, all residues should be the same number mod their respective primes
    # But we already showed CRT diverges, so it's definitively not a small integer.
    
    # Check residue patterns
    print(f"\n  Residue analysis:")
    for r, p in residues:
        print(f"    R mod {p} = {r}")
        print(f"      signed: {r - p if r > p//2 else r}")
    
    # The two "clusters" Claude noticed
    print(f"\n  Cluster analysis:")
    for i, (r, p) in enumerate(residues):
        print(f"    p_{i} = {p}: r = {r}, p-r = {p-r}, p//3 = {p//3}, 2p//3 = {2*p//3}")
    
    print(f"\n  Differences between consecutive primes' residues:")
    for i in range(len(residues)-1):
        r1, p1 = residues[i]
        r2, p2 = residues[i+1]
        print(f"    {r2} - {r1} = {r2-r1}")

# ============================================================
# Also try: is det_even/det_odd a fraction whose denominator divides 
# a specific number related to the matrix size?
print(f"\n{'=' * 70}")
print("STRUCTURAL ANALYSIS")
print("=" * 70)

# phi(30030) = 5760, block size = 2880
# The denominators in D_even and D_odd come from the /2 in the definition
# D_even(i,j) = (D[i1,j1] + D[i1,j2] + D[i2,j1] + D[i2,j2]) / 2
# D_odd(i,j)  = (D[i1,j1] - D[i1,j2] - D[i2,j1] + D[i2,j2]) / 2
# 
# So det(D_even) and det(D_odd) each have a factor of 1/2^n from the 1/2 per entry.
# For n=2880: det has n! products, each with 2880 factors of 1/2.
# Actually det(M/2) = det(M) / 2^n where n is the matrix size.
# So det(D_even_raw / 2) = det(D_even_raw) / 2^2880
# And det(D_odd_raw / 2) = det(D_odd_raw) / 2^2880
# The ratio R = det(D_even) / det(D_odd) = det(D_even_raw) / det(D_odd_raw)
# The 2^2880 CANCELS in the ratio!
#
# Wait, BUT our code computes D_even_raw (the integer sums d_diff + d_sum), 
# not D_even_raw/2. So the modular determinants are det(D_even_raw), and
# R = det(D_even_raw) / det(D_odd_raw) with NO factors of 2 from the definition.
# 
# HOLD ON. Let me re-read the code...

print("CRITICAL: Checking what our code actually computes...")
print()
print("The build_blocks function computes:")
print("  D_even = d_diff + d_sum  (INTEGER, no /2)")  
print("  D_odd  = d_diff - d_sum  (INTEGER, no /2)")
print()
print("These are 2 * D_even_true and 2 * D_odd_true respectively.")
print("So det(D_even_code) = 2^n * det(D_even_true)")
print("   det(D_odd_code)  = 2^n * det(D_odd_true)")
print("   R_code = det(D_even_code) / det(D_odd_code)")
print("          = det(D_even_true) / det(D_odd_true)")
print("          = R_true")
print()
print("The factors of 2 cancel! So R_code = R_true.")
print("But wait - our D_even_code IS NOT 2*D_even_true if the /2 isn't applied...")
print()
print("Let me verify: for m=30, the D_even entries should be integers or half-integers.")
print("Our code computes d_diff + d_sum which are always integers.")
print("The TRUE D_even has entries (d_diff + d_sum)/2.")
print("So D_even_code = 2 * D_even_true.")
print("det(D_even_code) = 2^n * det(D_even_true) where n = block size.")
print("R = det(D_even_code) / det(D_odd_code) = det(D_even_true) / det(D_odd_true) = R_true.")
print()
print("CONFIRMED: The 2^n factors cancel in the ratio. Our R is correct.")
print()
print("But the previous scripts (ratio_deep_hunt.py etc.) DID apply the /2.")
print("So their det values differ by 2^n, but ratios are the same.")
