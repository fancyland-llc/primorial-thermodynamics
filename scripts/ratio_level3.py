"""
ratio_level3.py — Crack the level 3 polynomial
================================================

Need R(3·5·7·q) for q=11 (known: -108927) and q=13 (to compute).
m=1365 = 3·5·7·13, φ=576, blocks=288×288.

Then determine: -2·R(3·5·7·q) = a₃·q² + b₃·q - 1

Known leading coefficients: 1, 3, 55, a₃=?
"""

from math import gcd
from sympy import Matrix, factorint
from fractions import Fraction
import time
import sys

m = 1365  # 3·5·7·13
print(f"Target: m = {m} = 3·5·7·13")
coprimes = sorted([x for x in range(1, m) if gcd(x, m) == 1])
phi = len(coprimes)
print(f"φ({m}) = {phi}")
pairs = [(x, m - x) for x in coprimes if x < m - x]
half = len(pairs)
print(f"Block size: {half}×{half}")
print(f"Building D_sym...")

idx = {x: i for i, x in enumerate(coprimes)}

# Build distance matrix (full)
t0 = time.time()
D = [[0]*phi for _ in range(phi)]
for i in range(phi):
    for j in range(i+1, phi):
        d = (coprimes[j] - coprimes[i]) % m
        d = min(d, m - d)
        D[i][j] = d
        D[j][i] = d
print(f"D_sym built in {time.time()-t0:.1f}s")

# Build palindromic blocks
print(f"Building {half}×{half} palindromic blocks...")
t0 = time.time()
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
dt = time.time() - t0
print(f"Blocks built in {dt:.1f}s")

# Free the huge D matrix
del D

# Compute det(D_odd) first (negative definite, might be simpler)
print(f"\n{'='*60}")
print(f"Computing det(D_odd) [{half}×{half} exact integer]...")
print(f"  Started at: {time.strftime('%H:%M:%S')}")
sys.stdout.flush()

t0 = time.time()
Mo = Matrix(D_odd)
det_odd = Mo.det()
t1 = time.time()
del Mo

print(f"  Finished in {t1-t0:.1f}s")
print(f"  det(D_odd) has {len(str(abs(det_odd)))} digits")
print(f"  sign = {'+' if det_odd > 0 else '-'}")
sys.stdout.flush()

# Compute det(D_even)
print(f"\nComputing det(D_even) [{half}×{half} exact integer]...")
print(f"  Started at: {time.strftime('%H:%M:%S')}")
sys.stdout.flush()

t0 = time.time()
Me = Matrix(D_even)
det_even = Me.det()
t1 = time.time()
del Me

print(f"  Finished in {t1-t0:.1f}s")
print(f"  det(D_even) has {len(str(abs(det_even)))} digits")
print(f"  sign = {'+' if det_even > 0 else '-'}")
sys.stdout.flush()

# Compute ratio
print(f"\n{'='*60}")
print("RESULTS")
print(f"{'='*60}")

if det_odd != 0 and det_even % det_odd == 0:
    ratio = det_even // det_odd
    print(f"Ratio = det(D_even)/det(D_odd) = {ratio}")
    print(f"Integer: YES")
    print(f"|ratio| = {abs(ratio)}")
    print(f"Factorization of |ratio|: {dict(factorint(abs(ratio)))}")
    
    # p_max exclusion
    print(f"\np_max = 13 exclusion:")
    print(f"  det_even % 13 = {det_even % 13}")
    print(f"  det_odd % 13 = {det_odd % 13}")
    if det_even % 13 != 0 and det_odd % 13 != 0:
        print(f"  13 EXCLUDED ✓")
    else:
        print(f"  13 present!")
    
    # Now determine level 3 polynomial
    # -2·R({3,5,7,q}) = a·q² + b·q - 1
    # q=11: -2·(-108927) = 217854 = 121a + 11b - 1 → 121a + 11b = 217855
    # q=13: -2·ratio = 121a + 13b - 1 (if ratio is for q=13)
    neg2R_11 = 2 * 108927  # = 217854
    neg2R_13 = -2 * int(ratio)
    
    print(f"\n{'='*60}")
    print("LEVEL 3 POLYNOMIAL: -2·R(3·5·7·q) = a·q² + b·q - 1")
    print(f"{'='*60}")
    print(f"  q=11: -2R = {neg2R_11 + 1 - 1} → 121a + 11b = {neg2R_11 + 1}")
    print(f"  q=13: -2R = {neg2R_13} → 169a + 13b = {neg2R_13 + 1}")
    
    # Solve: 121a + 11b = 217855 ... [1]
    #        169a + 13b = neg2R_13 + 1 ... [2]
    # From [1]: 11(11a + b) = 217855, so 11a + b = 217855/11
    val1 = neg2R_11 + 1  # = 217855
    val2 = neg2R_13 + 1
    
    # 13·[1] - 11·[2]: 13·121a - 11·169a + 13·11b - 11·13b = 13·val1 - 11·val2
    # (1573 - 1859)a = 13·val1 - 11·val2
    # -286a = 13·val1 - 11·val2
    # a = (11·val2 - 13·val1) / 286
    
    a_num = 11 * val2 - 13 * val1
    a3 = Fraction(a_num, 286)
    
    # b from [1]: b = (val1 - 121·a) / 11
    b3 = Fraction(val1 - 121 * a3, 11)
    
    print(f"\n  a₃ = {a3} = {float(a3):.4f}")
    print(f"  b₃ = {b3} = {float(b3):.4f}")
    
    if a3.denominator == 1 and b3.denominator == 1:
        a3_int = int(a3)
        b3_int = int(b3)
        print(f"\n  BOTH INTEGER!")
        print(f"  -2·R(3·5·7·q) = {a3_int}·q² + ({b3_int})·q - 1")
        
        # Verify
        check_11 = a3_int * 121 + b3_int * 11 - 1
        check_13 = a3_int * 169 + b3_int * 13 - 1
        print(f"  Verify q=11: {check_11} → R = {-check_11//2} (expect -108927)")
        print(f"  Verify q=13: {check_13} → R = {-check_13//2} (expect {ratio})")
        
        # Factor a₃
        print(f"\n  a₃ = {a3_int} = {dict(factorint(a3_int))}")
        
        # Complete hierarchy
        print(f"\n  {'='*60}")
        print(f"  COMPLETE POLYNOMIAL HIERARCHY")
        print(f"  {'='*60}")
        print(f"  Level 0: -2R(q)       = 1·q² - 2·q - 1       (a=1)")
        print(f"  Level 1: -2R(3·q)     = 3·q² - 4·q - 1       (a=3)")
        print(f"  Level 2: -2R(3·5·q)   = 55·q² - 84·q - 1     (a=55)")
        print(f"  Level 3: -2R(3·5·7·q) = {a3_int}·q² + ({b3_int})·q - 1  (a={a3_int})")
        
        # Leading coefficient sequence: 1, 3, 55, a₃
        print(f"\n  Leading coefficient sequence: 1, 3, 55, {a3_int}")
        print(f"  Ratios: 3/1=3, 55/3={55/3:.4f}, {a3_int}/55={a3_int/55:.4f}")
        
        # Is a₃ = a₂·p₃² + b₂·p₃?  (where p₃=7, the prime just added to get level 3)
        # a₂·49 + b₂·7 = 55·49 + (-84)·7 = 2695 - 588 = 2107
        recurrence_test = 55 * 49 + (-84) * 7
        print(f"\n  Recurrence test: a₂·p₃² + b₂·p₃ = 55·49 + (-84)·7 = {recurrence_test}")
        print(f"  Does this equal a₃ + 1? {recurrence_test} vs {a3_int + 1} → {recurrence_test == a3_int + 1}")
        print(f"  Does this equal a₃? {recurrence_test == a3_int}")
        
        # General recurrence: a_k = a_{k-1}·p_k² + b_{k-1}·p_k
        # Check level 0→1: a₁ = a₀·3² + b₀·3 = 1·9 + (-2)·3 = 9-6 = 3 ✓
        # Check level 1→2: a₂ = a₁·5² + b₁·5 = 3·25 + (-4)·5 = 75-20 = 55 ✓
        print(f"\n  RECURRENCE: a_k = a_{{k-1}}·p_k² + b_{{k-1}}·p_k")
        print(f"    Level 0→1: a₁ = 1·9 + (-2)·3 = {1*9 + (-2)*3}")
        print(f"    Level 1→2: a₂ = 3·25 + (-4)·5 = {3*25 + (-4)*5}")
        print(f"    Level 2→3: a₃ = 55·49 + (-84)·7 = {55*49 + (-84)*7} (actual a₃={a3_int})")
        
        # Similarly for b_k:
        # b₀=-2, b₁=-4, b₂=-84
        # If b_k = a_{k-1}·(-2p_k) + b_{k-1}·(something)?
        # b₁ = ? from b₀, a₀, p₁=3
        # -4 = ? 
        
        # Actually, -2R(S∪{q}) evaluated AT q = p_{k+1} gives -2·R(S∪{p_{k+1}})
        # which is the NEXT primorial ratio!
        
        # Predicted next level: predict R(3·5·7·11·q) for q=13
        # Need a₄, b₄ from the level 4 polynomial.
        # But we only have one point for level 4 right now.
        
        # PREDICT R(30030) = R(primorial 13):
        # R(30030) = R(3·5·7·11·13)
        # This is level 3 polynomial evaluated at q=??? — NO!
        # Level 3 is R(3·5·7·q), so q=11 gives R(2310) and q=13 gives R(1365·11/13... no)
        # Actually, R(3·5·7·11) = R(2310). That uses level 3 at q=11.
        # R(3·5·7·11·13) requires level 4: R(3·5·7·11·q) at q=13.
        
        # But can I GET the level 4 polynomial?
        # Need: -2R(3·5·7·11·q) = a₄q² + b₄q - 1
        # If a₄ = a₃·11² + b₃·11 (following recurrence):
        a4_pred = a3_int * 121 + b3_int * 11
        print(f"\n  IF recurrence continues: a₄ = a₃·121 + b₃·11 = {a3_int}·121 + ({b3_int})·11 = {a4_pred}")
        
        # And b₄ from... we need the b recurrence too.
        # Check: b₁ = ? from (a₀, b₀, p₁=3)
        # Need: b_k = f(a_{k-1}, b_{k-1}, p_k)
        # b₀=-2, a₀=1, p₁=3 → b₁=-4
        # b₁=-4, a₁=3, p₂=5 → b₂=-84
        # b₂=-84, a₂=55, p₃=7 → b₃=???
        
        # Try linear: b_k = α·a_{k-1}·p_k + β·b_{k-1}·p_k + γ
        # -4 = α·1·3 + β·(-2)·3 + γ = 3α - 6β + γ
        # -84 = α·3·5 + β·(-4)·5 + γ = 15α - 20β + γ
        # Diff: -80 = 12α - 14β
        # b₃: α·55·7 + β·(-84)·7 + γ = 385α - 588β + γ
        # From diff eq: α = (-80 + 14β)/12 = (-40 + 7β)/6
        # Many solutions. Need b₃ to pin down.
        
        print(f"\n  b recurrence test:")
        print(f"    b₀=-2, b₁=-4, b₂=-84, b₃={b3_int}")
        print(f"    Checking b_k = x·a_{{k-1}}·p_k² + y·b_{{k-1}}·p_k + z:")
        # -4 = 9x - 6y + z
        # -84 = 75x - 20y + z  → diff: -80 = 66x - 14y
        # b₃ = 49·55·x + 7·(-84)·y + z = 2695x - 588y + z
        # From first: z = -4 - 9x + 6y
        # Check with 2nd: 75x - 20y + (-4-9x+6y) = 66x - 14y - 4 = -84 → 66x - 14y = -80 ✓
        # For 3rd: 2695x - 588y + (-4-9x+6y) = 2686x - 582y - 4 = b₃
        # So: 2686x - 582y = b₃ + 4
        # And: 66x - 14y = -80
        # From 2nd: y = (66x + 80)/14 = (33x + 40)/7
        # Sub into 1st: 2686x - 582·(33x+40)/7 = b₃+4
        # 2686x - 82·(33x+40) + 582·(33x+40)·(1-1)... wait let me redo
        # 582/7 = 83.14... not integer. So this form doesn't work with integer coefficients x,y,z.
        
        # Try b_k = x·a_{k-1}·p_k + y·b_{k-1} + z
        # -4 = 3x - 2y + z
        # -84 = 15x - 4y + z  → diff: -80 = 12x - 2y → y = 6x + 40
        # b₃ = 385x - 84y + z
        # Sub z = -4 - 3x + 2y: 385x - 84y + (-4-3x+2y) = 382x - 82y - 4
        # = 382x - 82(6x+40) - 4 = 382x - 492x - 3280 - 4 = -110x - 3284
        # So b₃ = -110x - 3284
        
        # From y = 6x+40, and z = -4-3x+2(6x+40) = -4-3x+12x+80 = 9x+76
        # b₃ = -110x - 3284
        # If x=0: b₃ = -3284, y=40, z=76
        # If x=1: b₃ = -3394, y=46, z=85
        # If x=-1: b₃ = -3174, y=34, z=67
        # If x=-28: b₃ = -3284+3080 = -204, y=-128, z=-176
        
        # Without knowing b₃ in advance, I can't solve.
        # But I DO have b₃! It's the value computed above.
        # b₃ = -110x - 3284 → x = -(b₃ + 3284)/110
        x_val = Fraction(-(b3_int + 3284), 110)
        y_val = 6 * x_val + 40
        z_val = 9 * x_val + 76
        print(f"    Assuming b_k = x·a_{{k-1}}·p_k + y·b_{{k-1}} + z:")
        print(f"    x = {x_val} = {float(x_val):.4f}")
        print(f"    y = {y_val} = {float(y_val):.4f}")
        print(f"    z = {z_val} = {float(z_val):.4f}")
        if x_val.denominator == 1 and y_val.denominator == 1 and z_val.denominator == 1:
            print(f"    ALL INTEGER: x={int(x_val)}, y={int(y_val)}, z={int(z_val)}")
            b4_pred = int(x_val) * a3_int * 11 + int(y_val) * b3_int + int(z_val)
            print(f"    Predicts b₄ = {int(x_val)}·{a3_int}·11 + {int(y_val)}·{b3_int} + {int(z_val)} = {b4_pred}")
        
else:
    # Not integer
    print(f"  Ratio is NOT integer: {float(det_even)/float(det_odd):.6f}")
    print(f"  det_even = {det_even}")
    print(f"  det_odd = {det_odd}")

# Save results
with open("det_1365_exact.txt", "w") as f:
    f.write(f"det_even = {det_even}\n")
    f.write(f"det_odd = {det_odd}\n")
    if det_odd != 0 and det_even % det_odd == 0:
        f.write(f"ratio = {det_even // det_odd}\n")

print(f"\nResults saved to det_1365_exact.txt")
