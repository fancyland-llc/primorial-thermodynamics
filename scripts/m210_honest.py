"""Honest analysis of m=210 vs m=30 confinement predictions."""
import numpy as np

C_30 = 8 / np.sqrt(5)
C_210 = 48 / np.sqrt(7)

M = 2.5; nu = 0.1; eps = 0.32

# Crossover analysis
K = (C_210/C_30) * M**(1/7 - 1/5) * nu**(-2/47 + 2/7) * eps**(-7/24 + 5/4)
rho_cross = K**(-48/23)

print(f"Ratio constant K = {K:.6f}")
print(f"rho* exponent difference = 23/48 = {23/48:.6f}")
print(f"Crossover rho* = {rho_cross:.6f}")
print(f"  m=210 wins when rho* > {rho_cross:.4f}")
print(f"  ITER: rho* ~ 0.002-0.005")
print(f"  Small tokamaks: rho* ~ 0.01-0.05")
print()

# Ratio scan
print("tau_210/tau_30 vs rho*:")
for rho in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
    ratio = K * rho**(23/48)
    winner = "m=210" if ratio > 1 else "m=30"
    print(f"  rho*={rho:.3f}: ratio={ratio:.3f}x  ({winner} wins)")
print()

# BURN STABILITY
print("=== BURN STABILITY (sensitivity analysis) ===")
print("Sensitivity |d(ln tau)/d(ln x)| = |exponent|")
print()
exponent_data = [
    ("rho*", -5/8, -7/48),
    ("M", 1/5, 1/7),
    ("nu*", -2/7, -2/47),
    ("eps", -5/4, -7/24),
]
for name, e30, e210 in exponent_data:
    reduction = 1 - abs(e210)/abs(e30)
    print(f"  {name:<8s}: m=30 |{e30:+.4f}|={abs(e30):.4f}  "
          f"m=210 |{e210:+.4f}|={abs(e210):.4f}  "
          f"({reduction*100:.0f}% flatter)")

s30 = sum(abs(e) for _, e, _ in exponent_data)
s210 = sum(abs(e) for _, _, e in exponent_data)
print()
print(f"Total sensitivity: m=30={s30:.4f}, m=210={s210:.4f}")
print(f"Sensitivity reduction: {(1 - s210/s30)*100:.0f}%")
print(f"m=210 burn is {s30/s210:.1f}x more robust to parameter perturbations")
print()

# Effective burn quality
print("=== EFFECTIVE BURN QUALITY ===")
print("(tau_E at given params, but NOW consider: flatter scaling = wider operating window)")
print()
print("Absolute tau_E comparison:")
for rho in [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]:
    tau30 = C_30 * rho**(-5/8) * M**(1/5) * nu**(-2/7) * eps**(-5/4)
    tau210 = C_210 * rho**(-7/48) * M**(1/7) * nu**(-2/47) * eps**(-7/24)
    print(f"  rho*={rho:.3f}: tau30={tau30:>10.1f}  tau210={tau210:>10.1f}  "
          f"ratio={tau210/tau30:.3f}x")
print()

# The REAL win: operating window width
# If you perturb rho* by +/-20%, how much does tau_E change?
print("=== OPERATING WINDOW (20% perturbation in each variable) ===")
rho_base = 0.01
delta = 0.20

print(f"\nAt rho*={rho_base}, M={M}, nu*={nu}, eps={eps}:")

for regime, C, exp_rho, exp_M, exp_nu, exp_eps in [
    ("m=30", C_30, -5/8, 1/5, -2/7, -5/4),
    ("m=210", C_210, -7/48, 1/7, -2/47, -7/24),
]:
    tau_base = C * rho_base**exp_rho * M**exp_M * nu**exp_nu * eps**exp_eps
    # Perturb each variable by +20%
    tau_rho_up = C * (rho_base*1.2)**exp_rho * M**exp_M * nu**exp_nu * eps**exp_eps
    tau_rho_dn = C * (rho_base*0.8)**exp_rho * M**exp_M * nu**exp_nu * eps**exp_eps
    swing_rho = (tau_rho_up - tau_rho_dn) / tau_base * 100

    tau_nu_up = C * rho_base**exp_rho * M**exp_M * (nu*1.2)**exp_nu * eps**exp_eps
    tau_nu_dn = C * rho_base**exp_rho * M**exp_M * (nu*0.8)**exp_nu * eps**exp_eps
    swing_nu = (tau_nu_up - tau_nu_dn) / tau_base * 100

    tau_eps_up = C * rho_base**exp_rho * M**exp_M * nu**exp_nu * (eps*1.2)**exp_eps
    tau_eps_dn = C * rho_base**exp_rho * M**exp_M * nu**exp_nu * (eps*0.8)**exp_eps
    swing_eps = (tau_eps_up - tau_eps_dn) / tau_base * 100

    print(f"\n  {regime}: tau_base = {tau_base:.1f}")
    print(f"    rho* +/-20%: swing = {abs(swing_rho):.1f}%")
    print(f"    nu*  +/-20%: swing = {abs(swing_nu):.1f}%")
    print(f"    eps  +/-20%: swing = {abs(swing_eps):.1f}%")
    total_swing = abs(swing_rho) + abs(swing_nu) + abs(swing_eps)
    print(f"    Total swing = {total_swing:.1f}%")

print()
print("=== THE PREDICTION ===")
print(f"The m=210 regime does NOT give higher absolute tau_E at ITER-like rho*.")
print(f"What it gives is {s30/s210:.1f}x flatter scaling = {s30/s210:.1f}x wider burn window.")
print(f"")
print(f"For compact high-rotation devices (rho* ~ 0.02-0.05),")
print(f"the flatter scaling means:")
print(f"  1. More robust ignition (less sensitive to edge perturbations)")
print(f"  2. Wider operating space for sustained burn")
print(f"  3. Compact geometry becomes viable (eps penalty drops from 1.25 to 0.29)")
print(f"")
print(f"The impedance condition Omega*tau_E ~ 7 sets the rotation rate.")
