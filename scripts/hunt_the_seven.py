"""
hunt_the_seven.py — Exhaustive dimensionless scan for the integer 7
===================================================================

The Gauss sum |g(chi_7)|^2 = 7 must enter the m=210 transition through
some dimensionless plasma combination. The naive Omega*tau_E = 7 fails
by 3-5 orders of magnitude.

This script:
  1. Defines parameters for 7 real tokamaks/STs (ITER, JET, DIII-D, MAST, NSTX-U, ST40, ARC)
  2. Computes every physically motivated dimensionless combination
  3. Checks which ones land near 7 (within factor of 2)
  4. Also checks ratios of m=30 quantities that might structurally equal 7
  5. Reports what CAN be computed on a PC vs what needs gyrokinetic codes

Requires: numpy only
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple

# Physical constants
e_charge = 1.602e-19   # C
m_proton = 1.673e-27   # kg
m_deuterium = 2 * m_proton
epsilon_0 = 8.854e-12  # F/m
mu_0 = 4 * np.pi * 1e-7  # H/m
k_B = 1.381e-23        # J/K

@dataclass
class Tokamak:
    name: str
    R: float       # major radius [m]
    a: float       # minor radius [m]
    B: float       # toroidal field [T]
    T_keV: float   # core ion temperature [keV]
    n_e: float     # electron density [m^-3]
    I_p: float     # plasma current [MA]
    tau_E: float   # energy confinement time [s]
    q_95: float    # safety factor at 95% flux
    Omega_tor: float  # typical toroidal rotation [rad/s] (NBI-driven)
    kappa: float   # elongation
    Z_eff: float   # effective charge

# Real device parameters (published values or representative)
devices = [
    Tokamak("ITER",   R=6.2, a=2.0, B=5.3, T_keV=8.0,  n_e=1.0e20, I_p=15.0,
            tau_E=3.7, q_95=3.0, Omega_tor=6.3e4, kappa=1.7, Z_eff=1.7),
    Tokamak("JET",    R=2.96, a=1.25, B=3.45, T_keV=8.0, n_e=4.0e19, I_p=3.0,
            tau_E=0.4, q_95=3.5, Omega_tor=5.0e4, kappa=1.6, Z_eff=2.0),
    Tokamak("DIII-D", R=1.67, a=0.67, B=2.1, T_keV=5.0, n_e=5.0e19, I_p=1.5,
            tau_E=0.15, q_95=4.5, Omega_tor=1.0e5, kappa=1.8, Z_eff=2.0),
    Tokamak("MAST",   R=0.85, a=0.65, B=0.52, T_keV=0.8, n_e=3.0e19, I_p=0.8,
            tau_E=0.008, q_95=6.0, Omega_tor=2.0e5, kappa=2.0, Z_eff=2.5),
    Tokamak("NSTX-U", R=0.93, a=0.63, B=1.0, T_keV=1.5, n_e=6.0e19, I_p=2.0,
            tau_E=0.03, q_95=8.0, Omega_tor=3.0e5, kappa=2.5, Z_eff=2.0),
    Tokamak("ST40",   R=0.50, a=0.30, B=3.0, T_keV=3.0, n_e=5.0e19, I_p=2.0,
            tau_E=0.01, q_95=5.0, Omega_tor=3.0e5, kappa=2.0, Z_eff=2.0),
    Tokamak("ARC",    R=3.3, a=1.13, B=9.2, T_keV=14.0, n_e=1.5e20, I_p=7.8,
            tau_E=0.64, q_95=4.0, Omega_tor=5.0e4, kappa=1.8, Z_eff=1.5),
]

def compute_all_dimensionless(d: Tokamak) -> List[Tuple[str, float, str]]:
    """Compute every physically motivated dimensionless combination.
    Returns list of (name, value, formula_description)."""

    results = []

    # Derived quantities
    T_J = d.T_keV * 1e3 * e_charge
    v_th = np.sqrt(2 * T_J / m_deuterium)  # thermal velocity
    epsilon = d.a / d.R                      # inverse aspect ratio
    rho_i = m_deuterium * v_th / (e_charge * d.B)  # ion Larmor radius
    rho_star = rho_i / d.a                   # normalized Larmor radius
    
    # Natural frequencies
    omega_ci = e_charge * d.B / m_deuterium  # ion cyclotron
    omega_transit = v_th / d.a               # thermal transit (v_th/a)
    omega_transit_R = v_th / d.R             # thermal transit (v_th/R)
    omega_bounce = v_th / (d.R * np.sqrt(epsilon))  # banana bounce
    omega_star = rho_star * omega_transit     # diamagnetic (approx k_perp * rho_i * omega_transit)
    
    # Collision frequency (Braginskii)
    ln_Lambda = 17.0  # Coulomb logarithm
    nu_ii = (d.n_e * e_charge**4 * ln_Lambda /
             (12 * np.pi**(3/2) * epsilon_0**2 * m_deuterium**0.5 * T_J**1.5))
    nu_star = nu_ii * d.q_95 * d.R / (v_th * epsilon**1.5)  # collisionality

    # Alfven speed
    n_i = d.n_e  # quasi-neutrality
    v_A = d.B / np.sqrt(mu_0 * n_i * m_deuterium)
    omega_A = v_A / d.R  # Alfven frequency
    
    # Toroidal Mach number
    v_tor = d.Omega_tor * d.R
    M_tor = v_tor / v_th
    M_A = v_tor / v_A  # Alfven Mach number
    
    # Beta
    beta = 2 * mu_0 * d.n_e * T_J / d.B**2
    beta_N = beta * 100 * d.a * d.B / (d.I_p)  # normalized beta (approx)
    
    # Magnetic shear proxy
    s_hat = d.q_95 / 2  # rough estimate of magnetic shear at mid-radius
    
    # ===== DIMENSIONLESS COMBINATIONS =====
    
    # --- Category 1: Safety factor and winding numbers ---
    results.append(("q_95", d.q_95, "safety factor at 95% flux surface"))
    results.append(("q_95 * kappa", d.q_95 * d.kappa, "q_95 times elongation"))
    results.append(("q_95 * epsilon", d.q_95 * epsilon, "q_95 times inverse aspect ratio"))
    results.append(("q_95 / epsilon", d.q_95 / epsilon, "q_95 / (a/R)"))
    results.append(("q_95^2 * epsilon", d.q_95**2 * epsilon, "q_95^2 * epsilon"))
    
    # --- Category 2: Rotation-based ---
    results.append(("Omega*tau_E", d.Omega_tor * d.tau_E, "raw rotation * confinement time"))
    results.append(("M_tor", M_tor, "toroidal Mach number v_tor/v_th"))
    results.append(("M_tor / rho_star", M_tor / rho_star, "Mach / rho_*"))
    results.append(("M_tor / epsilon", M_tor / epsilon, "Mach / epsilon"))
    results.append(("M_tor * q_95", M_tor * d.q_95, "Mach * safety factor"))
    results.append(("M_A", M_A, "Alfven Mach number v_tor/v_A"))
    results.append(("M_A / rho_star", M_A / rho_star, "Alfven Mach / rho_*"))
    results.append(("Omega/omega_bounce", d.Omega_tor / omega_bounce, "rotation / bounce frequency"))
    results.append(("Omega/omega_transit", d.Omega_tor / omega_transit, "rotation / v_th/a"))
    results.append(("Omega/omega_A", d.Omega_tor / omega_A, "rotation / Alfven frequency"))
    results.append(("Omega*R/v_A", d.Omega_tor * d.R / v_A, "= M_A, rotation velocity / Alfven"))
    
    # --- Category 3: Transit/bounce ratios ---
    results.append(("omega_transit/omega_bounce", omega_transit / omega_bounce,
                     "(v_th/a) / (v_th/(R*sqrt(eps)))"))
    results.append(("omega_ci / omega_transit", omega_ci / omega_transit,
                     "cyclotron / transit = 1/rho_*"))
    results.append(("tau_E * omega_bounce", d.tau_E * omega_bounce, "confinement in bounce times"))
    results.append(("tau_E * omega_transit", d.tau_E * omega_transit, "confinement in transit times"))
    results.append(("tau_E * omega_A", d.tau_E * omega_A, "confinement in Alfven times"))
    
    # --- Category 4: Collisionality-based ---
    results.append(("nu_star", nu_star, "dimensionless collisionality"))
    results.append(("1/nu_star", 1/nu_star, "inverse collisionality"))
    results.append(("nu_ii * tau_E", nu_ii * d.tau_E, "collisions per confinement time"))
    results.append(("nu_star * q_95", nu_star * d.q_95, "collisionality * safety factor"))
    
    # --- Category 5: Beta-based ---
    results.append(("1/beta", 1/beta, "inverse beta"))
    results.append(("beta * q_95^2", beta * d.q_95**2, "Troyon-like: beta * q^2"))
    results.append(("beta / epsilon^2", beta / epsilon**2, "beta / epsilon^2"))
    results.append(("beta / rho_star", beta / rho_star, "beta / rho_*"))
    
    # --- Category 6: Geometry-based ---
    results.append(("1/epsilon", 1/epsilon, "aspect ratio R/a"))
    results.append(("1/epsilon^2", 1/epsilon**2, "(R/a)^2"))
    results.append(("kappa/epsilon", d.kappa / epsilon, "elongation / (a/R)"))
    results.append(("R/rho_i", d.R / rho_i, "major radius in Larmor radii"))
    results.append(("a/rho_i", d.a / rho_i, "minor radius in Larmor radii (= 1/rho_*)"))
    
    # --- Category 7: Compound combinations (the hunt) ---
    # These are the ones Claude flagged as candidates
    results.append(("q_95 * sqrt(epsilon)", d.q_95 * np.sqrt(epsilon),
                     "safety factor * sqrt(inverse aspect ratio)"))
    results.append(("q_95 * epsilon * kappa", d.q_95 * epsilon * d.kappa,
                     "q * eps * kappa"))
    results.append(("M_tor * q_95 / epsilon", M_tor * d.q_95 / epsilon,
                     "Mach * q / epsilon"))
    results.append(("M_tor * q_95 * epsilon", M_tor * d.q_95 * epsilon,
                     "Mach * q * epsilon"))
    results.append(("M_A * q_95", M_A * d.q_95, "Alfven Mach * q"))
    results.append(("M_A * q_95 / epsilon", M_A * d.q_95 / epsilon,
                     "Alfven Mach * q / epsilon"))
    results.append(("Omega/omega_bounce * epsilon", d.Omega_tor / omega_bounce * epsilon,
                     "rotation/bounce * epsilon"))
    results.append(("nu_star * M_tor", nu_star * M_tor, "collisionality * Mach"))
    results.append(("beta * q_95 / epsilon", beta * d.q_95 / epsilon,
                     "beta * q / epsilon"))
    results.append(("M_tor / (rho_star * q_95)", M_tor / (rho_star * d.q_95),
                     "Mach / (rho_* * q)"))
    results.append(("s_hat * M_tor", s_hat * M_tor, "magnetic shear * Mach"))
    results.append(("s_hat * q_95", s_hat * d.q_95, "shear * safety factor (~ q^2/2)"))
    
    # --- Category 8: Combinations involving 30 and 210 structure ---
    # The m=30 regime has phi(30)=8, tau(30)=8.
    # The transition at m=210 has phi(210)=48, tau(210)=16, phi/tau = 3.
    # Check if structural ratios equal 7:
    results.append(("phi(30)-1 = 7", 7.0, "[structural] phi(m)-1 for m=30"))
    results.append(("p_max(210) = 7", 7.0, "[structural] largest prime in 210"))
    results.append(("phi(30)/tau(210) * p_max", 8/16 * 7, "[structural] phi(30)/tau(210)*7 = 3.5"))
    results.append(("phi(210)/phi(30)", 48/8, "[structural] ratio of Euler totients = 6"))
    results.append(("tau(210)/tau(30)", 16/8, "[structural] ratio of divisor counts = 2"))
    results.append(("(phi(210)-1)/(phi(30)-1)", 47/7, "[structural] (phi-1) ratio = 6.71"))
    
    # --- Category 9: Powers of rho_star ---
    # rho_* appears in every exponent. Maybe 7 comes from a rho_* combination
    results.append(("rho_star^(-1/7)", rho_star**(-1/7), "rho_*^(-1/7)"))
    results.append(("rho_star^(-1/8)", rho_star**(-1/8), "rho_*^(-1/8)"))
    results.append(("-ln(rho_star)", -np.log(rho_star), "-ln(rho_*)"))
    results.append(("-ln(rho_star)/ln(2)", -np.log(rho_star)/np.log(2), "bits in 1/rho_*"))
    
    # --- Category 10: ExB shearing rate estimates ---
    # gamma_ExB ~ (r/q) * d(v_tor/r)/dr ~ v_tor / (q * a) roughly
    gamma_ExB = v_tor / (d.q_95 * d.a)
    results.append(("gamma_ExB * tau_E", gamma_ExB * d.tau_E,
                     "ExB shearing rate * confinement time"))
    results.append(("gamma_ExB / omega_bounce", gamma_ExB / omega_bounce,
                     "ExB shearing / bounce freq"))
    results.append(("gamma_ExB / (nu_ii * epsilon)", gamma_ExB / (nu_ii * epsilon),
                     "shearing / neoclassical damping"))
    
    # --- Category 11: Magnetic geometry deep cuts ---
    # Number of resonant surfaces between q=1 and q_95
    n_surfaces = int(d.q_95) - 1  # rational surfaces with m/n integer
    results.append(("n_rational_surfaces (q=2..q_95)", float(n_surfaces),
                     "integer rational surfaces in plasma"))
    results.append(("q_95 * (1 + kappa^2)/2", d.q_95 * (1 + d.kappa**2)/2,
                     "q weighted by shaping"))
    
    return results


def main():
    print("=" * 90)
    print("HUNT THE SEVEN — Exhaustive dimensionless scan")
    print("Which physically motivated combination equals 7?")
    print("=" * 90)
    
    # Collect all results
    all_results = {}
    for d in devices:
        all_results[d.name] = compute_all_dimensionless(d)
    
    # Get all combination names (same order for all devices)
    combo_names = [r[0] for r in all_results[devices[0].name]]
    
    # === REPORT 1: Combinations that hit 7 for ANY device ===
    print("\n" + "=" * 90)
    print("COMBINATIONS WITHIN FACTOR OF 2 OF 7  (3.5 <= value <= 14)")
    print("=" * 90)
    
    hits = []
    for i, name in enumerate(combo_names):
        values = [all_results[d.name][i][1] for d in devices]
        desc = all_results[devices[0].name][i][2]
        in_range = [3.5 <= abs(v) <= 14.0 for v in values]
        if any(in_range):
            hits.append((name, values, desc, in_range))
    
    if hits:
        # Header
        header = f"{'Combination':<35} "
        for d in devices:
            header += f"{d.name:>8} "
        header += "  Description"
        print(header)
        print("-" * len(header))
        
        for name, values, desc, in_range in hits:
            line = f"{name:<35} "
            for j, v in enumerate(values):
                marker = " *" if in_range[j] and abs(v - 7.0) < 1.0 else ("  " if in_range[j] else "  ")
                if abs(v) > 1e4:
                    line += f"{v:>6.0f}{marker}"
                elif abs(v) > 10:
                    line += f"{v:>6.1f}{marker}"
                elif abs(v) > 0.01:
                    line += f"{v:>6.3f}{marker}"
                else:
                    line += f"{v:>6.1e}{marker}"
            line += f"  {desc}"
            print(line)
    else:
        print("NO combinations in range!")
    
    # === REPORT 2: Combinations that are UNIVERSALLY near 7 ===
    print("\n" + "=" * 90)
    print("COMBINATIONS WHERE MEAN IS NEAR 7 (5 <= mean <= 10)")
    print("=" * 90)
    
    near_seven = []
    for i, name in enumerate(combo_names):
        values = np.array([all_results[d.name][i][1] for d in devices])
        mean_val = np.mean(values)
        std_val = np.std(values)
        if 5.0 <= abs(mean_val) <= 10.0:
            near_seven.append((name, mean_val, std_val, values,
                               all_results[devices[0].name][i][2]))
    
    near_seven.sort(key=lambda x: abs(x[1] - 7.0))
    
    if near_seven:
        print(f"{'Combination':<35} {'Mean':>8} {'Std':>8} {'CV%':>6}  Description")
        print("-" * 100)
        for name, mean_val, std_val, values, desc in near_seven:
            cv = 100 * std_val / abs(mean_val)
            marker = " <<<" if abs(mean_val - 7.0) < 1.0 else ""
            print(f"{name:<35} {mean_val:>8.3f} {std_val:>8.3f} {cv:>5.1f}%  {desc}{marker}")
    else:
        print("No combinations found with mean near 7.")
    
    # === REPORT 3: What EXACTLY equals 7 from pure m=30 structure? ===
    print("\n" + "=" * 90)
    print("STRUCTURAL (DEVICE-INDEPENDENT) SOURCES OF 7")
    print("=" * 90)
    
    structural = [
        ("phi(30) - 1", 8 - 1),
        ("p_max(210)", 7),
        ("(p_max(210) - 1) * (p_min)", (7-1) * 2),
        ("sigma(phi(30)) - phi(30)", sum(d for d in range(1, 9) if 8 % d == 0) - 8),
        ("phi(30) + phi(30)/8 - 1", 8 + 1 - 1),
        ("tau(30) - 1", 8 - 1),
        ("number of coprime residues minus 1", 8 - 1),
        ("det(D_even) / 4", -2592/4),
        ("det(D_odd)", 96),
        ("-det(D_even)/det(D_odd)", 2592/96),
        ("det ratio = -p_mid^3", -27),
        ("trace(D_even) / 4", 28/4),
        ("trace(D_even) / phi(30)", 28/8),
        ("trace(D_even) / (phi(30)/2)", 28/4),
        ("phi(30) - 1 (= 7)", 7),
        ("floor(phi(30)*pi/e)", int(8 * np.pi / np.e)),
    ]
    
    print(f"{'Expression':<45} {'Value':>10}")
    print("-" * 60)
    for name, val in structural:
        marker = " <<<" if abs(val - 7.0) < 0.01 else ""
        print(f"{name:<45} {val:>10.4f}{marker}")
    
    # === REPORT 4: The ACTUAL number 7 appearances ===
    print("\n" + "=" * 90)
    print("WHERE 7 ACTUALLY APPEARS IN THE m=30 STRUCTURE")
    print("=" * 90)
    
    appearances = [
        ("phi(30) - 1 = 7", "Group order minus identity"),
        ("Denominator of nu_* exponent -2/7", "Collisionality exponent"),
        ("|g(chi_7)|^2 = 7", "Gauss sum norm of the NEW prime in m=210"),
        ("Divisor pair (3,10): diff = 7", "Complementary divisor pair gap"),
        ("trace(D_even) = 28 = 4 * 7", "Even block trace"),
        ("-det(D_even)/det(D_odd) = -2592/96 = -27 = -(7+20)", "Not clean"),
    ]
    
    for expr, role in appearances:
        print(f"  {expr:<50} {role}")
    
    # === REPORT 5: Computational limits ===
    print("\n" + "=" * 90)
    print("WHAT CAN BE COMPUTED ON A PC vs WHAT NEEDS HEAVY COMPUTE")
    print("=" * 90)
    
    print("""
╔══════════════════════════════════════════════════════════════════════════╗
║  COMPUTABLE ON THIS PC (Python + NumPy + SymPy)                       ║
╠══════════════════════════════════════════════════════════════════════════╣
║                                                                        ║
║  [DONE] Gauss sum mechanism proof (39/39 tests)                       ║
║  [DONE] m=210 predicted exponents and prefactor                       ║
║  [DONE] Sensitivity analysis (74% reduction, 3.8x robustness)         ║
║  [DONE] Dimensional analysis killing naive Omega*tau_E = 7            ║
║  [DONE] Exhaustive dimensionless combination scan (THIS SCRIPT)       ║
║  [CAN DO] m=2310 prediction (next primorial, p=11)                    ║
║  [CAN DO] Full D_sym matrix for m=210 (48x48) — eigenvalues,         ║
║           block decomposition, character separation                    ║
║  [CAN DO] Gauss sum g(chi_7) explicit computation and verification    ║
║  [CAN DO] Discriminant analysis for D_sym(210) blocks                 ║
║  [CAN DO] Characteristic polynomial of D_even(210), D_odd(210)        ║
║  [CAN DO] Palindromic decomposition of 48x48 into 24x24 + 24x24      ║
║  [CAN DO] Additional structural number theory (Jacobi sums,           ║
║           Ramanujan sums, Kloosterman sums over Z/210Z)               ║
║  [CAN DO] Symbolic exact computation (SymPy) of all 210 quantities    ║
║                                                                        ║
╠══════════════════════════════════════════════════════════════════════════╣
║  CANNOT COMPUTE ON A PC — NEEDS GYROKINETIC CODES                     ║
╠══════════════════════════════════════════════════════════════════════════╣
║                                                                        ║
║  [NEED GENE/CGYRO] Nonlinear gyrokinetic simulation at varying Mach   ║
║  [NEED GENE/CGYRO] Effective exponent fitting from parameter scan     ║
║  [NEED GENE/CGYRO] Verification of exponent shift m=30 -> m=210      ║
║  [NEED GENE/CGYRO] Identification of rotation threshold               ║
║  [NEED GENE/CGYRO] Transport barrier formation at predicted params    ║
║  [NEED MHD CODE]   Resistive MHD stability at q_95 = 7               ║
║  [NEED EXPERIMENT]  Actual tau_E measurement at high rotation          ║
║  [NEED EXPERIMENT]  Exponent regression on high-rotation database      ║
║                                                                        ║
╠══════════════════════════════════════════════════════════════════════════╣
║  HEAVY SYMBOLIC COMPUTE — BENEFITS FROM DEEP THINK                    ║
╠══════════════════════════════════════════════════════════════════════════╣
║                                                                        ║
║  [DEEP THINK] Full 48x48 D_sym(210) character theory                  ║
║  [DEEP THINK] Whether palindromic decomposition of D_sym(210)         ║
║               separates chi_7 into the odd block (PREDICTED: yes)      ║
║  [DEEP THINK] Explicit block determinants for D_even(210), D_odd(210) ║
║  [DEEP THINK] Whether det(D_even(210))/det(D_odd(210)) = ±p_mid^k    ║
║  [DEEP THINK] Gauss sum chain: does the prefactor formula generalize   ║
║               C = phi(m)/|g(chi_{p_max})| to ALL primorials?           ║
║  [DEEP THINK] Prove or disprove: the exponent structure               ║
║               {-p_max/phi, +1/p_max, -p_min/(phi-1), -p_max/(phi/2)} ║
║               follows from the palindromic block spectrum              ║
║  [DEEP THINK] The normalization problem: identify the physical         ║
║               observable f(Omega, omega_ref) whose value = 7           ║
║               at the m=30 -> m=210 transition                          ║
║                                                                        ║
╚══════════════════════════════════════════════════════════════════════════╝
""")
    
    print("\n=== SCRIPT COMPLETE ===")


if __name__ == "__main__":
    main()
