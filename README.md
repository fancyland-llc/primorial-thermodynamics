# Primorial Thermodynamics: Gauss Sums, the IPB98 Scaling Law, and the Hierarchy Break at m = 30030

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19188924.svg)](https://doi.org/10.5281/zenodo.19188924)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Author:** Antonio P. Matos ([ORCID 0009-0002-0722-3752](https://orcid.org/0009-0002-0722-3752))
**Affiliation:** Independent Researcher; Fancyland LLC
**Date:** March 23, 2026
**MSC 2020:** 11L05 (primary), 11A15, 82B05, 82D10 (secondary)

## Summary

We prove that the IPB98(y,2) empirical scaling law for energy confinement time in tokamak plasmas is determined by the Gauss sum structure of the multiplicative group (Z/30Z)\*, where 30 = 2·3·5 is the primorial of the first three primes. The prefactor is C = φ(30)/|g(χ₅)| = 8/√5 = 3.57771…, matching the experimental fit to 0.005% — a factor of 173× closer than the previously assumed √(4π).

### Key results

| Result | Status |
| --- | --- |
| Uniqueness of m = 30 as the primorial matching IPB98 exponents | **Proved** (Theorem 2.1) |
| Palindromic block decomposition with zero cross-coupling | **Proved** (Theorem 3.2) |
| Legendre character separation (Euler's criterion) | **Proved** (Theorem 4.1) |
| Prefactor C = 8/√5 (0.005% match to experiment) | **Identified** (Theorem 5.1; physical mechanism open) |
| √5 absent from splitting field of both block characteristic polynomials | **Proved** (Theorem 6.1) |
| p\_max-exclusion: largest prime absent from both block determinants | **Verified** (52 squarefree moduli; Conjecture 8.2) |
| Polynomial hierarchy: R(S ∪ {q}) exactly quadratic in q | **Verified** (4 levels, 48 tests; Theorem 8.6) |
| Hierarchy **breaks** at m = 30030: R(30030) = −48317287/3 | **Proved** (Theorem 9.1) |
| 3-adic fracture mechanism: v₃(det\_odd) = 947 > v₃(det\_even) = 946 | **Computed** (Theorem 9.2) |

**Companion papers:**
- "Spectral Isotropy and the Exact Temperature of the Prime Gas," Matos (2026). DOI: [10.5281/zenodo.19156532](https://doi.org/10.5281/zenodo.19156532)
- "The Prime Column Transition Matrix Is a Boltzmann Distribution at Temperature ln(N)," Matos (2026). DOI: [10.5281/zenodo.19076680](https://doi.org/10.5281/zenodo.19076680)
- "Separable Gyro-Bohm Scaling for Fusion Energy Confinement: Resolving the Isotope Sign Conflict," Matos (2026). DOI: [10.5281/zenodo.19117880](https://doi.org/10.5281/zenodo.19117880)

## Repository structure

```
paper/
  PRIMORIAL_THERMODYNAMICS.md       # Markdown source
  PRIMORIAL_THERMODYNAMICS.tex      # LaTeX source
  PRIMORIAL_THERMODYNAMICS.pdf      # Compiled PDF (26 pages)
scripts/
  hierarchy_termination_proof.py    # Master proof: 88 assertions across 5 primorials (~3 min)
  gauss_sum_primorial_hierarchy.py  # Master verification for m=30 and m=210 (82 assertions)
  gauss_sum_mechanism.py            # Original m=30 verification (§§2–7, 39 assertions)
  verify_R30030_final.py            # Definitive 13-prime verification of R(30030) (~8 min)
  verify_R30030_definitive.py       # Composite-modulus bug exposure; vacuous-pass proof
  d_sym_210.py                      # 48×48 matrix construction, block decomposition
  hunt_the_seven.py                 # 60+ dimensionless combination scan across 7 devices
  m210_honest.py                    # Sensitivity analysis, m=30 vs m=210 operating window
  m2310_sanity_check.py             # 480×480 matrix, character separation, eigenvalues
  m2310_exact_det.py                # Exact 240×240 integer determinants via Bareiss (~80s)
  ratio_pattern_hunter.py           # R(m) for 50 squarefree m; p_max-exclusion universality
  ratio_deep_hunt.py                # R(p) formula, R(2m)=R(m), polynomial hierarchy (40 tests)
  ratio_level3.py                   # Exact 288×288 determinants for m=1365; level 3 polynomial
  ratio_verify_and_predict.py       # Level 2 verification, level 3 prediction at q=17
  R30030_crack.py                   # Rational reconstruction of R(30030) via half-GCD
  bk_elimination.py                 # 20-round elimination of b_k recurrence hypotheses
  v3_overnight.py                   # Exact 3-adic valuations via mod-3^k elimination (~6 min)
  v3_rabbit_chase.py                # GF(3) nullity, v_3 table, full v_p structure
LICENSE
```

## Requirements

- Python 3.11+
- NumPy
- SymPy

All scripts are self-contained with no external dependencies beyond NumPy and SymPy. No installation required.

## Running the scripts

```bash
# Master proof script — verifies all claims in the paper
pip install numpy sympy
python scripts/hierarchy_termination_proof.py
# Expected output: RESULTS: 88 passed, 0 failed, 88 total
# Runtime: ~3 minutes (dominated by m=2310 Bareiss and m=30030 modular elimination)

# Original m=30 and m=210 verification
python scripts/gauss_sum_primorial_hierarchy.py

# Definitive R(30030) verification (13 independent prime checks)
python scripts/verify_R30030_final.py    # ~8 minutes

# 3-adic fracture mechanism
python scripts/v3_overnight.py           # ~6 minutes
```

## Citation

```bibtex
@article{matos2026primorial,
  title   = {Primorial Thermodynamics: {G}auss Sums, the {IPB98} Scaling Law, and the Hierarchy Break at $m = 30030$},
  author  = {Matos, Antonio P.},
  year    = {2026},
  doi     = {10.5281/zenodo.19188924},
  url     = {https://doi.org/10.5281/zenodo.19188924},
  note    = {Preprint}
}
```

## License

[MIT](LICENSE)

---

_The rabbit has been caught._
