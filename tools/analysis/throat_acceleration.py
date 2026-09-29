#!/usr/bin/env python3
"""Throat acceleration parameter K for the two validation experiments and the M1.

Backs the numbers quoted in docs/validation.md section 5b.  Run it with no
arguments; it prints the table that section cites.

WHY K.  K = (nu_e / u_e^2)(du_e/dx) measures how strong a favourable pressure
gradient is relative to the viscous scale.  Above roughly 3e-6 turbulence
production in a boundary layer is suppressed faster than it regenerates and the
layer reverts toward laminar, so heat transfer falls below the
equilibrium-turbulent value a correlation like Bartz assumes (Moretti & Kays,
Int. J. Heat Mass Transfer 8, 1187, 1965; Back, Cuffel & Massier, AIAA J. 7,
730, 1969).  It is the first thing to check when asking whether an
over-prediction at the throat is a relaminarization effect.

WHY NOT A NUMERICAL GRADIENT.  Taking du/dx from the 1-D area-Mach relation by
finite differences is useless at the throat: dM/d(A/A*) is singular at M = 1,
so the discrete gradient grows without bound as the grid is refined.  Expanding
both sides about M = 1 gives a finite closed form.  For a throat of radius r_t
whose wall is a circular arc of radius R_c,

    A/A* = 1 + x^2 / (r_t R_c)          (arc tangent at the throat)
    A/A* = 1 + (g+1)/4 (M-1)^2          (area-Mach relation about M = 1)
 => dM/dx = 2 / sqrt((g+1) r_t R_c)

and with u = M a and da/dx = -a (g-1)/(g+1) dM/dx at M = 1,

    (du/dx)_t = 4 a* / [ (g+1)^1.5 sqrt(r_t R_c) ]

VISCOSITY.  Air is Sutherland.  The combustion-product cases use Ignis's own
Chapman-Enskog mixture viscosity, tabulated by running

    ignis_equilibrium -c configs/methane_nominal.yaml --problem tp --temperature T

over the nozzle temperature range and interpolated in log-log.  Chapman-Enskog
is a dilute-gas theory with no density correction, which is the right regime
here; see include/ignis/thermo/Transport.hpp.
"""
import numpy as np

K_CRIT = 3.0e-6            # relaminarization threshold
IN = 0.0254
PSI = 6894.757293168361
RANKINE = 5.0 / 9.0
R_UNIV = 8314.462618

# Ignis Chapman-Enskog mixture viscosity for the M1 products, Pa s.
MU_T = np.array([2100, 2400, 2700, 3000, 3200, 3356, 3400, 3515], float)
MU_V = np.array([7.01923e-05, 7.72883e-05, 8.40764e-05, 9.06913e-05,
                 9.50721e-05, 9.84917e-05, 9.94576e-05, 1.01987e-04], float)


def sutherland(T):
    """Air viscosity, Pa s."""
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def mu_products(T):
    """Combustion-product viscosity, Pa s, from the Ignis table above."""
    return float(np.exp(np.interp(np.log(T), np.log(MU_T), np.log(MU_V))))


def throat_k(label, Tt, pt, gamma, R_gas, r_t, R_c, mu_of, note=""):
    """Print and return K at the throat for one condition."""
    T = Tt / (1 + 0.5 * (gamma - 1))
    p = pt * (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    rho = p / (R_gas * T)
    a = np.sqrt(gamma * R_gas * T)
    dudx = 4 * a / ((gamma + 1) ** 1.5 * np.sqrt(r_t * R_c))
    K = mu_of(T) / rho / a ** 2 * dudx
    print(f"{label:<36} {T:6.0f} {rho:7.3f} {a:7.0f} {dudx:9.0f} "
          f"{K:9.2e} {K / K_CRIT:7.2f}x  {note}")
    return K


def main():
    head = (f"{'case':<36} {'T* K':>6} {'rho*':>7} {'a* m/s':>7} "
            f"{'du/dx 1/s':>9} {'K':>9} {'K/Kcrit':>8}")
    print(f"relaminarization threshold K_crit = {K_CRIT:.1e}\n")
    print(head)
    print("-" * (len(head) + 22))

    print("JPL TR 32-415   air, Tt = 1500 R, r_t = 0.9015 in, R_c = 1.800 in")
    for test, approach, pt_psia, err in [(319, 0, 35.9, 2.66),
                                         (318, 0, 51.0, 2.49),
                                         (271, 18, 175.0, 1.59),
                                         (273, 18, 201.0, 1.44)]:
        throat_k(f"  test {test}  l={approach:>2} in  pt={pt_psia:>5.1f} psia",
                 1500 * RANKINE, pt_psia * PSI, 1.34, 287.05,
                 0.9015 * IN, 1.800 * IN, sutherland,
                 f"Bartz/measured = {err:.2f}x")

    print("\nNASA TN D-2832  LOX/GH2, T_c ~ 3300 K, r_t = 2.5 in")
    print("  (the report does not give R_c, so both plausible bounds are shown;")
    print("   the measured over-prediction at the throat is 1.72x at every p_c)")
    for pc in (157.0, 473.0, 966.0):
        for ratio in (1.0, 2.0):
            throat_k(f"  p_c={pc:>5.0f} psia   R_c={ratio:.0f} r_t",
                     3300.0, pc * PSI, 1.21, R_UNIV / 13.5,
                     2.5 * IN, ratio * 2.5 * IN, mu_products)

    print("\nIgnis-M1        LOX/CH4, p_c = 5.5 MPa, r_t = 70 mm")
    r_t = 0.070
    for label, R_c in [("R_c = Rd   (0.382 r_t)", 0.382 * r_t),
                       ("R_c = mean (0.941 r_t)", 0.5 * (1.5 + 0.382) * r_t),
                       ("R_c = Ru   (1.5 r_t)", 1.5 * r_t)]:
        throat_k(f"  {label}", 3515.4, 5.5e6, 1.14, R_UNIV / 22.4,
                 r_t, R_c, mu_products)

    print("\nThe M1 sits below the threshold -- and so does the rocket at "
          "966 psia,\nwhich is still over-predicted by 1.72x.  Being below it "
          "is not protection.")


if __name__ == "__main__":
    main()
