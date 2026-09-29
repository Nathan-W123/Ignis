#!/usr/bin/env python3
"""Can the JPL dataset attribute its heat-transfer error to any one cause?

Backs the collinearity claim in docs/validation.md section 5b.  Short answer:
no.  The four tests in JPL TR 32-415 differ in cooled approach length AND in
stagnation pressure at the same time, so approach length, throat Reynolds
number and the acceleration parameter K all track the error -- and each other
-- at |r| > 0.98.  With four points that is not separable, and the validation
documentation says so rather than picking the flattering explanation.
"""
import numpy as np

IN = 0.0254
PSI = 6894.757293168361
RANKINE = 5.0 / 9.0
K_CRIT = 3.0e-6

GAMMA, R_GAS, TT = 1.34, 287.05, 1500 * RANKINE
R_T, R_C = 0.9015 * IN, 1.800 * IN

# test, cooled approach [in], stagnation pressure [psia], Bartz/measured median
TESTS = [(319, 0, 35.9, 2.66), (318, 0, 51.0, 2.49),
         (271, 18, 175.0, 1.59), (273, 18, 201.0, 1.44)]


def sutherland(T):
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def main():
    T = TT / (1 + 0.5 * (GAMMA - 1))
    a = np.sqrt(GAMMA * R_GAS * T)
    mu = sutherland(T)
    dudx = 4 * a / ((GAMMA + 1) ** 1.5 * np.sqrt(R_T * R_C))

    print(f"{'test':>5} {'approach':>9} {'pt psia':>8} {'Re_D*':>10} "
          f"{'K':>9} {'K/Kcrit':>8} {'Bartz/meas':>11}")
    print("-" * 66)
    data = []
    for test, approach, pt_psia, err in TESTS:
        rho = (pt_psia * PSI * (2 / (GAMMA + 1)) ** (GAMMA / (GAMMA - 1))
               / (R_GAS * T))
        Re = rho * a * 2 * R_T / mu
        K = mu / rho / a ** 2 * dudx
        data.append((approach, pt_psia, Re, K, err))
        print(f"{test:>5} {approach:>7} in {pt_psia:>8.1f} {Re:>10.2e} "
              f"{K:>9.2e} {K / K_CRIT:>7.2f}x {err:>10.2f}x")

    approach, pt, Re, K, err = map(np.array, zip(*data))

    print("\nEach candidate cause against the error (n = 4):")
    for name, v in [("cooled approach length", approach),
                    ("stagnation pressure", pt),
                    ("throat Reynolds number", Re),
                    ("acceleration parameter K", K)]:
        x = np.log(v) if v.min() > 0 else v
        print(f"  {name:<26} Pearson r = {np.corrcoef(x, err)[0, 1]:+.4f}")

    print("\nThe candidate causes against each other:")
    for name, v in [("log(pt)", np.log(pt)), ("log(Re)", np.log(Re)),
                    ("log(K)", np.log(K))]:
        r = np.corrcoef(approach, v)[0, 1]
        print(f"  approach length vs {name:<8} r = {r:+.4f}")

    print("\nEvery 0-in. test is a low-pressure test and every 18-in. test is a")
    print("high-pressure one.  The dataset cannot tell these causes apart; the")
    print("rocket data in NASA TN D-2832 is what settles it (see 5b).")


if __name__ == "__main__":
    main()
