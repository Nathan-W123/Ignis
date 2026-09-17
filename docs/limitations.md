# Known limitations

This document exists because a simulator that does not state its limits is a
simulator that will be misused. Everything below is a real, deliberate
restriction of the model as implemented — not a wish list.

Read it alongside [`validation.md` §6](validation.md), which says what has and
has not been validated.

---

## 0. The headline

**Ignis-M1 and Ignis-H1 are conceptual engines invented for this project.**
They are not models of any real hardware. Their predicted thrust, specific
impulse, heat flux, wall temperatures and pressure drops have never been
compared with a test stand. Nothing in this repository is flight-ready, and
nothing here should be used to justify a hardware decision.

The tool is validated to reproduce **NASA CEA's idealisation** of rocket
performance to a few parts in 10⁴ ([`validation.md`](validation.md)). Real
engines depart from that idealisation by amounts that Ignis parameterises
(η<sub>c\*</sub>, λ<sub>div</sub>, η<sub>nozzle</sub>, correlation multipliers)
but does not predict.

---

## 1. Chemistry

| Limitation | Consequence |
|---|---|
| **Gas phase only.** No condensed species are carried. | Fuel-rich hydrocarbon cases that would form solid carbon are wrong. For LOX/CH<sub>4</sub> this matters below roughly O/F 1.5; the shipped sweeps stay at O/F ≥ 2. Metallised or chlorine-bearing propellants are out of scope entirely. |
| **Equilibrium only — no finite-rate kinetics.** | Chamber composition assumes complete mixing and infinite residence time. The real recombination lag in a nozzle lies *between* the frozen and shifting limits; both bounds are computed and reported, but the true answer is not. |
| **Ideal-gas equation of state for the products.** | At 20 MPa and 3600 K the compressibility factor of the product mixture departs from 1 by well under a percent, so this is a good approximation for the pressures here. It would not be for a very-high-pressure staged-combustion chamber. |
| **40 species, restricted by element set.** | Species not in `data/thermo/ignis_nasa7.yaml` simply do not exist for the solver. Adding one is a data edit, not a code change, but it is an edit. |
| **NASA TM-4513 (1993) 7-coefficient data.** | A few parts in 10³ of flame-temperature difference against the CEA 2002 9-coefficient set; this is the dominant term in the CEA comparison and is measured, not assumed. |
| **No ionisation.** | Irrelevant below ~5000 K, which covers every case here, but it is a hard ceiling. |

## 2. Nozzle flow

| Limitation | Consequence |
|---|---|
| **Quasi-1D.** | Properties vary only along the axis. No radial profiles, no boundary layer, no shock structure, no Mach-disc geometry. Divergence loss enters as the analytic λ = (1 + cos θ<sub>e</sub>)/2 factor, which is a cone result applied to a bell. |
| **Inviscid core.** | There is no boundary-layer displacement thickness, so the effective throat area is the geometric one. Real engines lose ~1 % of c<sub>F</sub> here. |
| **Separation is predicted, not modelled.** | The Summerfield and Schmucker criteria are correlations of test data. Ignis reports *where* separation is predicted and flags it; it does **not** alter the inviscid solution, so a separated nozzle's reported thrust is the attached-flow thrust and is optimistic. The two criteria disagree with each other by a noticeable margin, which is itself informative. |
| **Internal normal shocks are 1-D.** | When an over-expanded nozzle can support an internal shock, it is placed by a 1-D Rankine–Hugoniot solve. Real over-expanded nozzles form oblique shock systems and separate first; the report says so explicitly rather than presenting the 1-D answer as the physical one. |
| **Contour is analytic, not method-of-characteristics.** | The bell is the Rao-type parabolic approximation (Huzel & Huang construction) with C¹ joins, not an MOC-designed contour. An MOC nozzle of the same length and area ratio would perform slightly better. |
| **No film or transpiration cooling in the flow solution.** | A film-cooled wall changes the near-wall gas composition and temperature; Ignis models neither. |

## 3. Thermal and cooling

This is the least certain part of the model, and it is the part most likely to
be mistaken for a prediction.

| Limitation | Consequence |
|---|---|
| **Bartz correlation for the hot-gas coefficient.** | Bartz has a documented scatter of roughly ±20–30 % against measured rocket heat flux, and it is worst exactly where it matters most — near the throat where curvature effects are strong. The Monte Carlo campaign disperses it explicitly (`cooling.bartz_multiplier`, lognormal σ<sub>ln</sub> = 0.15) because pretending it is exact would be dishonest. It is, unsurprisingly, the dominant driver of wall temperature (SRC 0.994) and pressure drop (SRC 0.915). |
| **1-D wall, no axial conduction.** | Each station conducts only radially. Near the throat, where the flux gradient is steepest, axial conduction genuinely redistributes heat and the real peak is a little lower and broader than reported. |
| **Rectangular-fin channel model.** | Fin efficiency uses the straight-fin relation. Real channels have fillets, varying aspect ratio and curvature-induced secondary flow. |
| **Dittus–Boelter / Gnielinski coolant Nusselt numbers.** | Both are smooth-tube correlations for fully developed turbulent flow. Near-critical methane and supercritical hydrogen show property variation across the boundary layer that neither captures; the two disagree by up to 20 % on the shipped case, and `cooling.nusselt_multiplier` exists so that disagreement can be propagated. |
| **No thermal-stress or life analysis.** | The wall temperature limit is a plain material number, not a low-cycle-fatigue criterion. Real chamber life is set by strain ratcheting ("dog-house" failure), which is not modelled. |
| **Steady thermal state only.** | The jacket is solved at equilibrium. There is no soak-back, no transient wall heating, no start-up thermal shock. |
| **Wall conductivity is extrapolated at the cold end.** | The CuCrZr fit is valid over 250–900 K; the coolant enters at 111.66 K, so the first few jacket segments extrapolate. Ignis **warns about this by name and gives the range it reached** rather than hiding it. The effect is conservative there and the peak-flux station is inside the fitted range. |
| **Radiation is a gray-gas estimate.** | Disabled by default (`gas_emissivity: 0.0`). When enabled it uses a single emissivity, not a band model; for LOX/CH<sub>4</sub> the real H<sub>2</sub>O/CO<sub>2</sub> radiation is a few percent of total flux. |

## 4. Feed system

| Limitation | Consequence |
|---|---|
| **No cycle balance.** | Ignis sizes the pressures a feed system must supply and reports an *equivalent* pump power. It does not close a gas-generator, staged-combustion or expander cycle: there is no turbine, no preburner, no power balance. |
| **Incompressible, isothermal lines.** | Propellant density and viscosity are constants read from the library. Cavitation, two-phase flow and line dynamics are absent. |
| **No injector model.** | The injector is a discharge coefficient and an area. Element type, mixing quality, atomisation and their effect on η<sub>c\*</sub> are all outside the model — which is exactly why η<sub>c\*</sub> is an input, not an output. |

## 5. Transient

| Limitation | Consequence |
|---|---|
| **Zero-dimensional.** | One pressure, one temperature, one composition for the whole chamber. No wave dynamics, no acoustic modes, and therefore **no combustion-stability analysis of any kind**. |
| **Products-based: valid from ignition onward.** | The chamber is filled with equilibrium combustion products at all times. It cannot represent the cold pre-ignition fill, the ignition overpressure spike, or a hard start. The start it computes begins the moment heat release begins. |
| **Ignition is a prescribed heat-release fraction.** | `ignition.values` is an input schedule, not a model of an igniter. The shipped plateau (0.874) was *calibrated* so the steady part of the run reproduces the steady model's 5.5 MPa to 0.02 %; that makes the two a cross-check on each other, but it does not predict ignition. |
| **Quasi-steady nozzle.** | The transient's thrust estimate scales the converged steady axial solution by the instantaneous chamber pressure. The nozzle flow is not re-solved at every instant. Both the report and the animation say so. |
| **Integration order degrades across p<sub>c</sub> = p<sub>amb</sub>.** | The sub-critical outflow relation has a √ singularity there, crossed once during start-up. The measured order over the shipped case is ≈ 1.5; on a fully choked problem it is 4.15/4.05/4.17 as RK4 requires. The adaptive integrator handles it by taking small steps, and the conservation residuals stay at 1e-15/1e-13, but the formal order claim only holds away from that crossing. |
| **Tabulated EOS discretisation.** | Chamber properties come from a 41 × 97 × 13 cubic-Hermite table, not a live equilibrium solve. The measured error is ≤ 0.13 % in temperature; it is reported on every run so it cannot be forgotten. |

## 6. Optimisation and uncertainty

| Limitation | Consequence |
|---|---|
| **Nelder–Mead is a local method.** | Latin-hypercube multi-start mitigates this, but nothing here proves global optimality. The reported design is the best feasible point found in 2500 evaluations, and it is described that way. |
| **Derivative-free, so no optimality certificate.** | There are no KKT multipliers to check. What *is* checked is that the returned design re-analyses to the reported objective and satisfies every constraint to 1e-9. |
| **Monte Carlo propagates the uncertainties you declare.** | The campaign disperses ten inputs. It cannot discover an uncertainty that was not declared — model-form error above all. A 2000-sample campaign with tight distributions will report tight outputs whether or not the model is right. |
| **Standardised regression coefficients assume near-linearity.** | Each ranking carries its linear-model R², and they are 0.96–1.00 for the shipped case, so the rankings are meaningful there. For a strongly non-linear response the SRC ranking would be misleading and the R² would say so. |
| **Input correlations are ignored.** | All dispersed inputs are drawn independently. Real manufacturing tolerances are correlated. |

## 7. Software

| Limitation | Consequence |
|---|---|
| **Single process, shared memory.** | No MPI, no GPU. The largest thing here (a 2000-sample campaign) takes a minute on four threads, so this has not been a constraint. |
| **Double precision throughout, no interval arithmetic.** | Residuals near 1e-14 are at the level where round-off dominates; tightening tolerances further is not meaningful. |
| **Linux-focused.** | CI covers GCC and Clang on Ubuntu. The code is standard C++17 with no POSIX-only calls, and `constants::pi` exists rather than `M_PI` precisely so MSVC would work, but Windows and macOS are untested. |
| **No plugin or scripting interface.** | Extending the model means writing C++ and rebuilding. [`extending.md`](extending.md) describes where. |

---

## What would have to change for this to be a design tool

Honestly stated, in rough order of importance:

1. **Validation against measured engine data** — heat flux, wall temperature and
   delivered I<sub>sp</sub> from an instrumented test article. Everything else
   on this list is secondary to that.
2. **A boundary-layer solution** coupled to the core flow, replacing both the
   inviscid-throat assumption and the Bartz correlation.
3. **A conjugate 2-D or 3-D thermal solution** of the wall and channel,
   replacing the 1-D fin model, with a thermal-stress and life criterion.
4. **Finite-rate nozzle kinetics**, replacing the frozen/shifting bracket with
   an actual recombination calculation.
5. **A closed cycle balance** — turbine, preburner, power balance — so that
   chamber pressure becomes an output of the cycle rather than an input.
6. **Method-of-characteristics contour design**, replacing the Rao
   approximation.
7. **A multi-zone or CFD-coupled transient**, to say anything at all about
   ignition transients and stability.
