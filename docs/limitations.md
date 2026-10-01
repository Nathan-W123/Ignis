# Known limitations

This document exists because a simulator that does not state its limits is a
simulator that will be misused. Everything below is a real, deliberate
restriction of the model as implemented — not a wish list.

Read it alongside [`validation.md` §6](validation.md), which says what has and
has not been validated.

---

## 0. The headline

**Ignis-M1, Ignis-H1 and Ignis-K1 are conceptual engines invented for this project.**
They are not models of any real hardware. Their predicted thrust, specific
impulse, heat flux, wall temperatures and pressure drops have never been
compared with a test stand. Nothing in this repository is flight-ready, and
nothing here should be used to justify a hardware decision.

The tool is validated to reproduce **NASA CEA's idealisation** of rocket
performance — worst case 0.18 % in flame temperature and 0.06 % in c\* over
156 cases ([`validation.md`](validation.md)). Real engines depart from that
idealisation by amounts that Ignis parameterises (η<sub>c\*</sub>,
λ<sub>div</sub>, η<sub>nozzle</sub>, correlation multipliers) but does not
predict. Reproducing an idealisation accurately is not the same thing as
predicting an engine, and this document is the difference.

One of those parameters is no longer an assumption. The gas-side heat transfer
has been compared against two measurements
([`validation.md` §5b](validation.md)). Bartz's correlation over-predicts a
heated-air nozzle by about 45 % at high pressure and 150 % at low pressure,
and a LOX/GH<sub>2</sub> rocket's throat by 1.72 while matching its chamber to
1 %. Its error is positional, because Bartz carries no boundary-layer history.
Ignis's default is now an integral turbulent boundary layer marched from where
the layer starts ([`theory.md` §12.1](theory.md)). It lands on the
high-pressure air data (median 1.00 against Bartz's 1.45) and halves Bartz's
shape error in the rocket (rms 0.134 against 0.320).

Three things follow, and they are the honest summary of what is and is not
known about the thermal side of this tool:

* The boundary layer is not accurate everywhere. On the low-pressure air tests
  it is as wrong as Bartz (×2.3–2.6), because the layer there is laminarising,
  which a turbulent closure cannot represent; the run flags it. In the rocket
  it reproduces only half of the measured throat dip. The source authors warn
  that injection and combustion can put any mean-flow prediction off by up to
  a factor of two in a real engine.
* The M1's wall-temperature margin (692 K against CuCrZr's 800 K, with the 3 %
  fuel film the boundary layer showed it needs) is smaller than that
  uncertainty, and so is not a demonstration of a feasible chamber.
* The Monte Carlo dispersion on the hot-gas multiplier (lognormal,
  σ<sub>ln</sub> = 0.15) was too narrow for Bartz: the two measured factors sat
  2.5 σ and 3.6 σ out. Against the boundary layer the high-pressure air and the
  rocket's shape error sit inside it, and the laminarising tests and the
  factor-of-two caution do not. It is kept as a judgement, not a measurement.

---

## 1. Chemistry

| Limitation | Consequence |
|---|---|
| **A nozzle is a wall radius and nothing else.** Geometry enters as `r(x)`, whether from the analytic parameterisation or an imported contour. | This is what a quasi-1D internal-flow model is, not a gap in the importer. A CAD model's manifolds, injector face, channel routing and mounting hardware have nowhere to go and are discarded. Importing a real engine gives you its flow passage, not its hardware. |
| **Gas phase only.** No condensed species are carried. | Fuel-rich hydrocarbon cases that would form solid carbon are wrong. For LOX/CH<sub>4</sub> this matters below roughly O/F 1.5; the shipped sweeps stay at O/F ≥ 2. The exception the code makes on purpose is a fuel-rich hydrocarbon gas generator or preburner (O/F ≈ 0.3), whose gas is taken as equilibrium at the turbine temperature: its temperature is the one asked for, but its composition (no soot, no cracked fuel) is approximate, and the cycle says so. Metallised or chlorine-bearing propellants are out of scope entirely. |
| **Equilibrium only — no finite-rate kinetics.** | Chamber composition assumes complete mixing and infinite residence time. The real recombination lag in a nozzle lies *between* the frozen and shifting limits; both bounds are computed and reported, but the true answer is not. |
| **Ideal-gas *thermal* equation of state for the products** (p = ρRT, Z ≡ 1). The *caloric* behaviour is fully variable — c<sub>p</sub>(T), γ(T) and the molar mass all change with the state — which is why the nozzle and shock solvers are described as "variable-property equilibrium-gas" rather than "real-gas". No compressibility factor, fugacity or non-ideal mixing rule is implemented. | At the highest condition used here (20 MPa, 3600 K) the product mixture sits at ≈ 14.5 kg/m³, a molar volume of 1.5 L/mol and a reduced temperature above 15 — deep in the ideal-gas regime, so the approximation is sound at these pressures. It would need revisiting for a very-high-pressure staged-combustion chamber. |
| **40 species, restricted by element set.** | Species not in `data/thermo/ignis_nasa7.yaml` simply do not exist for the solver. Adding one is a data edit, not a code change, but it is an edit. |
| **NASA TM-4513 (1993) 7-coefficient data.** | A few parts in 10³ of flame-temperature difference against the CEA 2002 9-coefficient set; this is the dominant term in the CEA comparison and is measured, not assumed. |
| **No ionisation.** | Irrelevant below ~5000 K, which covers every case here, but it is a hard ceiling. |

## 2. Nozzle flow

| Limitation | Consequence |
|---|---|
| **Quasi-1D.** | Properties vary only along the axis. No radial profiles, no shock structure, no Mach-disc geometry. Divergence loss enters as the analytic λ = (1 + cos θ<sub>e</sub>)/2 factor, which is a cone result applied to a bell. |
| **The boundary layer is an integral correction to an inviscid core.** | Its displacement thickness narrows the throat (a discharge coefficient) and the exit (a smaller effective expansion ratio). Its momentum thickness is a thrust deficit ([`theory.md` §9.3](theory.md)). That costs about 1 % of vacuum I<sub>sp</sub> on the shipped engines. The layer does not reshape the core beyond the area ratio, takes no part in the shock and separation logic, is turbulent from its start, and away from a cooling jacket runs on an assumed wall temperature (`uncooled_wall_temperature`, 1000 K by default), which it reports. |
| **Separation is predicted, not modelled.** | The Summerfield and Schmucker criteria are correlations of test data. Ignis reports *where* separation is predicted and flags it; it does **not** alter the inviscid solution, so a separated nozzle's reported thrust is the attached-flow thrust and is optimistic. The two criteria disagree with each other by a noticeable margin, which is itself informative. Both are fits for *free-shock* separation and can be conservative: Summerfield flags the RS-25 at 104.5 % power as separating at sea level (p<sub>e</sub>/p<sub>a</sub> = 0.20), where the real engine runs with its thrust-optimised bell flowing full ([`data/engines/rs25/README.md`](../data/engines/rs25/README.md)). |
| **Internal normal shocks are 1-D.** | When an over-expanded nozzle can support an internal shock, it is placed by a 1-D Rankine–Hugoniot solve. Real over-expanded nozzles form oblique shock systems and separate first; the report says so explicitly rather than presenting the 1-D answer as the physical one. |
| **Contour is analytic, not method-of-characteristics.** | The bell is the Rao-type parabolic approximation (Huzel & Huang construction) with C¹ joins, not an MOC-designed contour. An MOC nozzle of the same length and area ratio would perform slightly better. |
| **A film is in the thermal solution, not the flow solution.** | A wall film sets the wall's driving temperature (Hatch & Papell) and its cost in specific impulse is bracketed between fully mixed (the headline) and an unmixed two-stream limit. The core flow does not carry the film's composition or temperature, and where between the two limits a real film lands is not predicted. No transpiration cooling. |

## 3. Thermal and cooling

This is the least certain part of the model, and it is the part most likely to
be mistaken for a prediction.

| Limitation | Consequence |
|---|---|
| **A turbulent integral boundary layer for the hot-gas coefficient** (Bartz's correlation as an option). | 1/7-law and Crocco–Busemann profiles, an Eckert reference state, Blasius skin friction and a Colburn analogy, started at Preston's minimum turbulent Re<sub>θ</sub>. Against measurements ([`validation.md` §5b](validation.md)) it lands on high-pressure air data but is as wrong as Bartz (×2.3–2.6) where the layer laminarises. The run flags the acceleration parameter but does not model relaminarisation. In a rocket it gets about half the measured throat dip. Injection and combustion secondary flows, which it does not represent, can put any such prediction off by up to a factor of two. The Monte Carlo campaign disperses it (`cooling.hot_gas_multiplier`, lognormal σ<sub>ln</sub> = 0.15), and it is the dominant driver of wall temperature and pressure drop. |
| **Film cooling is a 1959 correlation used beyond its data.** | Hatch & Papell's slot-film effectiveness was fitted to helium films in heated air; it reproduces its own data to 4 % in its fitted range ([`validation.md` §5c](validation.md)). Ignis applies it to a rocket's fuel film through the coupled jacket solve. The film's velocity ratio is clamped to the fitted range [0.45, 33.3] with a warning. There is no liquid film, no film evaporation or combustion, and no injector-scale mixing. |
| **RP-1 is cooled as n-dodecane.** | No reference equation of state exists for RP-1 in CoolProp, so the jacket uses the single-component surrogate: 7 % light in density, and valid only to 700 K. Coking is a wall-temperature threshold (`coolant_wall_limit`, from NASA's RP-2 deposit onset), not a deposit model, and the shipped K1 example passes it on purpose. |
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
| **One equivalent turbine on all the pumps.** | The `cycle` block closes gas-generator, staged-combustion and expander power balances ([`theory.md` §14](theory.md)). A twin-shaft engine balances each shaft; Ignis checks their sum. Split fuel paths, second preburners and boost turbopumps (the RS-25 has all three) are not represented; the boost pumps appear only as a pump inlet pressure. |
| **Efficiencies and temperatures are inputs.** | Pump, turbine and mechanical efficiencies, turbine inlet temperature and injector stiffnesses are design inputs. There are no pump or turbine maps, no speeds, no NPSH or cavitation check, no bearing-coolant, igniter or pressurisation bleeds and no start transient. The shipped engines' values are labelled assumptions. |
| **Components validated; closures verified.** | Pumps, turbines and the preburner energy balance reproduce the RS-25's published turbopump data to 0.1–4.5 % ([`validation.md` §5d](validation.md)). The closures themselves have been checked for internal consistency, not against a whole engine's balance. |
| **Ideal-gas turbine expansion.** | Generator and preburner gas expands as an ideal-gas mixture. At the RS-25's fuel-preburner state (994 K, 33 MPa), pure hydrogen has Z = 1.064, and its real isentropic work over that turbine's 1.50 pressure ratio is 5.6 % more than the ideal gas's. Those figures are from CoolProp's hydrogen equation of state; the steam in the real gas is not included, and Ignis does not compute either. Ignis's turbine there comes out 4.5 % low, which fits. Expander turbines use the real-fluid tables and do not share this. |
| **Pumps stop at 60 MPa.** | The real-fluid tables end there. A cycle whose pumps would need more is reported as not closing, with the ceiling named. |
| **Incompressible, isothermal lines.** | Propellant density and viscosity are constants read from the library. Cavitation, two-phase flow and line dynamics are absent. |
| **No injector model.** | The injector is a discharge coefficient and an area. Element type, mixing quality, atomisation and their effect on η<sub>c\*</sub> are all outside the model — which is exactly why η<sub>c\*</sub> is an input, not an output. |

## 5. Transient

| Limitation | Consequence |
|---|---|
| **Zero-dimensional.** | One pressure, one temperature, one composition for the whole chamber. No wave dynamics, no acoustic modes, and therefore **no combustion-stability analysis of any kind**. |
| **Products-based: valid from ignition onward.** | The chamber is filled with equilibrium combustion products at all times. It cannot represent the cold pre-ignition fill, the ignition overpressure spike, or a hard start. The start it computes begins the moment heat release begins. |
| **Ignition is a prescribed heat-release fraction.** | `ignition.values` is an input schedule, not a model of an igniter. The shipped plateau (0.874) was *calibrated* so the steady part of the run reproduces the pressure the steady model implies for the same commanded flows, to 0.031 %; that makes the two a cross-check on each other, but it does not predict ignition. |
| **Quasi-steady nozzle.** | The transient's thrust estimate scales the converged steady axial solution by the instantaneous chamber pressure. The nozzle flow is not re-solved at every instant. Both the report and the animation say so. |
| **Integration order degrades across p<sub>c</sub> = p<sub>amb</sub>.** | The sub-critical outflow relation has a √ singularity there, crossed once during start-up. The measured order over the shipped case is ≈ 1.5; on a fully choked problem it is 4.15/4.05/4.17 as RK4 requires. The adaptive integrator handles it by taking small steps, and the conservation residuals stay at 1e-15/1e-13, but the formal order claim only holds away from that crossing. |
| **Tabulated EOS discretisation.** | Chamber properties come from a 41 × 97 × 13 cubic-Hermite table, not a live equilibrium solve. The measured error is ≤ 0.13 % in temperature; it is reported on every run so it cannot be forgotten. |

## 6. Optimisation and uncertainty

| Limitation | Consequence |
|---|---|
| **Nelder–Mead is a local method.** | Latin-hypercube multi-start mitigates this, but nothing here proves global optimality. The reported design is the best feasible point found in the declared evaluation budget, and it is described that way. |
| **A derivative-free search needs the feasible set to be findable.** | The shipped trade study's feasible set is narrow — a 20 kN-wide thrust band intersected with a wall-temperature limit and a jacket budget. Inside the declared engineering bounds, 12 starts find it reliably and converge to the same optimum. Widening the variable box to 15 MPa and a 40 mm throat makes the search **fail outright**: 6000 evaluations, no feasible point found. That is a property of the method, not of the physics, and Ignis reports it as an `InfeasibleError` rather than returning an infeasible design. It also means the declared variable bounds are doing real work and have to be defensible engineering limits, not arbitrary box edges. |
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
   on this list is secondary to that. Partially begun: the gas-side correlation
   now has a measured error against both a heated-air nozzle and a LOX/GH<sub>2</sub>
   heat-sink rocket ([`validation.md` §5b](validation.md)). What those two do
   *not* cover is the rest of the chain — neither has a regeneratively cooled
   wall, so the coupled gas/wall/coolant solve is still unvalidated end to end;
   neither constrains η<sub>c\*</sub>, because both take the combustion gas as
   given rather than predicting it; and neither measures delivered
   I<sub>sp</sub> at all. A firing of an instrumented regeneratively cooled
   engine remains the thing needed.
2. **A boundary-layer solution** coupled to the core flow. Partly done: an
   integral turbulent layer now sets the hot-gas coefficient and the nozzle's
   viscous losses. What it lacks is relaminarisation, which both
   low-pressure air tests and the rocket's throat point to, and a viscous–
   inviscid interaction that reshapes the core.
3. **A conjugate 2-D or 3-D thermal solution** of the wall and channel,
   replacing the 1-D fin model, with a thermal-stress and life criterion.
4. **Finite-rate nozzle kinetics**, replacing the frozen/shifting bracket with
   an actual recombination calculation.
5. **A closed cycle balance.** Partly done: gas-generator, staged-combustion
   and expander balances close for a given chamber pressure and say when they
   cannot. Turbomachinery maps, separate shafts and a validation against a
   whole engine's balance are what is missing, and chamber pressure is still
   an input.
6. **Method-of-characteristics contour design**, replacing the Rao
   approximation.
7. **A multi-zone or CFD-coupled transient**, to say anything at all about
   ignition transients and stability.
