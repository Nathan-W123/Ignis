# Configuration reference

Every Ignis run is driven by one YAML file. The applications add no physics of
their own: what is in the file is what is solved. Unknown top-level keys are
rejected, and every error message names the offending YAML path and file, for
example

```
configs/methane_nominal.yaml: nozzle.expansion_ratio: expected a value > 1, got 0.5
```

All quantities are **SI** — Pa, K, m, kg/s, W — except where a key name says
otherwise (`..._angle` keys are degrees, `..._ratio` and `..._fraction` keys
are dimensionless).

## Top-level sections

| Section | Read by | Required |
|---|---|---|
| `species` | all | no (defaults to C/H/O) |
| `propellants` | all | no (defaults to LOX/LCH4) |
| `chamber` | all | no |
| `nozzle` | all but `ignis_equilibrium` | yes in practice |
| `performance` | `ignis_engine`, `ignis_nozzle`, studies | no |
| `cooling` | `ignis_engine`, studies | no (off by default) |
| `feed` | `ignis_engine`, studies | no (off by default) |
| `transient` | `ignis_transient` | yes for that tool |
| `sweep` | `ignis_sweep` | yes for that tool |
| `optimization` | `ignis_sweep --optimize` | yes for that mode |
| `monte_carlo` | `ignis_mc` | yes for that tool |
| `output` | all | no (defaults to `results/`) |

---

## `species`

```yaml
species:
  database: data/thermo/ignis_nasa7.yaml   # optional; the shipped DB is found automatically
  elements: [C, H, O]                      # keep every species built only from these
  include: [CO2, H2O, CO, H2, O2, OH, H, O]   # explicit list; overrides `elements`
```

`elements` is the normal choice: it keeps the solve to the chemistry that can
actually occur. Adding N to a LOX/CH4 case only adds dead columns. `include`
is for deliberately frozen or reduced-mechanism studies.

## `propellants`

```yaml
propellants:
  library: data/propellants/ignis_propellants.yaml   # optional
  oxidizer: LOX            # LOX, GOX
  fuel: LCH4               # LCH4, LH2, GCH4, GH2
  oxidizer_temperature: 90.18   # K; omit or 0 => the library reference state
  fuel_temperature: 111.66      # K
  mixture_ratio: 3.4            # O/F by mass
```

The propellant library carries the **absolute** enthalpy of each propellant in
its stated state (liquid at its normal boiling point, for the cryogens), so the
chamber solve is a genuine constant-enthalpy problem from the tanks, not from
gaseous reactants at 298.15 K. Setting a temperature away from the reference
value applies the tabulated liquid c<sub>p</sub>; going far outside the
tabulated range raises `RangeError` rather than extrapolating silently.

## `chamber`

```yaml
chamber:
  pressure: 5.5e6          # Pa, stagnation
  eta_c_star: 0.96         # empirical combustion efficiency, multiplies c*
  composition: equilibrium # equilibrium (shifting) | frozen
```

`eta_c_star` is an **empirical factor**, not a physical model: it scales the
ideal characteristic velocity and hence the mass flow at a given pressure. All
outputs report both the ideal and the corrected value so the two are never
confused.

## `nozzle`

```yaml
nozzle:
  throat_radius: 0.070          # m   (or throat_area, m^2 — give exactly one)
  contraction_ratio: 2.8        # Ac/At (or chamber_radius, m)
  chamber_length: 0.22          # m, cylindrical section
  converging_half_angle: 30.0   # deg
  chamber_fillet_ratio: 0.5     # R1/Rc
  throat_upstream_ratio: 1.5    # Ru/Rt
  throat_downstream_ratio: 0.382 # Rd/Rt   (Rao's value)
  expansion_ratio: 20.0         # Ae/At
  type: bell                    # bell | conical
  bell_length_fraction: 0.8     # of the equivalent 15 deg cone (bell only)
  bell_initial_angle: 33.0      # deg (bell only)
  bell_exit_angle: 8.0          # deg (bell only)
  cone_half_angle: 15.0         # deg (conical only)
  stations: 400                 # contour stations
  bartz_curvature: mean         # mean | throat | none — which radius feeds Bartz
```

## `performance`

```yaml
performance:
  altitude: 0.0                 # m; sets the ambient pressure from USSA-1976
  ambient_pressure: 101325.0    # Pa; used when `altitude` is absent
  auto_divergence: true         # lambda = (1 + cos theta_e)/2
  eta_nozzle: 1.0               # extra multiplicative nozzle efficiency
  separation: summerfield       # summerfield | schmucker | none
  resolve_internal_shocks: true # search for an internal normal shock when over-expanded
```

`altitude` and `ambient_pressure` are mutually exclusive in effect — whichever
appears last in the parameter registry wins when a study drives them.

## `cooling`

```yaml
cooling:
  enabled: true
  coolant: methane              # methane | hydrogen | oxygen (tabulated reference EOS)
  material: CuCrZr              # CuCrZr | Copper | Inconel718 | SS316
  num_channels: 300
  width_mode: fraction_of_pitch # fraction_of_pitch | fixed
  width_fraction: 0.50          # channel width / local pitch
  channel_width: 2.0e-3         # m, when width_mode: fixed
  channel_height: 5.0e-3        # m
  wall_thickness: 0.6e-3        # m, hot-gas-side wall
  roughness: 5.0e-6             # m, absolute
  coolant_fuel_fraction: 1.0    # fraction of the fuel flow through the jacket
  inlet_temperature: 111.66     # K
  inlet_pressure: 15.0e6        # Pa
  counterflow: true             # coolant enters at the nozzle end
  x_end_area_ratio: 10.0        # jacket ends where A/At falls to this value
  nusselt_correlation: dittus-boelter  # dittus-boelter | gnielinski
  nusselt_multiplier: 1.0       # correlation uncertainty knob
  bartz_multiplier: 1.0         # correlation uncertainty knob
  gas_emissivity: 0.0           # 0 disables the radiation term
  wall_conductivity: 0.0        # W/(m K); 0 => the material's k(T) fit
  segments: 240                 # jacket march segments
```

`x_end_area_ratio` matters whenever ε is swept or optimised: without it the
jacket would either stop short of, or run past, the end of the contour when the
nozzle changes length. `nusselt_multiplier` and `bartz_multiplier` exist so
that the Monte Carlo campaign can disperse the **correlation** uncertainty
explicitly rather than pretending the correlations are exact.

## `feed`

```yaml
feed:
  enabled: true
  mode: sizing                  # sizing | capability
  injector_stiffness: 0.20      # dp_inj / p_c, the design target
  minimum_stiffness: 0.15       # reported as a violation below this
  pump_efficiency: 0.65
  pump_inlet_pressure: 3.0e5    # Pa
  oxidizer: { line_length: 1.5, line_diameter: 0.10, fitting_k: 2.5,
              viscosity: 1.9e-4, injector_cd: 0.78 }
  fuel:     { line_length: 1.8, line_diameter: 0.07, fitting_k: 3.0,
              viscosity: 1.2e-4, injector_cd: 0.75 }
```

## `transient`

```yaml
transient:
  enabled: true
  chamber_volume: 0.0           # m^3; 0 => from the nozzle geometry
  ambient_pressure: 101325.0
  initial_pressure: 101325.0
  initial_temperature: 300.0
  initial_mixture_ratio: 3.4
  oxidizer_schedule: { open: 0.004, ramp_up: 0.035, value: 37.12,
                       close: 0.600, ramp_down: 0.025 }   # kg/s, s
  fuel_schedule:     { open: 0.002, ramp_up: 0.035, value: 10.92,
                       close: 0.602, ramp_down: 0.026 }
  ignition:
    times:  [0.000, 0.004, 0.020, 1.0e9]
    values: [0.00,  0.05,  0.874, 0.874]   # fraction of q_comb released
  t_end: 0.80
  dt: 2.0e-6                    # fixed step (rk4) / initial step (rk45)
  dt_output: 1.0e-3
  integrator: rk45              # rk4 | rk45 (Cash-Karp)
  rtol: 1.0e-9
  atol: 1.0e-12
  dt_max: 2.0e-4
  table:                        # the pre-tabulated equilibrium EOS
    mr_min: 0.25, mr_max: 12.0, mr_points: 41
    t_min: 200.0, t_max: 4000.0, t_points: 97
    p_min: 5.0e4, p_max: 2.0e7, p_points: 13
  interpolation_samples: 200    # random interior points used to measure table error
  interpolation_seed: 20260917
```

`ignition.values` is the **released fraction of the heat of combustion**, not
the steady c\* efficiency; the two are related but not equal (see
[`theory.md` §16](theory.md)). Ignis reports the table's measured
interpolation error on every transient run so the discretisation is never
invisible.

## `sweep`

```yaml
sweep:
  axes:
    - { parameter: propellants.mixture_ratio, min: 2.0, max: 5.0, points: 31 }
    - { parameter: nozzle.expansion_ratio,    values: [10, 20, 40, 80] }
  metrics: [chamber.temperature, performance.isp, cooling.max_heat_flux]
```

Axes take either `min`/`max`/`points` (linear) or an explicit `values` list.
Multiple axes form the full tensor product. A point that raises
`InfeasibleError` or `ConvergenceError` is recorded as failed with its message
and excluded from the CSV, never written as a plausible-looking number.

## `optimization`

```yaml
optimization:
  objective: performance.isp_vacuum
  sense: maximize               # maximize | minimize
  starts: 6                     # Latin-hypercube multi-start
  max_evaluations: 2500
  outer_iterations: 8           # augmented-Lagrangian outer loops
  simplex_tolerance: 1.0e-6
  initial_penalty: 10.0
  penalty_growth: 5.0
  seed: 20260917
  variables:
    - { parameter: chamber.pressure, min: 3.0e6, max: 11.0e6, start: 5.5e6 }
  constraints:
    - { metric: cooling.max_wall_temperature, op: "<=", bound: 800.0 }
```

`op` is `<=`, `>=` or `==`. Infeasible evaluations are given a finite, large
penalty rather than being discarded, so the simplex can walk back into the
feasible region.

## `monte_carlo`

```yaml
monte_carlo:
  samples: 2000
  seed: 20260917
  threads: 0                    # 0 => hardware concurrency
  record_samples: true          # write the full sample matrix
  fd_step: 1.0e-3               # relative step for the local elasticities
  inputs:
    - { parameter: chamber.pressure, distribution: normal, mean: 5.5e6, sigma: 1.1e5 }
    - { parameter: chamber.eta_c_star, distribution: triangular,
        low: 0.930, mode: 0.960, high: 0.980 }
    - { parameter: cooling.bartz_multiplier, distribution: lognormal,
        median: 1.0, sigma_log: 0.15, min: 0.5, max: 2.0 }
    - { parameter: propellants.fuel_temperature, distribution: uniform,
        low: 108.0, high: 115.0 }
  outputs: [performance.thrust, performance.isp, cooling.max_wall_temperature]
```

Distributions: `normal` (`mean`, `sigma`), `lognormal` (`median`,
`sigma_log`), `uniform` (`low`, `high`), `triangular` (`low`, `mode`, `high`).
Any distribution accepts optional `min`/`max` truncation bounds, applied by
redrawing. Sampling is **thread-count independent**: sample *i* is seeded with
`SplitMix64(seed, i)`, so a 1-thread and an 8-thread run give byte-identical
output. `scripts/run_all.sh` checks this on every full run.

## `output`

```yaml
output:
  directory: results/methane_nominal
  prefix: m1
```

The command-line `-o/--output` overrides `directory`. Files are named
`<prefix>_<kind>.<ext>`.

---

## Parameter names

`sweep`, `optimization` and `monte_carlo` address the configuration through a
registry of named parameters, so a study never has to duplicate the schema.

**Inputs** (`parameter:`)

```
chamber.pressure                chamber.eta_c_star
propellants.mixture_ratio       propellants.oxidizer_temperature
propellants.fuel_temperature
nozzle.throat_radius            nozzle.throat_area
nozzle.expansion_ratio          nozzle.contraction_ratio
nozzle.chamber_length           nozzle.bell_length_fraction
nozzle.bell_exit_angle          nozzle.cone_half_angle
performance.ambient_pressure    performance.altitude
cooling.channel_height          cooling.width_fraction
cooling.num_channels            cooling.wall_thickness
cooling.wall_conductivity       cooling.nusselt_multiplier
cooling.bartz_multiplier        cooling.inlet_temperature
cooling.inlet_pressure          cooling.coolant_fuel_fraction
cooling.roughness               feed.injector_stiffness
```

**Outputs** (`metric:` / `metrics:` / `outputs:`)

```
chamber.temperature        chamber.pressure           chamber.molar_mass
chamber.density            chamber.gamma_s            chamber.gamma_frozen
chamber.cp                 chamber.sound_speed        chamber.c_star_ideal
chamber.c_star             chamber.equivalence_ratio  chamber.throat_temperature
chamber.element_residual   chamber.gibbs_residual     chamber.enthalpy_residual
performance.thrust         performance.thrust_ideal   performance.thrust_momentum
performance.thrust_pressure performance.isp           performance.isp_ideal
performance.isp_vacuum     performance.cf             performance.cf_ideal
performance.c_effective    performance.mdot           performance.exit_pressure
performance.exit_temperature performance.exit_mach    performance.exit_velocity
performance.pressure_ratio performance.expansion_ratio
performance.mass_flow_residual performance.energy_residual
geometry.l_star            geometry.residence_time    geometry.throat_area
geometry.exit_area
cooling.max_wall_temperature cooling.max_heat_flux    cooling.total_heat_load
cooling.pressure_drop      cooling.outlet_temperature cooling.temperature_rise
cooling.energy_balance_residual cooling.flux_residual
feed.oxidizer_tank_pressure feed.fuel_tank_pressure   feed.total_pump_power
```

The residual metrics are deliberately exposed as sweepable quantities: a sweep
can be asked to show that the element balance stays at 1e-13 across its whole
range, which is a far stronger statement than checking it once.

Running `ignis_sweep --list-parameters` prints this registry with
units, generated from the same table the solver uses.
