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
| `cycle` | `ignis_engine`, studies | no (off by default) |
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
  fuel: LCH4               # LCH4, LH2, RP1, GCH4, GH2
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

`RP1` is the NASA CEA library's RP-1, written per carbon atom as CH<sub>1.95</sub>
with CEA's enthalpy (Mehta et al., AIAA 95-2962), at the density NIST measured
(Outcalt, Laesecke & Brumback 2009). Its cooling-jacket table is
`coolant: dodecane`, the n-dodecane surrogate. Run with `species: { elements:
[C, H, O] }`; `configs/kerosene_nominal.yaml` is the worked example.

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

### Importing a real engine's geometry

Everything above describes a nozzle by parameters. A real engine is usually
described by its wall, so the nozzle block also accepts one directly:

```yaml
nozzle:
  contour_file: my_engine_contour.csv   # x,r pairs in metres
  throat_upstream_ratio: 1.5            # only a fallback, see below
  stations: 400
```

The file is two columns, `x` then `r`, both in **metres**, one pair per line,
comma or whitespace separated, with `#` comments and an optional header line
skipped. Nothing richer, because nothing richer is what tools emit.

When `contour_file` is present every other shape key is ignored, and the
quantities they would have set are **measured from the contour** instead:
throat radius and position, contraction and expansion ratio, chamber length,
L\*, the exit wall angle that sets the divergence loss, and the throat
curvature radius Bartz needs (a least-squares circle through the samples within
15 % of the throat radius). `expansion_ratio` and `chamber_length` therefore
stop being required. `throat_upstream_ratio` survives only as a fallback for
the curvature if that circle fit fails, which happens on a sharp throat where
the samples are collinear; the run reports which of the two it used.

The wall between consecutive samples is a **straight line**. The geometry is
exactly the polyline you supply, so its accuracy is the file's resolution —
the throat area comes out slightly large if the throat is sampled coarsely,
and that is not silently corrected.

A contour that is not a nozzle is rejected rather than repaired. Samples must
increase in x, radii must be positive, there must be at least nine of them, the
throat must lie inside the table with room on both sides, and the wall must
contract to the throat and expand after it. Each failure names the station that
caused it.

To get a contour out of a CAD model:

```bash
python3 tools/contour_from_stl.py engine.stl --axis x --scale 0.001 -o contour.csv
```

That reduces a triangle mesh to the revolve profile it was drawn from. It takes
the inner (gas-side) surface by default — `--outer` if your model is of the
flow volume rather than the solid — and it checks how axisymmetric the mesh
actually is about the axis you named, because a model that still has its
flanges on produces a contour that looks reasonable and means nothing.

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
  boundary_layer_losses: true   # march the wall layer and correct C_d, exit and thrust
  uncooled_wall_temperature: 1000.0  # K, the layer's wall where no jacket is solved
  kinetics:                     # optional: finite-rate recombination in the nozzle
    enabled: true
    apply: true                 # scale the delivered thrust by the kinetic efficiency
    mechanism: ""               # default data/kinetics/gri30_nozzle.yaml (GRI-Mech 3.0 subset)
    start_frozen_mach: 1.10     # where the march leaves shifting equilibrium
    rate_multiplier: 1.0        # scales every rate (0 = frozen from the start)
    relative_tolerance: 1.0e-6
  ascent_profile:               # optional; enables performance.isp_ascent
    altitudes: [0, 5000, 10000, 20000, 40000]   # m, geometric
    weights:   [0.30, 0.25, 0.20, 0.15, 0.10]   # time weights, normalised internally
```

`ascent_profile` declares the trajectory over which `performance.isp_ascent` is
averaged. It is an **input to a trade study, not a prediction**: Ignis has no
vehicle model. Thrust is exactly linear in ambient pressure for a full-flowing
nozzle, so the weighted mean is exact rather than a quadrature. Omitting the
section leaves `performance.isp_ascent` at zero.

`altitude` and `ambient_pressure` are mutually exclusive in effect — whichever
appears last in the parameter registry wins when a study drives them.

`kinetics` integrates the H/O/CO recombination along the divergent nozzle at
GRI-Mech 3.0's rates ([`theory.md` §9.4](theory.md)). It needs
`chamber.composition: equilibrium`, since it starts from shifting equilibrium
just past the throat. Its kinetic efficiency (finite-rate over shifting vacuum
impulse) scales the delivered thrust unless `apply: false`, and the march is
written to `<prefix>_kinetics.csv`. A run takes a few tenths of a second more,
so the sweep, optimisation and Monte Carlo configs leave it off.

## `cooling`

```yaml
cooling:
  enabled: true
  coolant: methane              # methane | hydrogen | oxygen | dodecane (tabulated reference EOS)
  material: CuCrZr              # CuCrZr | Copper | Inconel718 | SS316
  num_channels: 300
  width_mode: fraction_of_pitch # fraction_of_pitch | fixed
  width_fraction: 0.50          # channel width / local pitch
  channel_width: 2.0e-3         # m, when width_mode: fixed
  channel_height: 5.0e-3        # m
  wall_thickness: 0.6e-3        # m, hot-gas-side wall
  roughness: 5.0e-6             # m, absolute
  coolant_fuel_fraction: 1.0    # fraction of the fuel flow through the jacket
  inlet_temperature: 111.66     # K; omit: the fuel's storage temperature, or
                                # the fuel pump's outlet under a cycle
  inlet_pressure: 15.0e6        # Pa; omit: 1.6 p_c, or the fuel pump's
                                # discharge under a cycle
  counterflow: true             # coolant enters at the nozzle end
  x_end_area_ratio: 10.0        # jacket ends where A/At falls to this value
  nusselt_correlation: dittus-boelter  # dittus-boelter | gnielinski
  nusselt_multiplier: 1.0       # correlation uncertainty knob
  hot_gas_model: boundary_layer # boundary_layer (default) | bartz
  hot_gas_multiplier: 1.0       # hot-gas coefficient uncertainty knob
                                # (older name bartz_multiplier still read)
  upstream_wall_temperature: 0  # K, wall upstream of a jacket that starts
                                # downstream of the injector; 0 => first station's
  gas_emissivity: 0.0           # 0 disables the radiation term
  segments: 240                 # jacket march segments
  coolant_wall_limit: 0         # K; warn when the coolant-side wall passes this
                                # (e.g. 728 K, RP-2's measured coking onset); 0 => none
  taper:                        # optional: channel size along the jacket
    x:      [0.00, 0.30, 0.36, 0.45]        # m from the injector face
    height: [5.0e-3, 3.0e-3, 3.0e-3, 5.0e-3] # m; linear between, constant beyond
    width:  [0.50, 0.40, 0.40, 0.50]         # fraction of pitch (m if width_mode: fixed)
  film:                         # optional: wall film from the injector
    fuel_fraction: 0.03         # of the fuel flow (or mass_flow: kg/s)
    slot_height: 0.5e-3         # m, equivalent slot of the film orifices
    x: 0.0                      # m, injection station
    temperature: 0.0            # K; 0 => the jacket's coolant outlet temperature
    coolant: ""                 # property table; empty => the jacket coolant
```

`x_end_area_ratio` matters whenever ε is swept or optimised: without it the
jacket would either stop short of, or run past, the end of the contour when the
nozzle changes length. `nusselt_multiplier` and `hot_gas_multiplier` exist so
that the Monte Carlo campaign can disperse the **correlation** uncertainty
explicitly rather than pretending the correlations are exact.

The wall's conductivity is the material's k(T) fit; there is no YAML key to
override it, but studies can scale its reference value through the
`cooling.wall_conductivity` parameter below.

`hot_gas_model` chooses where the hot-gas film coefficient comes from: the
integral turbulent boundary layer marched from the injector face (the default),
or Bartz's closed form ([`theory.md` §12](theory.md)). Under the boundary layer
the coupled solve takes a few passes, reported in the summary.

`film` injects part of the fuel along the wall (Hatch & Papell, [`theory.md`
§12.4](theory.md)). Its cost in specific impulse is reported as a bracket
between the fully mixed engine (the headline numbers) and an unmixed two-stream
limit. `taper` makes the channel height and/or width vary along the axis; give
`height`, `width` or both, one value per `x`.

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

## `cycle`

```yaml
cycle:
  enabled: true
  type: gas_generator           # gas_generator | staged_combustion | expander
  pump_efficiency: 0.70         # both main pumps (or pump_efficiency_oxidizer / _fuel)
  boost_pump_efficiency: 0.70   # staged combustion's boost stage; default pump_efficiency
  turbine_efficiency: 0.60
  mechanical_efficiency: 0.98   # bearings, seals, gearing
  pump_inlet_pressure: 3.0e5    # Pa, both pumps (or pump_inlet_pressure_oxidizer / _fuel)
  injector_stiffness: 0.20      # main injector dp / p_c; default: the feed block's
  line_loss_fraction: 0.05      # valves and lines, as a fraction of p_c
  turbine_inlet_temperature: 900.0  # K, gas generator or preburner
  fuel_rich: true               # generator / preburner side of stoichiometric
  # gas generator
  gas_generator_pressure: 0     # Pa; 0 => the chamber pressure
  turbine_pressure_ratio: 20.0
  # staged combustion
  preburner_injector_stiffness: 0.15  # dp / p_preburner
  hot_gas_injector_stiffness: 0.10    # (turbine outlet - p_c) / p_c
  preburner_flow_fraction: 1.0        # of the preburner's major propellant
  # expander
  turbine_bypass_fraction: 0.0        # of the jacket flow, round the turbine
```

A pump-fed engine's turbopump power balance ([`theory.md`
§14](theory.md#14-feed-system-and-turbopump-cycles)). A gas generator solves
its flow in closed form and brackets its exhaust's thrust. A staged-combustion
cycle solves its turbine pressure ratio, and an expander its fuel pump's
discharge pressure. Either reports, rather than throws, when no operating point
balances. The jacket, when the fuel cools it, sits between the fuel pump and
everything downstream. With `cooling.inlet_pressure` and
`cooling.inlet_temperature` left out, the engine iterates jacket and cycle
until the jacket starts at the fuel pump's discharge. A configured inlet
pressure above what the injector needs is delivered by the pump, and the
excess is reported as throttled. Propellants with their own real-fluid table
(oxygen, methane, hydrogen) are pumped along it; RP-1, whose table is a
surrogate, is pumped as an incompressible liquid at its measured density.
`ignis_engine --no-cycle` skips the analysis. With both a `feed` and a `cycle`
block, the feed block describes a pressure-fed alternative.

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
    - { parameter: cooling.hot_gas_multiplier, distribution: lognormal,
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
cooling.hot_gas_multiplier      cooling.inlet_temperature
cooling.inlet_pressure          cooling.coolant_fuel_fraction
cooling.roughness               cooling.film_fuel_fraction
cooling.film_slot_height        feed.injector_stiffness
cycle.turbine_inlet_temperature cycle.turbine_efficiency
cycle.pump_efficiency           cycle.turbine_pressure_ratio
cycle.injector_stiffness         performance.kinetic_rate_multiplier
```

`cooling.bartz_multiplier` is accepted as the older name of
`cooling.hot_gas_multiplier`.

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
performance.thrust_sea_level performance.isp_sea_level
performance.isp_ascent     performance.separation_margin performance.kinetic_efficiency
geometry.l_star            geometry.residence_time    geometry.throat_area
geometry.exit_area         geometry.exit_radius       geometry.exit_diameter
geometry.total_length      geometry.divergent_length
cooling.max_wall_temperature cooling.max_heat_flux    cooling.total_heat_load
cooling.pressure_drop      cooling.outlet_temperature cooling.temperature_rise
cooling.energy_balance_residual cooling.flux_residual
cooling.max_coolant_side_wall_temperature
feed.oxidizer_tank_pressure feed.fuel_tank_pressure   feed.total_pump_power
cycle.isp_vacuum_delivered cycle.isp_delivered        cycle.isp_vacuum_dumped
cycle.gas_generator_flow_fraction cycle.fuel_pump_discharge_pressure
cycle.oxidizer_pump_discharge_pressure cycle.pump_power cycle.turbine_pressure_ratio
cycle.max_power_ratio      cycle.feasible
```

The residual metrics are deliberately exposed as sweepable quantities: a sweep
can be asked to show that the element balance stays at 1e-13 across its whole
range, which is a far stronger statement than checking it once.

Running `ignis_sweep --list-parameters` prints this registry with
units, generated from the same table the solver uses.
