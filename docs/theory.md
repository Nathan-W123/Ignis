# Mathematical formulation

Everything Ignis computes is written out here, with the sign and unit
conventions it uses and a citation for every model that is not derived from
first principles. Where a correlation is empirical it says so, in the text and
in the code.

## Contents

1. [Conventions](#1-conventions)
2. [Species thermodynamics](#2-species-thermodynamics)
3. [Chemical equilibrium](#3-chemical-equilibrium)
4. [Equilibrium derivative properties](#4-equilibrium-derivative-properties)
5. [Combustion chamber](#5-combustion-chamber)
6. [Nozzle geometry](#6-nozzle-geometry)
7. [Quasi-one-dimensional flow](#7-quasi-one-dimensional-flow)
8. [Normal shocks and separation](#8-normal-shocks-and-separation)
9. [Performance](#9-performance)
10. [Atmosphere](#10-atmosphere)
11. [Transport properties](#11-transport-properties)
12. [Hot-gas heat transfer and wall conduction](#12-hot-gas-heat-transfer-and-wall-conduction)
13. [Regenerative cooling](#13-regenerative-cooling)
14. [Feed system and turbopump cycles](#14-feed-system-and-turbopump-cycles)
15. [Transient chamber](#15-transient-chamber)
16. [Optimisation](#16-optimisation)
17. [Uncertainty and sensitivity](#17-uncertainty-and-sensitivity)

---

## 1. Conventions

All internal quantities are SI. The one exception is that *angles in
configuration files* are degrees, and every such key says so in its name.

| quantity | unit | quantity | unit |
|---|---|---|---|
| length | m | pressure | Pa |
| mass | kg | temperature | K |
| time | s | amount of substance | mol |
| molar mass | kg/mol | specific enthalpy | J/kg |
| molar enthalpy | J/mol | specific entropy | J/(kg K) |
| mole numbers `n_j` | mol per **kg of mixture** | heat flux | W/m² |

Two conventions matter throughout:

* **The per-kilogram basis.** Composition is carried as mole numbers `n_j` of
  each species *per kilogram of mixture*, the basis NASA CEA uses. By
  construction `Σ_j n_j M_j = 1 kg`, and the mixture molar mass is
  `M = 1 / Σ_j n_j`. The solver reports `Σ_j n_j M_j − 1` as a residual so the
  invariant is measured rather than assumed.
* **Absolute enthalpies.** Every species enthalpy includes its enthalpy of
  formation, and so does every propellant. Combustion therefore needs no heat
  of reaction anywhere: mixing cold propellants and solving for equilibrium at
  constant enthalpy *is* the adiabatic flame calculation.

Sign conventions: heat flux `q` is positive from gas to wall to coolant; the
axial coordinate `x` increases from the injector face (`x = 0`) to the nozzle
exit; thrust is positive opposite to the exhaust velocity.

Reference pressure: **p° = 1 bar = 10⁵ Pa**, which is what NASA TM-4513 and the
NASA CEA database use. This matters — CHEMKIN-derived files are often
interpreted at 1 atm instead, which shifts entropies by `R ln(1.01325) =
0.109 J/(mol K)` and moves a LOX/methane flame temperature by about 2 K. The
verification suite pins the convention explicitly.

---

## 2. Species thermodynamics

Each species carries NASA seven-coefficient polynomials on two temperature
intervals `[T_min, T_mid]` and `[T_mid, T_max]`:

```
cp(T)/R      = a1 + a2 T + a3 T² + a4 T³ + a5 T⁴
h(T)/(R T)   = a1 + a2 T/2 + a3 T²/3 + a4 T³/4 + a5 T⁴/5 + a6/T
s°(T)/R      = a1 ln T + a2 T + a3 T²/2 + a4 T³/3 + a5 T⁴/4 + a7
g°(T)        = h(T) − T s°(T)
```

`a6 R` carries the enthalpy of formation, so `h` is absolute. `s°` is the
standard-state entropy at `p°`; the pressure and mixing dependence is applied
by the mixture layer:

```
s(T, p, n) = Σ_j n_j [ s°_j(T) − R ln(X_j p / p°) ] ,   X_j = n_j / Σ_k n_k
μ_j(T, p, n) = g°_j(T) + R T ln(X_j p / p°)
```

**Data source.** 40 gas-phase species, valid 200–6000 K, taken verbatim from

> B. J. McBride, S. Gordon and M. A. Reno, *Coefficients for Calculating
> Thermodynamic and Transport Properties of Individual Species*, NASA TM-4513,
> October 1993.

through Cantera's machine-readable transcription. `tools/build_thermo_db.py`
records the SHA-256 of the source file in the generated database, and modifies
no coefficient.

**Modelling exclusions.** Gas phase only. No condensed phases (no liquid or
solid water, no graphite), no ionised species and no electrons, so the model
has nothing to say about soot formation, carbon deposition at very fuel-rich
mixture ratios, or the weakly ionised layer at the highest chamber
temperatures. For LOX/CH4 and LOX/H2 at the mixture ratios of interest none of
these is significant; at O/F below about 2 for methane, condensed carbon would
appear in reality and Ignis will not show it.

---

## 3. Chemical equilibrium

### 3.1 The problem

Minimise the total Gibbs energy of the mixture,

```
G(n, T, p) = Σ_j n_j μ_j(T, p, n)
```

over mole numbers `n_j ≥ 0`, subject to conservation of every element:

```
Σ_j a_ij n_j = b_i⁰        i = 1 … E
```

where `a_ij` is the number of atoms of element `i` in species `j` and `b_i⁰`
comes from the propellant composition. For an ideal-gas mixture `G` is convex
on the feasible set, so the minimum is unique.

### 3.2 Stationarity

With Lagrange multipliers `π_i = −λ_i/(R T)` the first-order conditions are

```
μ_j / (R T) = Σ_i a_ij π_i                              (stationarity)
Σ_j a_ij n_j = b_i⁰                                     (element balance)
```

Substituting the ideal-gas chemical potential and rearranging gives the form
the solver actually uses:

```
ln n_j = ln n + Σ_i a_ij π_i − g°_j/(R T) − ln(p/p°)
```

The test suite checks this species by species at the converged point.

### 3.3 Newton iteration

Ignis uses the descent formulation of

> S. Gordon and B. J. McBride, *Computer Program for Calculation of Complex
> Chemical Equilibrium Compositions and Applications*, NASA RP-1311 Part I
> (1994).

Working in `ln n_j` makes every iterate strictly positive, so nonnegativity is
structural and never enforced by clipping. Newton's method applied to the
stationarity and constraint equations, after eliminating the species
corrections analytically, leaves an `(E+1)` system in `(π_i, Δln n)`:

```
Σ_k [Σ_j a_ij a_kj n_j] π_k + [Σ_j a_ij n_j] Δln n
        = b_i⁰ − Σ_j a_ij n_j + Σ_j a_ij n_j (μ_j/RT)
Σ_i [Σ_j a_ij n_j] π_i + [Σ_j n_j − n] Δln n
        = n − Σ_j n_j + Σ_j n_j (μ_j/RT)
```

and the species corrections follow from

```
Δln n_j = Σ_i a_ij π_i + Δln n − μ_j/(R T)
```

The cost per iteration is therefore independent of the number of species: a
4×4 solve for a CHON mixture whether it carries 8 species or 40.

### 3.4 Unknown temperature

For the adiabatic (constant enthalpy, constant pressure) problem, `ln T` joins
the unknowns and the species correction gains a term:

```
Δln n_j = Σ_i a_ij π_i + Δln n + (h_j/RT) Δln T − μ_j/(R T)
```

with the energy equation closing the system:

```
Σ_i [Σ_j a_ij n_j h_j/RT] π_i + [Σ_j n_j h_j/RT] Δln n
   + [Σ_j n_j cp_j/R + Σ_j n_j (h_j/RT)²] Δln T
        = (h⁰ − h)/(R T) + Σ_j n_j (h_j/RT)(μ_j/RT)
```

The isentropic (constant entropy, constant pressure) problem replaces `h_j/RT`
by `S_j/R` in the closing row, where `S_j/R = s°_j/R − ln(X_j) − ln(p/p°)`.
Constant-`u`/constant-`v` and constant-`T`/constant-`v` problems are solved by
an outer iteration around the constant-pressure solver.

### 3.5 Step control and trace species

The Newton step is damped with the RP-1311 limiter

```
λ = min(1, 2 / max(5|Δln n|, 5|Δln T|, max_j |Δln n_j|))
```

where the species maximum runs only over the *major* set (mole fraction above
`10⁻⁸`).

Trace species are handled differently, and this matters. A damped step applied
to a species at a mole fraction of `10⁻³⁰` can fling it to a finite value, after
which it dominates the step limiter and the iteration stalls — this is the
classic failure mode of a naive CEA-style solver at low temperature. Instead,
species below the trace fraction are placed *directly* on their stationarity
value

```
ln(n_j/n) = Σ_i a_ij π_i − g°_j/(R T) − ln(p/p°)
```

capped at a mole fraction of one. That is exactly where the converged solution
puts them, it is bounded, and it lets a species grow back into the major set
naturally when the solution moves. With a trace threshold of `10⁻¹⁴` the
element balance still closes to better than `10⁻¹²`.

### 3.6 Reported residuals

Every solve reports, and the caller can reject on, four independent residuals:

| residual | definition | typical value |
|---|---|---|
| element balance | `max_i \|Σ_j a_ij n_j − b_i⁰\| / max_i b_i⁰` | 10⁻¹⁴ |
| Gibbs optimality | `max_j \|μ_j/RT − Σ_i a_ij π_i\|` over `X_j > 10⁻¹²` | 10⁻¹³ |
| mass normalisation | `Σ_j n_j M_j − 1` | 10⁻¹⁵ |
| state closure | `\|h − h⁰\|/\|h⁰\|` or `\|s − s⁰\|/\|s⁰\|` | 10⁻¹⁴ |

A solve that fails any of them throws rather than returning. Failures name the
element vector, the pressure, and the residuals of every restart attempted.

---

## 4. Equilibrium derivative properties

Nozzle flow needs the speed of sound and the heat capacity *of a mixture whose
composition shifts*, which are not the frozen values. Differentiating the
converged stationarity relation gives two linear systems that share the matrix
of the last Newton step (RP-1311 §2.4). For `(∂/∂ln T)_p`:

```
Σ_k A_ik (∂π_k/∂ln T) + b_i (∂ln n/∂ln T) = −Σ_j a_ij n_j h_j/(R T)
Σ_i b_i (∂π_i/∂ln T)  + 0                 = −Σ_j n_j h_j/(R T)
∂ln n_j/∂ln T = Σ_i a_ij (∂π_i/∂ln T) + ∂ln n/∂ln T + h_j/(R T)
```

and for `(∂/∂ln p)_T` the right-hand sides become `b_i` and `n`. From those:

```
(∂ln v/∂ln T)_p = 1 + ∂ln n/∂ln T
(∂ln v/∂ln p)_T = −1 + ∂ln n/∂ln p
cp   = Σ_j n_j cp_j + (1/T) Σ_j n_j h_j (∂ln n_j/∂ln T)_p
cv   = cp + (p v / T) (∂ln v/∂ln T)_p² / (∂ln v/∂ln p)_T
γ_s  = −(∂ln p/∂ln v)_s = −(cp/cv) / (∂ln v/∂ln p)_T
a    = sqrt(γ_s p v)
```

`γ_s` is the isentropic exponent and is **not** `cp/cv` for a shifting mixture.
For LOX/methane at 5.5 MPa the frozen ratio is 1.198 while `γ_s` is 1.130; using
the wrong one shifts `c*` by several percent. The test suite checks each of
these against a finite difference of independent equilibrium solves.

---

## 5. Combustion chamber

The chamber is an *infinite-area combustor*: propellants burn to equilibrium at
the commanded chamber pressure with no inlet momentum, so the chamber state is a
stagnation state. This is NASA CEA's default assumption and is exact in the
limit of a large contraction ratio. The geometry's contraction ratio affects
residence time, `L*` and heat transfer, never the chamber thermodynamics.

The reactant enthalpy of a mixture at mixture ratio `r = ṁ_ox/ṁ_f` is

```
h⁰ = [ r/(1+r) ] h_ox / M_ox + [ 1/(1+r) ] h_f / M_f      [J/kg]
```

with `h_ox`, `h_f` the absolute molar enthalpies of the propellants in their
storage state. Solving the constant-`h⁰`, constant-`p` equilibrium problem gives
the adiabatic flame temperature and the chamber composition.

Characteristic velocity comes from the *variable-property sonic solution*, not
from a constant-γ formula:

```
c*_ideal = p_c / (ρ* u*)
```

where `(ρ*, u*)` is the state at which the isentropic expansion reaches `M = 1`.
Combustion efficiency is applied separately and visibly:

```
c* = η_c* c*_ideal ,        ṁ = p_c A_t / c*
```

`η_c*` is an empirical input. Ignis never invents a value for it and never
folds it into the thermochemistry.

Residence time and characteristic length use the geometry's chamber volume
(injector face to throat, computed as an exact solid of revolution):

```
t_res = V_c ρ_c / ṁ ,      L* = V_c / A_t
```

---

## 6. Nozzle geometry

The contour follows the classical layout of Huzel & Huang (*Modern Engineering
for Design of Liquid-Propellant Rocket Engines*, AIAA Progress in Astronautics
and Aeronautics Vol. 147, 1992, Ch. 4) with Rao's parabolic bell approximation
(G. V. R. Rao, *Exhaust Nozzle Contour for Optimum Thrust*, Jet Propulsion 28,
377–382, 1958). Six analytic segments, C1 continuous by construction:

1. **Chamber** — cylinder of radius `R_c` over `0 ≤ x ≤ L_c`.
2. **Chamber fillet** — circular arc of radius `R_1`, centred at
   `(L_c, R_c − R_1)`:
   `x = L_c + R_1 sin θ`, `r = R_c − R_1(1 − cos θ)`, `θ: 0 → β`.
3. **Converging cone** — straight, half-angle `β`.
4. **Throat upstream arc** — radius `R_u` (default `1.5 R_t`), centred at
   `(x_t, R_t + R_u)`: `x = x_t − R_u sin θ`, `r = R_t + R_u(1 − cos θ)`.
5. **Throat downstream arc** — radius `R_d` (default `0.382 R_t`, Rao):
   `x = x_t + R_d sin θ`, `r = R_t + R_d(1 − cos θ)`, `θ: 0 → θ_n`.
6. **Divergent** — either a straight cone at `θ_cone`, or a quadratic Bézier
   from `N` (slope `tan θ_n`) to the exit `E` (slope `tan θ_e`) with control
   point at the intersection of the two tangents:
   `P(t) = (1−t)² N + 2t(1−t) Q + t² E`.

The bell length is set as a fraction of the equivalent 15° conical length,

```
L_15 = [ R_t (√ε − 1) + R_d (sec 15° − 1) ] / tan 15°
```

Geometric validity is checked, not assumed: the converging section must
physically fit (`r_A > r_B`), the Bézier control point must lie inside the
divergent, the throat must be the unique minimum of the sampled contour, and
the area must be monotone on each side of it. Each failure names the offending
input and the value that would fix it.

**This is not a manufacturing contour.** `θ_n` and `θ_e` are user inputs, not
values read from Rao's charts; a true optimum bell requires the method of
characteristics. There is no injector-face detail, no wall thickness taper and
no throat insert.

---

## 7. Quasi-one-dimensional flow

Steady, adiabatic, inviscid, one-dimensional flow from the chamber stagnation
state:

```
h(p) + u²/2 = h⁰        (energy)
s(p)        = s⁰        (isentropic)
ρ u A       = ṁ         (mass)
```

Ignis marches on **static pressure** rather than area, because that makes the
state determination direct:

* **Frozen composition** — the chamber mole numbers are held fixed and `T` is
  found from `s(n_c, T, p) = s⁰`.
* **Shifting equilibrium** — a constant-entropy, constant-pressure Gibbs
  minimisation is solved at each station.

Then `u = sqrt(2(h⁰ − h))`, `M = u/a` with `a` from `γ_s`, and the local area
ratio follows from mass conservation:

```
A/A_t = G* / G ,      G = ρ u
```

No ratio of specific heats is ever assumed constant. The throat is located by
solving `M(p) = 1` with bracketed Illinois regula falsi. Each commanded area
ratio is inverted on the appropriate branch of `G(p)`, which rises from zero at
`p = p⁰` to a maximum at the sonic point and falls again. Every function
evaluation is a Gibbs minimisation, so the inversion uses the slope each solved
state already carries. Along the isentrope `dp = −ρu du` and `dp = a² dρ`, so

```
d ln G / d ln p = p/(ρa²) − p/(ρu²) = (M² − 1) / (γ_s M²)
```

That is a safeguarded Newton iteration in `ln p`. It starts from the
perfect-gas estimate at the throat's `γ_s`, and any step that would leave the
bracket is replaced by a bisection. It stops when `G` is within 10⁻¹¹ of the
target, usually after three or four Gibbs solves. The bracketing secant search
it replaced needed about 35.

**Range limit.** The expansion cannot proceed past the point where the static
temperature reaches the 200 K floor of the polynomial fits. Ignis reports the
largest reachable area ratio and the pressure at which it occurs, and refuses
the request; it does not extrapolate.

---

## 8. Normal shocks and separation

A normal shock is treated as a discontinuity of zero thickness with frozen
composition (the residence time inside a shock is orders of magnitude shorter
than any chemical time). The jump is solved from the exact conservation
equations using the mixture's own caloric equation of state h(T), not from
constant-γ relations.

> **Terminology.** This is a *variable-property equilibrium-gas* shock, not a
> "real-gas" shock. The thermal equation of state stays ideal, p = ρRT, with
> Z = 1 throughout; what varies is the caloric behaviour — c<sub>p</sub>(T),
> γ(T) and the mean molar mass. Ignis implements no compressibility factor, no
> fugacity and no non-ideal mixing rule, and at the densities involved (§1 of
> [`limitations.md`](limitations.md)) it does not need to.

Parameterising by `v = u₂/u₁ = ρ₁/ρ₂`:

```
p₂ = p₁ + ρ₁ u₁² (1 − v)
h₂ = h₁ + u₁² (1 − v²)/2
T₂ = p₂ v / (ρ₁ R)
```

and requiring `h(T₂) = h₂` gives a one-dimensional root find. `v = 1` is the
trivial root; the physical root is the other sign change. Entropy is checked to
increase.

Which regime applies is decided by comparing the ambient pressure with three
exit pressures — the shock-free supersonic value, the value behind a normal
shock standing at the exit plane, and the subsonic branch at the same area
ratio:

| condition | regime |
|---|---|
| `p_a < p_e,sup` | under-expanded |
| `p_a ≈ p_e,sup` | ideally expanded |
| `p_e,sup < p_a < p_e,shock@exit` | over-expanded; compression outside through oblique shocks (not modelled) |
| `p_e,shock@exit ≤ p_a < p_e,sub` | normal shock inside the divergent; its position is found by bisection |
| otherwise | reported as unresolved — **no supersonic solution is returned** |

**Separation** is viscous and cannot be predicted by an inviscid quasi-1D
model. Two published criteria are offered, both clearly labelled empirical, and
both evaluated purely as diagnostics that never alter the inviscid solution:

* Summerfield: separation where the wall pressure falls below `0.4 p_a`
  (M. Summerfield, C. Foster and W. Swan, *Flow Separation in Overexpanded
  Supersonic Exhaust Nozzles*, Jet Propulsion 24, 319–321, 1954);
* Schmucker: `p_sep/p_a = (1.88 M − 1)^(−0.64)`
  (R. H. Schmucker, NASA TM-77396, 1984).

---

## 9. Performance

```
F_momentum = ṁ u_e
F_pressure = (p_e − p_a) A_e
F_ideal    = F_momentum + F_pressure
F          = λ_div η_nozzle F_momentum + F_pressure
I_sp       = F / (ṁ g₀) ,     c_eff = F / ṁ
C_F        = F / (p_c A_t) ,  c* C_F = I_sp g₀
```

with `g₀ = 9.80665 m/s²` used only to express an effective exhaust velocity in
seconds. Ideal and corrected values are reported side by side, and the loss
factors are listed explicitly in every output.

The divergence efficiency, when it is computed automatically, is the standard
geometric correction `λ = (1 + cos θ_e)/2` — an **empirical** factor that stands
in for the radial component of the exit momentum.

Vacuum specific impulse is always reported from the shock-free solution, which
is the conventional definition. With the boundary layer on (§9.3) the mass
flow, the exit state and the momentum term carry its corrections; the formulas
above are then the inviscid ones that "ideal" still means.

### 9.1 The ambient-pressure family

For a **full-flowing** nozzle every term above is fixed by the converged exit
state except the pressure term, which is linear in the ambient pressure:

```
F(p_a) = F_vac − p_a A_e ,     F_vac = λ_div η_nozzle ṁ u_e + p_e A_e
```

So the sea-level thrust, the vacuum thrust and any weighted average over a
trajectory all follow from one solve. Ignis reports

```
F_sl   = F(101 325 Pa) ,                     I_sp,sl = F_sl / (ṁ g₀)
I_sp,ascent = Σ w_i F(p_a(z_i)) / (Σ w_i ṁ g₀)
```

over a **declared** list of altitudes `z_i` and time weights `w_i`
(`performance.ascent_profile`). Because F is linear in p_a, this weighted mean
is exact — it equals `F` evaluated at the weighted-mean ambient pressure — and
is not a quadrature approximation. The linearity is verified against the
explicit sum in the test suite.

This is a **mission weighting supplied as an input**, not a trajectory
simulation: Ignis has no vehicle model. It exists because optimising a
single-altitude specific impulse drives the expansion ratio to whatever bound
it is given, whereas a trajectory-weighted objective has an interior optimum.

### 9.2 Separation margin

The Summerfield and Schmucker criteria give a wall pressure `p_sep(M, p_a)`
below which the flow is predicted to detach. Wall pressure falls and Mach
number rises monotonically through the divergent, so `p_wall − p_sep` is
monotone and its smallest value is at the exit plane. Ignis reports the
normalised exit-plane margin

```
m_sep = (p_e − p_sep(M_e, p_a)) / p_a
```

which is positive when the flow is attached to the lip, negative when the
criterion predicts separation inside the nozzle, and usable directly as an
optimisation constraint. It is a **diagnostic on an empirical criterion**: the
inviscid solution is not modified when it goes negative, so a separated
nozzle's reported thrust is optimistic.

### 9.3 Boundary-layer losses

The quasi-1D core is inviscid; the real flow carries a wall boundary layer
whose displacement thickness `δ*` narrows the passage the core sees and whose
momentum thickness `θ` is momentum the wall has taken out of it. By default
(`performance.boundary_layer_losses: true`) Ignis marches the layer of §12.1
from the injector face to the exit plane — with the cooling jacket's hot-wall
temperature where a jacket is solved and `performance.uncooled_wall_temperature`
(default 1000 K, an assumption) elsewhere — and applies the three standard
corrections of the boundary-layer method (Elliott, Bartz & Silver, JPL TR
32-387, Section V; the same structure as the JANNAF boundary-layer module):

```
A*_eff = π (r_t − δ*_t)² ,       C_d = A*_eff / A*,      ṁ = C_d p_c A* / c*
ε_core = π (r_e − δ*_e)² / A*_eff          (the core's exit area ratio)
F      = λ_div η_nozzle (ṁ u_e − 2π r_e ρ_e u_e² θ_e) + (p_e − p_a) A_e
```

with the edge state (`u_e`, `ρ_e`, `p_e`) that of the core at `ε_core`. The
pressure term acts on the whole geometric exit. The momentum deficit uses the
thin-layer form `2π r θ`; on the engines here the exit layer is 1–3 % of the
exit radius. "Ideal" values stay inviscid at the geometric area ratio, so they
remain comparable with CEA, and the inviscid corrected values are reported
alongside for the difference.

A cold wall makes the layer denser and `δ*` smaller (on an M1-like test engine a 150 K wall
gives a third of the throat displacement a 1500 K wall does); in the limit it
can go negative and `C_d` exceed one. Nothing clips that.

On the shipped M1 the layer gives `C_d` = 0.99865, a core that expands to 19.6
instead of the geometric 20, and a vacuum-Isp loss of 1.27 % (4.5 s); on the H1,
1.31 % (6.0 s). These are not validated against a measured loss: neither §5b
experiment measured thrust. The coolant flow of the jacket is set from the
inviscid mass flow, before the layer moves it by `1 − C_d` (0.13 % on the M1).
If an internal normal shock is found, the momentum deficit is not applied — the
layer was marched against the attached flow — but `C_d` still is.

### 9.4 Finite-rate recombination

Shifting equilibrium assumes the recombination of H, O, OH and CO keeps pace
with the expansion to the exit; frozen flow assumes it stops at the throat. The
real gas does neither. With `performance.kinetics` enabled, Ignis integrates the
species along the divergent nozzle at finite rate (`ignis/kinetics/`), the
one-dimensional kinetics (ODK) calculation of the JANNAF methodology.

**Mechanism.** The reactions of GRI-Mech 3.0 (Smith et al.,
http://combustion.berkeley.edu/gri-mech/) among H, H<sub>2</sub>, O, O<sub>2</sub>,
OH, H<sub>2</sub>O, HO<sub>2</sub>, H<sub>2</sub>O<sub>2</sub>, CO, CO<sub>2</sub> and HCO:
41 reactions, extracted by `tools/build_kinetics.py` into
`data/kinetics/gri30_nozzle.yaml`. They are elementary, three-body and
Lindemann/Troe fall-off reactions. Only forward rates are taken. Each reverse
rate follows by detailed balance from Ignis's own species data,

```
k_r = k_f / K_c,     K_c = exp(−ΔG°/RT) (p°/RT)^Δν
```

so the kinetics relax to exactly the equilibrium the solver of §3 computes.
Reactions naming species the run's database does not carry (CO in a
hydrogen engine, N<sub>2</sub> or Ar colliders) are left out and listed.
Species outside the mechanism ride along at their start values.

**Equations.** Steady quasi-1D flow through the prescribed area `A(x)`, with
the composition `n_k` in mol/kg:

```
ρ u A = ṁ,     h(T, n) + u²/2 = h₀,     dn_k/dx = ω̇_k / (ρ u)
du/dx = u (−A′/A + Ψ) / (1 − M_f²),    Ψ = Σ dn_k / Σ n_k − Σ H_k dn_k / (c_p,f T)
```

`M_f` is the frozen Mach number and `Ψ` the recombination's heat release and
mole change, which accelerates the gas. `T` is solved from the energy equation
at every evaluation, so energy is conserved exactly, and `ρ` from mass.

**Start.** The momentum equation is singular at `M_f = 1`, just downstream of
the equilibrium throat (the frozen sound speed exceeds the equilibrium one). The
march starts from the shifting-equilibrium solution where `M_f` reaches
`start_frozen_mach` (1.10). Upstream of that point, at throat temperature and
pressure, the gas is close to equilibrium. The result moves by 3 × 10⁻⁶ in
specific impulse between starts at 1.05 and 1.20.

**Integration.** The species equations are stiff: the chemistry relaxes in
nanoseconds near the throat while the gas crosses the nozzle in a fraction of a
millisecond. Each step is linearly implicit Euler with one Richardson
extrapolation, which is second order and L-stable, with a finite-difference
Jacobian. Linear invariants pass through the implicit solve unchanged, so every
element is conserved to round-off. The difference between one full step and two
half steps sets the step size.

**What it changes.** The march gives the inviscid vacuum impulse at finite rate
and, from the same geometry, the shifting-equilibrium and frozen-from-the-start
limits. The kinetic efficiency `η_kin = I_sp,vac(finite rate) / I_sp,vac(shifting)`
is the JANNAF ODK/ODE ratio. With `apply: true` (the default) it scales the
delivered vacuum thrust, and every ambient figure follows as
`F = η_kin F_vac − p_a A_e`. The boundary-layer loss of §9.3 is computed on the
shifting core and reported separately. Ideal values stay shifting equilibrium,
comparable with CEA. The kinetic efficiency is 0.9953 on the M1 (1.6 s of
vacuum impulse), 0.9965 on the H1 and 0.9979 on the K1. The march recovers
89 %, 86 % and 94 % of the impulse between frozen and shifting.

**Limits and checks.** With every rate multiplied by zero the march reproduces a
frozen isentropic expansion from its start to 2 × 10⁻⁶. With the rates
multiplied by a thousand it comes within 10⁻⁵ of shifting equilibrium, and the
gap closes roughly tenfold per decade of rate. Forward rate constants agree
with Cantera's to 3 × 10⁻¹². The whole march is compared with Cantera on
identical data in [`validation.md` §2.5](validation.md).

---

## 10. Atmosphere

The U.S. Standard Atmosphere 1976 (NOAA-S/T 76-1562), seven layers to 84.852 km
geopotential altitude:

```
L_b ≠ 0:  p = p_b [ T_b / (T_b + L_b (H − H_b)) ]^( g₀ M₀ / (R* L_b) )
L_b = 0:  p = p_b exp[ −g₀ M₀ (H − H_b) / (R* T_b) ]
H = r₀ z / (r₀ + z)
```

using the constants *defined by that document* (`R* = 8.31432 J/(mol K)`,
`M₀ = 28.9644 kg/kmol`, `g₀ = 9.80665 m/s²`, `r₀ = 6 356 766 m`) rather than
modern CODATA values, so the published table is reproduced. Layer base
pressures are integrated upward from 101 325 Pa instead of being transcribed,
which removes any possibility of a copying error; the model returns 0.37338 Pa
at 86 km, matching the standard. Above 86 km it switches to an exponential
extrapolation and flags the state.

---

## 11. Transport properties

Pure-species viscosity from Chapman-Enskog first-order kinetic theory
(Hirschfelder, Curtiss & Bird, *Molecular Theory of Gases and Liquids*, Wiley
1954, Eq. 8.2-18):

```
μ_k = (5/16) sqrt(π m_k k_B T) / (π σ_k² Ω^(2,2)*(T*)) ,   T* = k_B T / ε_k
```

with the collision integral from the correlation of Neufeld, Janzen & Aziz
(*J. Chem. Phys.* 57, 1100, 1972),

```
Ω^(2,2)* = 1.16145 T*^−0.14874 + 0.52487 e^(−0.77320 T*) + 2.16178 e^(−2.43787 T*)
```

valid for `0.3 ≤ T* ≤ 100` with a stated maximum error of 0.064 %. Polar
species use Brokaw's correction (Poling, Prausnitz & O'Connell, *The Properties
of Gases and Liquids*, 5th ed., Eq. 9-4.3):

```
Ω^(2,2)*_polar = Ω^(2,2)*_LJ(T*) + 0.2 δ*² / T* ,
δ* = μ_D² / (2 · 4πε₀ · ε σ³)
```

Thermal conductivity uses the modified Eucken correlation (Poling et al.,
Eq. 10-3.3), which reduces exactly to the monatomic Eucken value:

```
λ_k = (μ_k / M_k)(1.32 cv_k + 1.77 R)
```

Mixtures use Wilke's rule (C. R. Wilke, *J. Chem. Phys.* 18, 517, 1950) for
viscosity and the Wassiljewa form with the same interaction coefficients for
conductivity:

```
μ_mix = Σ_k X_k μ_k / Σ_j X_j Φ_kj
Φ_kj  = [1 + sqrt(μ_k/μ_j) (M_j/M_k)^(1/4)]² / sqrt(8(1 + M_k/M_j))
```

Lennard-Jones parameters come from the GRI-Mech 3.0 transport database. The
model is a **dilute-gas** theory with no density correction; species without LJ
data are excluded from the average and the covered mole fraction is reported.

---

## 12. Hot-gas heat transfer and wall conduction

Two models supply the hot-gas film coefficient `h_g`, selected by
`cooling.hot_gas_model`. The default is an integral turbulent boundary layer
(§12.1); Bartz's closed form (§12.2) remains available. Both are measured
against the same two experiments in [`validation.md` §5b](validation.md).

### 12.1 Integral turbulent boundary layer (default)

Bartz fixes the boundary layer's thickness by assumption, scaling everything
on the throat diameter, so it cannot know how far the layer has grown — and
the measured error has exactly that shape (§12.2). Marching the layer from
where it starts carries its history. The integral equations are those of
Elliott, Bartz & Silver (JPL TR 32-387, 1963; Bartz, *Advances in Heat
Transfer* 2, 1965), for the momentum thickness `θ` and the energy (enthalpy)
thickness `Δ`:

```
dθ/dx = Cf/2 − θ [ (2 + H − M_e²) d ln u_e/dx + d ln r/dx ]
dΔ/dx = St   − Δ [ d ln(G_e c_p* (T_aw − T_w))/dx + d ln r/dx ]
```

with `G_e = ρ_e u_e`, and the edge density gradient eliminated by the
isentropic relation `d ln ρ_e = −M_e² d ln u_e`. The closure was fixed from the
textbook before any comparison with data and has not been adjusted since:

* **Profiles.** `u/u_e = (y/δ)^(1/7)`, temperature by Crocco–Busemann,
  `T = T_w + (T_aw − T_w)(u/u_e) − (T_aw − T_e)(u/u_e)²`. The shape factor
  `H = δ*/θ` follows by quadrature (taken in `u`, which removes the
  `y^(1/7)` cusp at the wall; it reproduces the incompressible 9/7 to 4 × 10⁻⁶).
* **Reference state.** Properties at Eckert's
  `T* = T_e + 0.5 (T_w − T_e) + 0.22 (T_aw − T_e)`, the local static pressure
  and the local frozen composition.
* **Skin friction** (1/7-law, Blasius) and **Stanton number** (Colburn analogy
  on the energy thickness):

```
Cf/2 = 0.0128 (ρ*/ρ_e) (G_e θ / μ*)^(−1/4)
St   = 0.0128 (ρ*/ρ_e) (G_e Δ / μ*)^(−1/4) Pr*^(−2/3)
h_g  = St G_e c_p*
```

On a flat plate this is the textbook `θ/x = 0.0366 Re_x^(−1/5)`, and the
unit tests check the march against that and against the exact solutions that
exist for an accelerating stream, a widening duct and a varying wall
temperature (with properties held fixed, `θ^(5/4)` obeys a linear equation).

**Where the layer starts.** It is turbulent from the injector face, as in
JPL TR 32-387 ("growth from zero thickness") and in Wang & Luong's SSME
analysis (*J. Thermophysics* 8(3), 1994: "the turbulent boundary layer started
from the injector faceplate"). It cannot start at zero — `Cf` and `St` go as
thickness<sup>−1/4</sup> — so it starts as the thinnest layer that can be
turbulent at all: momentum-thickness Reynolds number 320, Preston's minimum
(*J. Fluid Mech.* 3, 373, 1958), with `Δ = θ`. On the engines here that is a
few hundredths of a millimetre with a virtual origin a few millimetres
upstream; the layer forgets it within millimetres, and a start ten times
thinner moves the downstream `h_g` by 0.1–0.2 %.

**Consequence: a peak at the injector face.** A layer that starts there has
its highest coefficient there. On the M1 `h_g` near the face is about twice
Bartz's, crosses it about 15 cm downstream, and is about 0.8 of it at the
throat. Without a film the first few centimetres of the M1's jacket run near
1000 K, which is why the M1 now carries one (§12.4). Real chambers usually
measure low flux near the faceplate, which Wang & Luong attribute to film
cooling; a chamber with *no* film and a fully burnt gas at the wall is the
case this model describes, and it is probably conservative there — the
injector's own recirculation and incomplete combustion are not modelled.

**Validity.** The closure is turbulent everywhere. The march reports where the
acceleration parameter `K = (ν/u²) du/dx` passes 3 × 10⁻⁶, above which a layer
starts to laminarise (Moretti & Kays 1965), and where `Re_θ` falls below 320.
Neither condition is modelled; both are flagged.

**Coupling.** `h_g` now depends on the wall temperature both locally (the
reference state) and through the layer's history (the energy integral's
driving difference), so the cooling solve is repeated — march the layer with
the current wall temperatures, solve the jacket, update — until the hot-wall
temperature moves less than 0.5 K (three or four passes on the shipped
designs). Within each pass every station is balanced exactly at its trial wall
temperature; only the layer's thickness lags by one pass.

**Resolution.** The coefficient varies on the scale of the virtual-origin
length near the face, so where the jacket starts at the injector the first
segment is 1/N of the uniform width and the segments grow geometrically by 1.5
until they reach it. The peak wall temperature then converges with segment
count (0.6 %, 0.14 %, 0.03 %, 0.007 % per doubling from 50 segments on the test
engine), as the heat load does.

### 12.2 Bartz (option)

D. R. Bartz, *A Simple Equation for Rapid Estimation of Rocket Nozzle
Convective Heat Transfer Coefficients*, Jet Propulsion 27, 49–51, 1957:

```
h_g = (0.026 / D_t^0.2) (μ^0.2 cp / Pr^0.6) (p_c / c*)^0.8 (D_t/R_curv)^0.1 (A_t/A)^0.9 σ

σ = 1 / { [0.5 (T_wg/T_c)(1 + (γ−1)M²/2) + 0.5]^0.68 [1 + (γ−1)M²/2]^0.12 }
```

with `μ`, `cp` and `Pr` evaluated at chamber stagnation conditions, as in
Bartz's derivation. Written in SI the correlation is dimensionally consistent
and needs no conversion factor; the test suite checks this.

Measured against two 1965 experiments ([`validation.md` §5b](validation.md)),
it is within 1 % of the measured constant in a rocket chamber and about 70 %
high at the throat — an error with a *shape* along the engine, which is what
the boundary layer above was brought in to address.

### 12.3 Recovery temperature, conduction and radiation

Recovery temperature, turbulent:

```
T_aw = T (1 + r (γ−1) M² / 2) ,     r = Pr^(1/3)
```

Radial conduction through a cylindrical shell with a temperature-dependent
conductivity evaluated at the mean wall temperature:

```
q_inner = k(T̄) (T_wg − T_wc) / [ r_i ln(r_o/r_i) ]
```

Radiation, when enabled, is kept strictly separate from convection:

```
q_rad = σ_SB ε_gas ε_wall (T_gas⁴ − T_wg⁴)
```

`ε_gas` is a user input with no default. A gray-gas emissivity is a crude
stand-in for the H₂O/CO₂ band structure; quantitative radiative flux needs a
band model.

### 12.4 Film cooling

A film of coolant injected along the wall at `film.x` (default the injector
face) lowers the temperature the gas drives the wall towards. The
effectiveness is Hatch & Papell's (NASA TN D-130, 1959, their eq. 12):

```
η = (T_aw − t_w) / (T_aw − t_c)                                   (adiabatic wall)
ln η = −[ H − 0.04 ] (S V_g / α_c)^0.125 f(V_g/V_c) ,   η = 1 for H < 0.04
f(r) = 1 + 0.4 atan(r − 1)          for r = V_g/V_c ≥ 1
     = (1/r)^(1.5 (1/r − 1))         for r ≤ 1
H    = ∫ h L dx / (ẇ c_p)_c ,   h = 0.0265 (k/D) Re^0.8 Pr^0.3 at (t_g + t_c)/2
```

`S` is the slot height, `L` the circumference, `ẇ` the film flow, `t_c`,
`α_c` and `V_c` the coolant's temperature, diffusivity and velocity at the
slot, `V_g` the gas velocity there. The report takes `h` and `L` constant;
in a nozzle they are not, and `H` is integrated from the slot — the form the
report's own derivation (its eq. 2) has before that simplification. That
generalisation is Ignis's, not the report's. The wall then sees

```
T_drive = T_aw − η (T_aw − t_c) ,     q_conv = h_g (T_drive − T_wg)
```

with `h_g` left as the hot-gas model gives it for the unfilmed wall — the
usual superposition, in which the film changes the driving temperature and not
the layer's conductance. The film temperature is an input or, by default, the
jacket's coolant outlet temperature (a fuel film drawn from the injector
manifold after the jacket), which is iterated with the passes of §12.1.

Against the report's own helium data (Table I, 27 runs) the implementation
reproduces the wall temperature to a mean 4.0 % inside the velocity-ratio range
the report fitted (0.45 – 33.3), where it claims about ±5 %; below that range
eq. (11) is extrapolated, the function reaches 300, and the error is 40 %
([`validation.md` §5c](validation.md)). Outside the fitted range the ratio used
in `f` is held at the nearer limit and a warning is raised. A rocket film —
dense supercritical fuel beside a burning 3500 K gas — is outside every range
the correlation was fitted over; it is the classical first estimate, not a
validated rocket film model.

**Its cost in specific impulse is reported as a bracket.** The headline
performance assumes the film mixes and burns completely (the overall mixture
ratio). The other limit keeps it as a separate unburnt stream: the core burns
at the mixture ratio the film leaves it, the film expands as a calorically
perfect gas from its injection state to the core's exit pressure, and the two
are mass-weighted in vacuum. On the M1 (3 % of the fuel) the bracket is
344.83 s to 342.61 s, at most 0.64 %.

**This is an engineering estimate, not a conjugate CFD solution.** There is no
injector streak, no soot layer, no transpiration cooling and no combustion in
the film. `hot_gas_multiplier` (formerly `bartz_multiplier`, still accepted)
makes the coefficient's uncertainty explicit and is one of the dispersed Monte
Carlo inputs.

---

## 13. Regenerative cooling

Each axial segment couples the three resistances through one unknown, the
hot-side wall temperature:

```
q_gas  = h_g (T_drive − T_wg) + q_rad(T_wg)
T_wc   = T_wg − q_gas R_wall(T̄)
q_cool = h_c (A_cool,eff / A_gas) (T_wc − T_bulk)
```

where `T_drive` is the recovery temperature, or less under a film (§12.4), and
`q_gas(T_wg) = q_cool(T_wg)` is solved to tolerance. Both sides are genuinely
nonlinear in `T_wg` — `h_g` through its reference state, `h_c` through the
coolant properties — so the balance is solved with a bracketed Illinois
iteration. The residual falls strictly with `T_wg`, so a bracket from a kelvin
below both the coolant and the driving temperature up to the recovery
temperature always holds exactly one root; that matters beside a film drawn
from the jacket, which can sit within millikelvin of the coolant behind the
wall. The residual `|q_gas − q_cool| / max(|q_gas|, 1 kW/m²)` is reported per
station.

The land between channels is a straight rectangular fin of height equal to the
channel height:

```
m = sqrt(2 h_c / (k_wall t_land)) ,   η_fin = tanh(m H)/(m H)
A_cool,eff per unit length = w_channel + 2 H η_fin
```

The outer closeout is taken as adiabatic, which is conservative.

The coolant is marched in **enthalpy**, not temperature, so the energy balance
closes identically however violently `cp` varies near the pseudo-critical line:

```
h_out = h_in + Q_segment / ṁ_coolant ,   T = T(h, p)
```

Pressure carries both friction and momentum terms, the latter significant when
a supercritical coolant expands by a factor of five along the jacket:

```
dp = −f (dx/D_h) (G²/2ρ) − G du ,     u = G/ρ
```

which for a constant channel is the familiar `−G² d(1/ρ)`.

**Tapered channels.** Channel height and width may be given as tables in `x`
(`cooling.taper`), linear between points and constant beyond them; the width is
a fraction of the local pitch or, in fixed-width mode, a length. The momentum
term above then also carries the acceleration of the coolant through a
narrowing channel.

Correlations, all labelled:

| quantity | correlation | reference |
|---|---|---|
| friction, turbulent | Colebrook-White | Colebrook, *J. Inst. Civ. Eng.* 11, 133 (1939) |
| friction, laminar | `f = 64/Re` | — |
| Nusselt, default | `Nu = 0.023 Re^0.8 Pr^0.4` | Dittus & Boelter, Univ. Calif. Publ. Eng. 2, 443 (1930) |
| Nusselt, option | Gnielinski | Gnielinski, *Int. Chem. Eng.* 16, 359 (1976) |
| Nusselt, laminar | `Nu = 4.36` | fully developed, constant heat flux |

The two turbulent correlations differ by **10.8 %** in peak wall temperature
for the shipped design (Dittus-Boelter 797.4 K against Gnielinski 711.3 K,
measured by the test suite); that spread is the real coolant-side uncertainty
and is why `nusselt_multiplier` is a dispersed Monte Carlo input.

**Coolant properties** come from tabulated reference equations of state
(Setzmann & Wagner 1991 for methane, Leachman et al. 2009 for hydrogen,
Schmidt & Wagner 1985 for oxygen) rather than a cubic equation of state. A
Peng-Robinson implementation is included so that the choice is justified by a
measurement rather than an assertion: against the same reference tables its
density is **2.2 %** out above 250 K and **12.5 %** out below, in exactly the
cold dense region a regenerative jacket inlet sits in. The test requires the
dense-region error to exceed the dilute-region one, so the comparison would
fail if the cubic were ever silently swapped in.

**Failure modes are reported, never absorbed:** coolant boiling (a sub-critical
liquid crossing its saturation line), a state outside the tabulated range
(including methane's melting line, which rises into the grid above about
15 MPa), a wall above the material limit, conductivity evaluated outside its
fitted range, channels that leave no land, and a coolant pressure driven to
zero.

---

## 14. Feed system and turbopump cycles

### 14.1 Pressure-fed supply

Pressure-fed, per propellant leg:

```
p_tank = p_c + Δp_injector + Δp_lines + Δp_dynamic
Δp_injector = ṁ² / [2 ρ (C_d A)²]
Δp_lines    = [f L/D + Σ K] ρ V²/2
Δp_dynamic  = ρ V²/2
```

Two questions can be asked: *sizing* (given `p_c` and `ṁ`, what tank pressure
and injector area are needed?) and *capability* (given fixed tank pressures and
injector areas, what flow is delivered?). The second is a bracketed solve
because every loss term rises monotonically with flow.

Injector stiffness `Δp_injector / p_c` is reported and checked against a
configurable floor, because a soft injector invites combustion instability.
**Ignis does not model combustion stability.**

The feed block also reports `P = ṁ Δp / (ρ η)` for the same duty, as a
comparison figure. A pump-fed engine is the `cycle` block's job (14.2–14.3);
when both blocks are present the feed block describes a pressure-fed
alternative. Pressurant mass is the isothermal ideal-gas bound
`m = p_tank V M / (R T)`, stated as a bound and not a prediction.

### 14.2 Pumps, turbines, generators

`ignis/cycle/Turbomachinery.hpp`.

**Pump.** The isentropic work per kg is `∫ dp / ρ` along the isentrope. Where
the propellant has its own real-fluid table (`data/coolants`: oxygen, methane,
hydrogen) the integral is marched along the table, 400 midpoint steps from the
inlet to the discharge pressure. On an isentrope `dh = v dp`, so marching `h`
with `dp/ρ(T, p)` and recovering `T` from `(h, p)` follows the isentrope
without the table carrying entropy. Liquid hydrogen is compressible enough
that this matters: at the RS-25's high-pressure fuel pump (1.7 → 41 MPa) its
density rises 17 %, and `Δp/ρ_inlet` overstates the work by 13 %. Where the
table is a surrogate (RP-1's n-dodecane, 7 % lighter than RP-1) the pump would
inherit the surrogate's density error, so the liquid is taken as
incompressible at the propellant's measured density instead; a kerosene's
compressibility over a few tens of MPa is under one percent. The pump's actual
enthalpy rise is the isentropic one over its efficiency, and its outlet
temperature follows from the table.

**Gas turbine.** Generator or preburner gas of composition `n` expands at
frozen composition: `s(T₂s, p₂) = s(T₁, p₁)` on the species
thermodynamics, work `η_T [h(T₁) − h(T₂s)]` per kg.

**Fluid turbine.** An expander's heated fuel expands along its real-fluid
table, the pump's march in reverse.

**Generator or preburner.** The mixture ratio whose adiabatic flame
temperature is the requested turbine inlet temperature, on the fuel-rich or
the oxidiser-rich side of stoichiometric. It is solved as the energy balance
it is,

```
h_reactants(MR) = h_eq(MR; T_turbine, p)
```

with fixed-temperature equilibria, bisected in `log MR`. Comparing adiabatic
flame temperatures with the target instead would ask for adiabatic states
colder than the species data reach: cryogenic hydrogen far from
stoichiometric burns below their 200 K floor. The reactants carry whatever the
pumps and the jacket have put into them, which matters for hydrogen
(14.3). The gas is chemical equilibrium at the turbine inlet temperature. For
hydrogen that is a good description. A fuel-rich hydrocarbon generator makes
soot and cracked fuel far from equilibrium, and Ignis carries no condensed
carbon, so for methane or kerosene the composition is approximate. The run
says so.

### 14.3 Cycle closures

`ignis/cycle/Cycle.hpp`. Each closure finds the operating point at which

```
η_mech P_turbine = P_pumps
```

with one equivalent turbine on all the pumps. A twin-shaft engine balances
each shaft separately; this checks their sum.

**Pressures.** A pump's discharge is what lies downstream of it: the main
injector `p_c (1 + s_inj)`, valves and lines `f_line p_c`, and on the fuel
side the jacket's pressure drop from the cooling solve. The jacket follows
the fuel pump directly, so the fuel pump's discharge *is* the jacket's inlet
pressure, and its outlet temperature the jacket's inlet temperature. Unless
`cooling.inlet_pressure` and `cooling.inlet_temperature` are configured, the
engine iterates the jacket and the cycle until they agree (0.1 % in pressure,
0.05 K). When the inlet pressure is configured above what the injector needs,
the pump delivers it and the excess is reported as throttled.

**Gas generator (open).** The generator is fed from the main pump discharges
and runs at the chamber pressure unless told otherwise; its turbine pressure
ratio is an input. The pumps also lift the generator's propellant, so with `w`
the work per kg each pump does, `x` the generator's propellant mass fractions
and `Δh_T` the turbine's work per kg,

```
ṁ_gg = (ṁ_ox w_ox + ṁ_f w_f) / (η_m Δh_T − x_ox w_ox − x_f w_f)
```

in closed form. A non-positive denominator means each kilogram through the
turbine yields less than it costs to pump, and no flow closes the cycle. The
exhaust leaves through its own duct, and its contribution is bracketed:
nothing (dumped), or an ideal frozen expansion to the main nozzle's exit
pressure. That expansion stops where the frozen isentrope reaches the species
data's lower limit, before water would condense. The delivered specific
impulse is main chamber plus exhaust over main plus generator flow, given at
both ends of the bracket.

**Staged combustion (closed).** All of one propellant passes through a
preburner, or a set fraction of it: the fuel for a fuel-rich preburner, the
oxidiser for an oxidiser-rich one. It burns with a little of the other,
drives the turbine and is injected into the main chamber. The unknown is the
turbine pressure ratio `PR`:

```
p_pb   = PR · p_c (1 + s_hg)                     s_hg: hot-gas injector stiffness
p_feed = p_pb (1 + s_pb) + f_line p_c            the preburner's supply
```

The major propellant's pump delivers `p_feed`, through the jacket first on a
fuel-rich cycle. The minor propellant reaches it through a boost stage from
its main pump's discharge, as on the RS-25. A larger `PR` gives the turbine
more work per kg but asks more of the pump feeding it, so the balance
`g(PR) = η_m P_T − P_pumps` rises from `−P_pumps` at `PR = 1`, peaks and falls.
Ignis takes its smallest root, found on a 48-point grid in `log PR` and closed
by bisection. `PR` is searched up to 10, or to the point where a pump reaches
its property table's ceiling. The preburner gas depends weakly on `PR`
(through the pump work it carries and its pressure), so gas and ratio
alternate until the mixture ratio settles to 1e-10. When `g` never reaches
zero the cycle does not close. The result reports the best power ratio
`η_m P_T / P_pumps` any admissible `PR` reaches, and is marked infeasible.

**Expander (closed).** The jacket's heated fuel drives the turbine and then
goes to the injector. The turbine expands from the jacket outlet to the
injector supply, `p_c (1 + s_inj + f_line)`. The unknown is the fuel pump's
discharge `p_d`, with `PR = (p_d − Δp_jacket) / p_c (1 + s_inj + f_line)`,
solved as above. An optional bypass sends part of the jacket flow round the
turbine. A film drawn from the jacket outlet never reaches it.

**What the closures leave out.** Pump and turbine maps and speeds, NPSH and
cavitation, bearing-coolant, igniter and tank-pressurisation bleeds, heat
exchangers, the start transient, and turbine-inlet temperature limits beyond
the one requested. Efficiencies are inputs, not predictions. The low-pressure
boost turbopumps of an engine like the RS-25 are represented only by the
pump inlet pressure.

---

## 15. Transient chamber

A single well-stirred control volume of fixed volume `V`, tracking the
oxidiser-derived and fuel-derived masses separately so the instantaneous
chamber mixture ratio is exact:

```
d m_ox / dt = ṁ_ox,in − ṁ_out (m_ox/m)
d m_f  / dt = ṁ_f,in  − ṁ_out (m_f/m)
d U    / dt = ṁ_ox,in h_ox + ṁ_f,in h_f − ṁ_out h_c
              − (1 − η_heat) Q̇_releasable + Q̇_wall
```

Because enthalpies are absolute, combustion needs no source term. The ignition
ramp is modelled by *withholding* part of the chemical energy: `η_heat(t)` is
the fraction of the heat of combustion released, and `Q̇_releasable` counts only
the part of the incoming stream that has a reaction partner (pure oxidiser or
pure fuel releases nothing). With `η_heat = 0` the injected stream is
energetically equivalent to injecting combustion products at 298.15 K, which is
a well-defined and bounded no-heat-release limit.

`η_heat` is **not** the steady combustion efficiency `η_c*`. Withholding
`(1 − η_heat) q_comb` lowers the chamber temperature by roughly
`ΔT = (1 − η_heat) q_comb / cp`, and since `c* ∝ sqrt(T)`,

```
η_c* ≈ sqrt( 1 − (1 − η_heat) q_comb / (cp T_ad) )
```

For LOX/methane `q_comb/(cp T_ad) ≈ 0.36`, so `η_c* = 0.96` corresponds to
`η_heat ≈ 0.874`, not 0.92. The shipped startup case uses the calibrated value
and reproduces the chamber pressure the steady model implies for the same
commanded flows to **0.031 %** — see [verification.md](verification.md) 1.3a
for the comparison and what the residual is made of.

**Closure.** Given `(m_ox, m_f, U)`, the temperature and pressure follow from
`u = U/m`, `ρ = m/V` and the instantaneous `O/F` by inverting a tabulated
equilibrium equation of state. The table is tensor-product cubic Hermite in
(mixture ratio, temperature, ln pressure); linear interpolation in any
direction would leave the right-hand side only C0 there, which by itself caps
the observed Runge-Kutta order at two. Interpolation error is *measured*, not
assumed: see [verification.md](verification.md).

**Outflow.** Choked flow uses the tabulated characteristic velocity evaluated at
the *actual* chamber temperature,

```
c*(O/F, T, p) = κ(O/F, p) sqrt(R T/γ_s) [(γ_s+1)/2]^((γ_s+1)/(2(γ_s−1)))
```

where `κ` is the ratio of the exact variable-property `c*` to the
constant-γ result, tabulated at the adiabatic flame condition so that it is
exactly one there. Sub-critical flow uses the compressible orifice relation
divided by the same `κ`, so the two branches agree exactly at the critical
pressure ratio. Without that, mass flow steps by a few percent as the chamber
chokes, and a single discontinuity in the right-hand side drops the observed
order of any Runge-Kutta scheme to one.

**Integrators.** Classical fixed-step RK4 and an adaptive embedded Cash-Karp
RK4(5) with PI step control. Conservation integrals are accumulated with the
same Runge-Kutta weights as the state, so comparing them against the change in
`(m, U)` checks that the conservative form survived the implementation.

**Validity.** The chamber contents are represented as equilibrium products at
all times, so the model is valid from ignition onwards, not through a cold
pre-ignition fill. A run that drives the chamber outside the tabulated mixture
ratio, or below the energy a product mixture can hold at the table's floor,
stops with a message naming the cause and the remedy.

---

## 16. Optimisation

```
maximise (or minimise)  f(x)          a named result metric
subject to              g_i(x) ≤ 0    named metrics with bounds
                        x_lo ≤ x ≤ x_hi
```

The objective is an engine analysis: expensive, not differentiable in closed
form, and able to fail outright for infeasible geometry. Ignis therefore uses

* an **augmented Lagrangian** outer loop (Hestenes-Powell-Rockafellar),

  ```
  L(x; λ, μ) = −s f(x) + (1/2μ) Σ_i [ max(0, λ_i + μ g_i(x))² − λ_i² ]
  λ_i ← max(0, λ_i + μ g_i(x))
  ```

  with `μ` increased when the worst violation fails to fall by a factor of four.
  This drives iterates *to* the constraint boundary instead of stopping short of
  it the way a fixed quadratic penalty does;
* a **Nelder-Mead** inner loop with the dimension-adaptive coefficients of
  F. Gao and L. Han (*Comput. Optim. Appl.* 51, 259–277, 2012), operating on
  variables scaled to the unit box and reflected at the bounds;
* **Latin-hypercube multi-start**, so a single local basin is not mistaken for
  the global optimum.

A failed analysis returns a large finite penalty and is counted, so the simplex
walks around infeasible regions rather than crashing. **No claim of global
optimality is made.**

---

## 17. Uncertainty and sensitivity

**Determinism.** Each sample draws from its own generator seeded by hashing the
campaign seed with the sample index through SplitMix64:

```
state_k = splitmix64( seed XOR splitmix64(k + φ) )
```

Sample `k` therefore sees the same random numbers regardless of thread count or
scheduling order. `ignis_mc --threads N` reproduces `--threads 1` byte for byte,
and both the test suite and `scripts/run_all.sh` check exactly that.

**Distributions.** Normal, lognormal, uniform and triangular, each with optional
hard bounds. A draw outside the bounds is *re-drawn*, so the realised
distribution is a genuine truncated distribution rather than a clipped one with
a spike at the bound.

**Sensitivity.** Three independent measures, because no single one is
trustworthy alone:

* **standardised regression coefficients** from a least-squares fit of the
  outputs on the inputs, `SRC_i = β_i σ_{x_i} / σ_y`, reported together with the
  model `R²` so the linearity assumption is visible;
* **Spearman rank correlation**, which survives monotone nonlinearity;
* **local elasticities** `(∂y/∂x)(x/y)` by central differences at the nominal
  point. Where a step leaves the feasible region the elasticity is reported as
  unavailable rather than invented.
