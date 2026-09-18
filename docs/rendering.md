# Rendering the flow

Ignis produces figures for reading and renders for looking at. This page covers
the renders: what they show, how they are computed, and what they are not.

The short version: **nothing in a render is painted on.** The colour at every
point is Planck's law at the temperature the solver computed, the brightness is
the density it computed, and the exhaust plume is a separate CFD solution
started from the nozzle exit state rather than an artistic impression of one.

Everything here lives in `python/ignis_viz/` and is driven by
`tools/make_cover.py`.

---

## 1. What is being integrated

The gas in a running engine is genuinely incandescent — 3500 K in the chamber,
still 2100 K at the exit plane — so the honest picture of the flow field is a
picture of its own thermal emission.

For an optically thin gray emitter the emission coefficient is

    j = kappa rho B(T),

so the radiance reaching the camera along a ray is the line integral of
`rho sigma T^4`, weighted by the blackbody chromaticity at the local
temperature. That is exactly what is integrated.

`ignis_viz.blackbody` turns a temperature into a colour the only way that is
defensible: Planck's law sampled across the visible band, integrated against
the CIE 1931 colour-matching functions (the analytic multi-lobe fits of Wyman,
Sloan and Shirley, 2013), converted to linear sRGB and desaturated to the
gamut. 5778 K comes out very slightly warm white, 9000 K blue, 1500 K deep
orange — which is the check that it is right.

The grey absorption coefficient `kappa` is unknown, so the absolute scale is
arbitrary and the image is normalised for display. The *relative* structure —
which part of the engine is brighter than which, and by how much — is the
solver's.

### The response curve

The scene covers an enormous dynamic range: between the chamber and the exit
plane the density falls by a factor of 85 while T⁴ falls by a factor of 8, so
the chamber is some seven hundred times brighter than the exhaust leaving the
bell. A linear mapping would show the chamber and nothing else. The renders
therefore use a logarithmic response, as a camera does with a scene like this,
and a filmic shoulder applied to **luminance** so that the highlight keeps its
hue instead of clipping to white one channel at a time.

---

## 2. The nozzle

The quasi-1D solution is a function of axial position alone, so revolving it
about the axis says exactly that and no more: every cross-section is uniform.

The shell is drawn as a **half cutaway** — the near half removed, the far half
kept and shaded. This is a choice, and the reason for it is worth stating: a
translucent shell dims the gas behind it, which puts a step in brightness
wherever the wall ends. At the exit plane the nozzle gas would sit three times
darker than the plume it flows into, for no physical reason at all. Cutting the
near half away leaves the gas column unattenuated from injector to plume.

Ray marching against a fixed geometry is the expensive part, so for animations
of a *changing* flow through *fixed* geometry it is done once: `build_geometry`
returns the path length each pixel spends in each axial slice, and a frame is
then a matrix multiply against the per-slice emission.

---

## 3. The plume

Ignis stops at the exit plane, which is the right place for a quasi-1D internal
flow model to stop. But the plume is the part of a rocket engine everyone
recognises, and its shock structure is the visible consequence of the exit
state the solver computed — so `ignis_viz.plume` computes one.

It solves the axisymmetric Euler equations

    d(U)/dt + d(F)/dx + d(G)/dr = -(1/r) (rho v, rho u v, rho v^2, v(E+p), rho v Y)

by a second-order finite-volume method: MUSCL reconstruction with a minmod
limiter on primitive variables, an HLLC approximate Riemann solver at the
faces, and a two-stage SSP Runge-Kutta step. The jet enters through the exit
plane at the Ignis exit state and expands into still ambient air. `Y` is a
passive scalar, 1 in jet gas and 0 in ambient air, carried so that thermal
emission is never attributed to the surrounding atmosphere.

This resolves the expansion fan or lip shock, the barrel shock, the Mach disc
where one forms, and the shock-cell train downstream.

### Validation

An over-expanded jet leaves the nozzle below ambient pressure and must be
compressed back up to it through an oblique shock attached to the lip. The
angle of that shock is not a correlation — it has an exact answer. Given the
upstream Mach number and the required pressure ratio, the normal-shock relation

    p2/p1 = 1 + 2 gamma / (gamma + 1) (Mn1^2 - 1)

fixes the normal Mach number, and the wave angle follows from
`Mn1 = M1 sin(beta)`.

For the Ignis-M1 at 3 km (`M1 = 3.4595`, `gamma = 1.178`, `pa/pe = 2.039`):

| | wave angle |
|---|---|
| oblique-shock theory | 23.88° |
| computed, 600 x 200 grid | 25.19° |
| difference | 1.31° |

A coarser 320 x 110 solve puts it at 24.77°, 0.89° from theory. The finer grid
is *further* from the exact answer, not closer, because it resolves more of the
shear-layer rollup that the straight-line fit to the shock has to average
through — the shock position itself is sharper, the band it is measured over is
noisier. Both are within about a degree and a half, which is the honest claim.

`tools/make_cover.py` prints this comparison on every run.

### What the plume model does not do

This is an inviscid model and it is **not part of the validated Ignis core**:

* **No viscosity and no turbulence model.** The jet shear layer spreads only by
  numerical diffusion, so the shock cells persist much further downstream than
  they would in a real plume, and the jet edge is far sharper.
* **One ratio of specific heats for the whole domain**, frozen at the nozzle
  exit value. The ambient air really has gamma = 1.4. The jet core structure is
  set by the jet's own gamma and by pe/pa, and the ambient contributes nothing
  to the emission, so a single gamma is used rather than a second species
  equation being carried.
* **Frozen chemistry.** The products leave the nozzle fuel-rich and a real
  plume at low altitude afterburns with entrained air, which is bright. None of
  that is modelled, so the computed plume is dimmer and cooler downstream than
  a real one would be.
* **No radiation transport.** Emission is treated as optically thin, as in the
  nozzle.
* The residual plateaus around 0.3 rather than converging to zero, because the
  inviscid shear layer is unstable and rolls up continuously. The shock cells
  themselves are steady; the solution is a snapshot of a statistically steady
  flow, not a fixed point.

### Choosing an operating point

The nominal M1 is over-expanded at sea level and Ignis predicts the flow
**separates inside the bell** there (separation margin −0.061, separation at an
area ratio of 17.6 out of 20). Rendering an attached sea-level plume would
contradict the solver, so `make_cover.py` refuses to do it and the default
operating point is 3 km, the lowest round altitude at which the altitude sweep
reports the nozzle running full.

---

## 4. Making the flow visible

A steady solution has nothing moving in it, and a render of one is a still life
however bright it is. `ignis_viz.tracers` carries tracer particles through the
solved field so that the flow reads as flow.

Each tracer moves at the local velocity the solver computed — the quasi-1D
axial velocity inside the nozzle, where it keeps its fractional radius `r/R(x)`
(the streamline of a quasi-1D flow), and the two-dimensional velocity of the
Euler plume outside it. The streak a tracer leaves in one frame is therefore
literally speed times exposure, as it would be on a camera: the chamber barely
smears while the exhaust draws long lines.

**Where the tracers are is a visualisation choice and is not physics.** Injected
at a steady rate from the injector face they would settle at a number density
proportional to the gas density, which falls by a factor of eighty-five through
the engine, leaving the exhaust almost empty. Instead each tracer is given a
finite life and respawns somewhere in the flow, which spreads them evenly over
the picture. Their density means nothing. Their direction, their length and
their colour are the solution.

---

## 5. Running it

```bash
pip install -r python/requirements-render.txt

# everything: hero still, flow animation, ignition animation, camera move
python3 tools/make_cover.py

# recompute the plume rather than using the cache (tens of minutes)
python3 tools/make_cover.py --solve

# just one of them, at a different altitude
python3 tools/make_cover.py --only hero --altitude 6000
```
Outputs land in `results/figures/`:

| file | what |
|---|---|
| `17_engine_render.png` | the still |
| `19_engine_flow.mp4` | the engine running, flow visible, 3000x slow motion |
| `20_engine_ignition.mp4` | the engine lighting up, from the start-up transient |
| `21_engine_move.mp4` | a start with the camera moving through it |

### Time in the animations

The ignition sequences are driven by the zero-dimensional start transient,
which gives chamber pressure and temperature against time; the nozzle field is
scaled onto that state quasi-steadily — the same approximation the transient's
own thrust estimate makes — before its emission is integrated.

Exposure is set once, from the engine at full thrust, and held. Keyed to the
opening frame instead it would be exposing for a chamber that has not lit yet,
and everything after ignition would be white; more to the point, a frame that
looks brighter should mean an engine that got brighter, not a curve that moved.

**The plume appears part-way through a start, abruptly.** That is the model
being honest rather than the render being lazy. During the ramp the chamber
climbs from 0.1 to 5.5 MPa, so the exit-to-ambient pressure ratio sweeps from
about 0.01 up to 0.49, and for most of that the nozzle is so far over-expanded
that Ignis's separation criterion says the flow separates inside the bell.
There is no attached plume until chamber pressure passes about 3.9 MPa, so none
is drawn until it does. Fading one in would look smoother and mean less.

### Camera moves

The scene is axisymmetric, so swinging the camera *around* the nozzle axis
gives an identical picture frame after frame — the one direction of motion that
buys nothing. What changes the view is the angle between the line of sight and
the axis, so the move swings downstream from a square side elevation to an
oblique three-quarter, rising as it goes and pulling back from the throat to
the whole plume.

A moving camera means new geometry every frame, so none of the march can be
reused; a static camera over a steady field marches once and every frame after
that is a matrix multiply. That is the whole difference in cost between
`19_engine_flow.mp4` and `21_engine_move.mp4`.
