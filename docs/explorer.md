# Ignis Engine Explorer

A desktop front end for the Ignis solver. Change a chamber pressure, a mixture
ratio, an expansion ratio or an altitude and see what it does to thrust,
specific impulse, the wall temperature and the flow field — with the
constraints and the model's own limitations on screen next to the answer.

![The Explorer comparing a sea-level and a vacuum nozzle](../results/figures/16_explorer.png)

```bash
pip install -r python/requirements-explorer.txt
./scripts/build.sh          # the Explorer needs the Ignis binaries
python3 explorer.py
```

`IGNIS_BUILD_DIR=/path/to/build python3 explorer.py` points it at an existing
build.

---

## What it is, and what it is not

**The Explorer owns no physics.** It writes a configuration file, runs the same
`ignis_engine`, `ignis_equilibrium` and `ignis_nozzle` binaries the command
line drives, and reads back their JSON and CSV. A number on screen is a number
the solver produced. When the solver refuses a design — a jacket that cannot
pass the flow, an expansion past the property range — the Explorer shows the
refusal in red and leaves the charts on the last good state. It never fills the
gap with an estimate.

That also means the Explorer inherits every limitation in
[`limitations.md`](limitations.md), which is why a condensed version of that
document is pinned open in the right-hand column rather than hidden behind a
menu.

## Layout

| Column | Contents |
|---|---|
| Left | Design A and, when comparison is on, Design B: propellant, chamber pressure, mixture ratio, throat radius, expansion ratio, altitude, η<sub>c\*</sub>, frozen vs shifting expansion, and the jacket's channel height and wall thickness. |
| Centre | Five chart tabs — contour, axial profiles, wall and coolant, chamber composition, altitude sweep. |
| Right | Eight performance cards, the constraint and warning panel, and the assumptions panel. |

**Comparison.** Ticking *Compare A / B* opens a second design slot and overlays
both on every chart; the cards gain a Δ line showing A − B in absolute and
percentage terms, coloured green or amber according to which direction is
better for that quantity. The two slots keep fixed hues everywhere — cyan for
A, amber for B — so identity never moves when a chart is rescaled.

**Presets.** Ignis-M1 at sea level, Ignis-M1 with the vacuum nozzle, and the
Ignis-H1 LOX/hydrogen upper stage. They match the shipped configurations, so
the Explorer and `scripts/run_all.sh` can be checked against each other.

## The charts

* **Contour** — axisymmetric half-section, true to scale, the wall coloured by
  the local Mach number. Under comparison the two nozzles are stacked on a
  shared colour scale so the difference in size is read directly.
* **Axial** — static pressure (log), static temperature, velocity and Mach
  number along the axis, with the sonic throat marked.
* **Thermal** — wall heat flux, and the hot-side wall against the coolant bulk
  temperature, with the material limit drawn and named.
* **Composition** — chamber equilibrium mole fractions, ranked, on a log axis.
* **Altitude** — thrust and specific impulse from sea level to 80 km, with a
  shaded band where the Summerfield criterion predicts separation. The curve
  inside that band is optimistic, and the caption says so.

The altitude curve comes from `ignis_nozzle`, not from a Python atmosphere:
there is one implementation of the U.S. Standard Atmosphere in this repository
and it is the validated C++ one.

## Constraints panel

Six live checks, each with its value and the bound it is judged against:

| Check | Green when |
|---|---|
| Separation margin | the exit pressure is above the Summerfield separation pressure |
| Expansion regime | the nozzle is ideally expanded |
| Peak wall temperature | below 90 % of the material limit |
| Jacket pressure drop | below 35 % of the coolant inlet pressure |
| Coolant boiling | no sub-critical boiling detected |
| Solver residuals | element and mass-flow residuals below 1e-10 |

Below them, the solver's own warnings verbatim.

## Colour

Chart colour is assigned by the job it does, and the categorical palette was
validated against this theme's chart surface rather than chosen by eye — worst
adjacent colour-vision separation 14.1 (OKLab ΔE×100, 8 is the target), worst
adjacent normal-vision separation 16.2, every slot inside the dark lightness
band and above 3:1 contrast against the surface. Magnitude (Mach number) uses a
single-hue sequential ramp whose darkest step still clears the surface;
constraint state uses a reserved status palette that is never reused for a
series. No chart has a second y-axis.

## Performance

One design costs three solver invocations — the engine analysis, the
equilibrium export and the altitude sweep — and lands around 0.4 s on the
machine in [`benchmarks.md`](benchmarks.md). Solving runs on a worker thread,
so the window stays responsive, and the Run button re-enables when every
visible slot has finished. Changing a control marks the design stale in the
status bar rather than re-solving on every keystroke; Run, or Ctrl+Enter,
solves.
