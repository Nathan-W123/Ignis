# RS-25 — the Space Shuttle Main Engine, Block IIA

The Explorer opens on this engine: NASA's 3-D model of the RS-25 nozzle bell,
drawn around an Ignis solution of the engine's published geometry and
operating point. This file says where every piece of that comes from, what
was inferred rather than read, and how Ignis's answer compares with the
manufacturer's.

## What is in this directory

| File | What it is |
|---|---|
| `nasa_space_shuttle_d_eng.glb` | The engine bell from NASA's *Space Shuttle (D)* model, **byte for byte** — only the file name is changed. SHA-256 `e31cfb57e3365f036d40cee79ac9e2f8b5f4c48849f6aa9e19aee4075846a0bc`. |
| `engine.json` | Units, the axis fitted to the mesh, the preset the model belongs to, and the published figures the Explorer compares against. |
| `NOSA-1.3.txt` | The NASA Open Source Agreement, see below. |

**Source.** [nasa/NASA-3D-Resources](https://github.com/nasa/NASA-3D-Resources),
`3D Models/Space Shuttle (D)/Space Shuttle (D) eng.glb`, at commit
`11ebb4ee043715aefbba6aeec8a61746fad67fa7`.

**Licence.** The repository's README says its assets "are free and without
copyright" and points to NASA's media usage guidelines; its `meta.json` lists
the *NASA Open Source Agreement Version 1.3*. The model is redistributed here
unmodified, with that agreement alongside, which satisfies either reading. No
NASA insignia is used, and nothing here implies NASA reviewed or endorses
Ignis.

## What the model is, and what it is not

It is the bell as it appears on the orbiter, modelled in **inches at real
scale**. It has two surfaces:

* `eng out` — the outside of the bell, with its nine hat bands and the three
  feed ducts and steerhorns. **This is what is drawn.**
* `eng in` — the inside, modelled as a coarse straight cone. It is not drawn
  (the solved contour is the inside), and it is used for one thing: fitting
  the axis.

It has no chamber, no throat and no powerhead — on the Shuttle those are
inside the orbiter's aft fuselage. The Explorer draws those parts from the
solved contour, in the model's colour, and says so in the viewport caption.

**It is display only.** No number Ignis computes reads the mesh. The model is
drawn only when the solved geometry is this engine's geometry (throat,
expansion ratio, contraction ratio, chamber length and bell length), so
editing any of those removes it rather than leaving it wrapped around a
nozzle it no longer describes.

### Placing it

`tools/fit_engine_model.py` fits circles to the two end rings of `eng in` and
runs the axis through their centres. The exit ring carries both edges of the
nozzle lip — 44.97 in and 47.5 in from the axis, against the 45.15 in and
47.0 in that Rocketdyne's 90.3 / 94 in inside/outside exit diameters give —
so the tool fits the two edges separately and uses the gas side, where the
residual is 0.0001 in. The exit plane is tilted 0.85° to the axis. The model
is scaled by 0.0254 m/in and set with its exit on the solved exit plane.

As an independent check of both: the model's outer skin lies **0.2 to 4 in
outside** Ignis's solved gas-side wall along the length of the bell (median
2–4 in — tube wall, jacket and hat bands). 156 of its 16,464 vertices, all
on the exit lip, sit up to 0.24 in inside it. The model's bell is 102.1 in
long, from the exit to where it meets the orbiter's heat shield; Rocketdyne
gives 121 in from throat to exit.

## The design Ignis solves

Every input is from Rocketdyne's *SSME Orientation* training manual [1],
Block IIA at 104.5 % of rated power level, unless the row says otherwise. The
preset is `RS25` in `python/ignis_explorer/solver.py`.

| Input | Value used | As published | Where |
|---|---|---|---|
| Propellants | LOX / LH2 | oxygen/hydrogen | [1] engine summary |
| Mixture ratio | 6.032 | 6.032 lb O₂ per lb H₂ | [1] engine highlights |
| Chamber pressure | 19.753 MPa | 2,865 psia, *throat stagnation pressure* | [1] MCC operating parameters |
| Throat radius | 138.21 mm | throat area 93.02 in² | [1] MCC geometry |
| Contraction ratio | 2.66 | 2.66:1 | [1] MCC geometry |
| Expansion ratio | 69 | 69:1 | [1] nozzle geometry |
| η<sub>c\*</sub> | 0.996 | "two-stage combustion approximately 99.6 % efficient" | [1] engine highlights |
| Coolant slots | 430 | 430 | [1] MCC geometry |
| Coolant flow | 18.71 % of the hydrogen | 29 lb/s of 155 lb/s | [1] MCC parameters; flow schematic |
| Coolant inlet | 52.04 K, 38.93 MPa | −366 °F, 5,647 psia | [1] MCC operating parameters |
| Jacket end | area ratio 4.48 | MCC expansion ratio to the nozzle attach flange 4.48:1 | [1] MCC geometry |

Ignis takes chamber pressure as a stagnation pressure (an infinite-area
combustor, `docs/theory.md` §5), so the manual's *throat stagnation pressure*
is exactly the quantity it needs.

**Inferred, not read:**

| Input | Value used | How |
|---|---|---|
| Cylindrical chamber length | 136.59 mm | the value that puts Ignis's throat **14.7 in** from the injector face, the published injector-end-to-throat length [1] |
| Bell length | 81.4 % of a 15° cone | the value that makes Ignis's nozzle **121 in** from throat to exit, the published length [1] (Ignis gives 121.005 in) |
| Channel depth | 4.32 mm (0.170 in) | not given for this chamber. The mean of the throat and aft depths NASA published for the earlier 390-slot liner, 0.093 and 0.247 in [2]. Ignis's jacket has one depth; the real slots taper |
| Hot-wall thickness | 0.7 mm | not published in these sources; the value the Ignis-H1 uses |
| Liner material | CuCrZr | the real liner is NARloy-Z, a copper–silver–zirconium alloy [1]; Ignis's material library has no NARloy-Z |
| Bell angles | 33° / 8° | Ignis's defaults; not published in these sources |
| Altitude | 6 km | see *Sea level* below |

## Ignis against Rocketdyne's figures

The Explorer shows this table in its *Against the real engine* panel whenever
the design is this operating point. It is a **sanity check, not a
validation**: one operating point, rounded figures from a training manual,
no uncertainties. Nothing in the design was adjusted to improve it.

| Quantity | Ignis | Published [1] | Ignis − published |
|---|---|---|---|
| Exit diameter | 2.296 m | 2.294 m (90.3 in inside) | +0.1 % |
| Vacuum specific impulse | 460.5 s | ≈ 452 s | +1.9 % |
| Propellant flow | 512.8 kg/s | 494.0 kg/s (155 + 934 lb/s) | +3.8 % |
| Vacuum thrust | 2,316 kN | 2,188 kN (491,900 lbf) | +5.8 % |
| Chamber coolant exit temperature | 419 K | 265 K (17 °F) | +58 % |
| Chamber coolant pressure drop | 5.74 MPa | 8.32 MPa (5,647 − 4,441 psia) | −31 % |
| Peak hot-gas wall temperature | 1,052 K | 811 K (1,000 °F) | +30 % |

What accounts for each gap, as far as these sources allow:

* **Specific impulse, +1.9 %.** Ignis expands an ideal gas mixture in shifting
  equilibrium with no boundary layer and no finite-rate chemistry. Both of
  those cost a real engine performance, so Ignis should read high; this is
  the direction expected.
* **Propellant flow, +3.8 %.** Ignis computes flow as p<sub>c</sub>A<sub>t</sub>/c\*. The
  manual's pressure, throat area and flow together imply c\* = 2,400 m/s,
  3.4 % *above* the ideal equilibrium c\* Ignis computes for this mixture
  (2,321 m/s), and no combustion efficiency exceeds one. So those three
  numbers cannot all be what an ideal-flow model takes them to be — the
  effective throat may be smaller than its geometric area, or the pressure
  defined differently. A discharge coefficient of 0.963 would reconcile
  them; that is an inference, not something the sources state, and Ignis
  does not apply one.
* **Vacuum thrust, +5.8 %,** is the two above together (1.038 × 1.019).
* **Coolant exit temperature, +58 %.** That is heat: from the coolant's
  enthalpy (`data/coolants/hydrogen.csv`), the published temperatures and flow
  imply 38 MW into the chamber's coolant; Ignis puts in 71 MW, 1.85 times as
  much. Ignis's hot-gas side is the Bartz correlation, which over-predicts in
  Ignis's own validation against JPL nozzle data (`docs/validation.md`), and
  SSME chambers are fuel-film cooled near the injector [3], which Ignis does
  not model.
* **Peak wall temperature, +241 K,** follows from that heat, and from the
  single channel depth: the real slots are shallowest at the throat (0.093 in
  [2] against the 0.170 in Ignis uses everywhere), so the real coolant runs
  fastest exactly where the heat flux peaks. The Explorer's constraint panel
  flags Ignis's 1,052 K against CuCrZr's 800 K limit; that is Ignis's
  prediction, and the table above is why it should not be read as the engine.
* **Coolant pressure drop, −31 %,** is the same simplification the other way:
  the real throat slots, at about half the depth Ignis uses, dominate the
  real pressure drop.

## Sea level

At this operating point the exit pressure is 20.6 kPa, and Ignis's
separation criterion — Summerfield's, separation below 0.4 of ambient —
predicts the nozzle separates at sea level. The real engine runs at this
power level at sea level, on the launch pad and on its test stands, with
its nozzle flowing full. Summerfield's criterion is an empirical fit for
*free-shock* separation; thrust-optimised bells like the RS-25's separate in
a different, restricted-shock pattern, which Ignis does not model.

The preset is evaluated at **6 km**, the lowest whole kilometre at which the
criterion agrees the nozzle flows full (exit-to-ambient pressure 0.435; at
5 km it is 0.38, and the Flow tab would refuse to march an attached plume).
There the jet is still over-expanded, so the plume forms the shock diamond
and Mach disc of the familiar low-altitude picture. None of the comparison
above depends on the altitude: it is vacuum performance and the chamber's own
cooling.

## References

1. *Space Shuttle Main Engine Orientation*, Space Transportation System
   Training Data, Boeing Rocketdyne Propulsion & Power, June 1998.
   <http://large.stanford.edu/courses/2011/ph240/nguyen1/docs/SSME_PRESENTATION.pdf>
   — engine summary and highlights; Block IIA propellant flow schematic at
   104.5 % RPL; main combustion chamber geometry and Block IIA operating
   parameters at 104.5 % RPL; nozzle geometry.
2. *SSME Main Combustion Chamber (MCC) "Hot Oil" Dewaxing*, NASA Technical
   Reports Server 19950025356 (N95-31777): "390 milled axial coolant channels
   … ranging from 0.035" x 0.093" at the throat to 0.060" x 0.247" at the aft
   end of the liner."
   <https://ntrs.nasa.gov/citations/19950025356>
3. T.-S. Wang and V. Luong, "Hot-Gas-Side and Coolant-Side Heat Transfer in
   Liquid Rocket Engine Combustors," *Journal of Thermophysics and Heat
   Transfer* 8(3), 1994, pp. 524–530 — Table 3: fuel film-coolant flow
   4.84 lb/s in the standard-throat SSME chamber and 1.07 lb/s in the
   large-throat (LT VPS AMCC) design the paper analyses; the same paper's
   standard-throat case puts the coolant at 6,461 psi in and 5,137 psi out.
