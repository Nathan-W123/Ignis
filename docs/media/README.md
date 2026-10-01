# Media

Screenshots, animations and graphs of Ignis made for the portfolio site
(nathan-w123/Website, `public/projects/ignis/`). They were produced at commit
`89231c3` (before `a356e14` regenerated `results/`), so numbers in them are that
commit's; `results/` and `docs/validation.md` are the reference.

| File | What it is | How it was made |
|---|---|---|
| `rs25_plume.{mp4,webm,webp}` | The RS-25 plume marched from rest, 3-D volume view inside NASA's model of the bell (151 frames, 21.1 ms of flow) | Explorer Flow tab, 300×90 grid, frames read back from the GL view under Xvfb + Mesa; the viewport's caption text is masked |
| `explorer_flow_tab.{mp4,webm,webp}` | The whole Explorer window during the same march | Window grab with the GL view composited, every frame |
| `explorer_thermal_tab.png` | Thermal tab, constraints and "Against the real engine" panel, RS-25 live solve | `Explorer.grab()` after the RS-25 preset solved with the compiled binaries |
| `m1_wall_temperature.png` | Ignis-M1 hot-wall temperature and heat flux: Bartz, integral boundary layer, boundary layer + 3 % film | `ignis_engine -c configs/methane_nominal.yaml` with `hot_gas_model` / `film.fuel_fraction` varied |
| `heat_transfer_vs_experiments.png` | JPL TR 32-415 medians and NASA TN D-2832 station ratios, boundary layer vs Bartz | Values from `docs/validation.md` §5b |
| `rs25_vs_published.png` | RS-25 engine figures (live solve) and the turbopump / preburner checks | Explorer panel values; `docs/validation.md` §5d |
| `monte_carlo_wall_before_after.png` | Peak wall temperature over 2,000 samples, v1.0.0 (Bartz) vs boundary layer + film | `results/monte_carlo` as shipped then, against `ignis_mc -c configs/monte_carlo.yaml` at `89231c3` |
