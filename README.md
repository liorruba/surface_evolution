# REGOLIT: REworking and Gardening of Lunar Impacted Terrains
## A surface evolution model for impacted surfaces to study mechanical gardening of ice and regolith
The model is based on [Richardson et al. 2009 CTEM model](https://www.sciencedirect.com/science/article/pii/S0019103509003194?casa_token=j6uDz1cdmAkAAAAA:UsgleZ2OBuARNfT8Gj0a2jaye59Fh9o4tzBj2rApSYEn_61GKxn3XCTfej-JPxHY2O2Un595JA)
and adds an efficient 3-D description of the subsurface using layers.

### Installation
1. Requirements: a C++17 compiler (g++ or clang++) and GNU make. There are no external dependencies.
2. To compile, simply run `make` in the `REGOLIT` directory (`make debug` builds without optimization and with debug symbols).
3. To run, execute `./build/apps/regolit_main.run` from the `REGOLIT` directory (it reads `config/` and writes `output/` and `log/`).
4. The crater scaling follows Holsapple (1993): `k1`, `mu`, the effective strength `Ybar` and the
   strength constant `k2` (optional, default 1 as in the original model; 0.26 for soils and regolith).
5. Secondary craters (`isEmplaceSecondaries 1`) are ejecta fragments, see "Secondary craters" below.
6. The output will be saved in the `output` directory: binary files containing the surface elevation, surface composition and subsurface layers, and `existing_craters.txt`, the list of craters still visible at the end of the run.

### Secondary craters

Secondaries are formed from the ejecta of every primary rather than scattered at random. The Z-model
already divides the transient cavity into launch annuli with an ejection speed and a landing ring
(Richardson 2009); the secondary model (`include/secondaries.hpp`, `src/secondaries.cpp`) uses them:

- Ejecta landing slower than `secondaryMinimumVelocity` (default 20 m/s) only builds the continuous
  blanket, which the grid already emplaces. Faster annuli are treated as fragment populations.
- Fragments follow a cumulative size-frequency distribution `N(>L) = (L_max / L)^slope_secondaries`
  per annulus, weighted by the annulus volume. The largest fragment is anchored where the secondary
  field begins, in the annulus landing about three primary radii out: there it makes a crater of
  `secondaryLargestFraction` times the primary diameter (about 0.05 on the Moon, Allen 1979; Melosh
  1989). Faster annuli have smaller largest fragments, `L_max(v) = L_anchor (v / v_anchor)^-secondaryVelocityExponent`
  (1 for spallation scaling, Melosh 1984; Vickery 1986, 1987); slower ones are capped at `L_anchor`.
  The fragments of an annulus never carry more volume than the annulus.
- Each fragment lands in its annulus' landing ring at a random azimuth and forms a crater with the
  same pi-scaling as the primaries, using the fragment's mass (target density) and landing speed,
  with `secondaryDepthToDiameter`. A fragment forms a distinct crater only where that crater is deeper
  than the primary's blanket at the landing distance; nearer the rim it is buried, which places the
  inner edge of the secondary field beyond the continuous ejecta without a size-dependent rule.
  Only craters at least two pixels across are formed, and only the largest
  `maximumSecondariesPerPrimary` of a primary (the size floor is raised so the population fits).
- Fragments that land outside the domain are dropped. With `isEmplaceDistantSecondaries 1` the
  primaries that form outside the domain, out to `secondaryMaximumRange` beyond its edge, are sampled
  too (in square rings of doubling width, each restricted to the primaries large enough to deliver a
  resolvable secondary at that distance, at the same flux as the domain), and the fragments they send
  into the domain are formed. This supplies the background of distant secondaries that dominates the
  small-crater population near large primaries.

A test mode forms a single prescribed crater instead of the random population: `testCraterDiameter`
(final diameter, m; 0 = off) at (`testCraterX`, `testCraterY`), with its ejecta, ghosts and
secondaries; the surface is printed before and after it (two output steps).

Ballistics are flat-surface at the Z-model launch angles; fragments' landing speed equals their
ejection speed (no atmosphere), and the crater scaling uses the full landing speed.

### Ice and soot

`isVolatiles 1` turns on the volatile processes: periodic ice deposition (`iceEmplacementInterval`,
`iceEmplacementThickness`), sublimation (`sublimationInterval`, `sublimationThickness`; the sublimation
model itself is still a stub, `Grid::sublimateIce`, so these are reserved) and the loss of
ice and soot from the ejecta (`ejectaVolatileRetention`, `ejectaSootRetention`). With `isVolatiles 0`
(the default) none of these act, and ice and soot in the initial layers are passive tracers that the
craters redistribute.

### Seismic shaking and slope collapse

Following Richardson et al. (2005) and Richardson (2009, Section 2.6), every impact shakes the
surface around it (`include/seismic.hpp`, `src/seismic.cpp`):

- The impact radiates `seismicEfficiency` of its kinetic energy as seismic energy that spreads in a
  thin hemispherical shell and is attenuated by scattering, E(l) = (η/12) π ρ_i v² d³ exp(−2πf l² /
  (K_s π² Q)) with K_s = v_s l_s / 3 (`prim_seis_freq`, `Q_factor`, `seis_wave_vel`, `seis_mean_free`).
  The peak acceleration a = 2πf √(2ε/ρ_t) defines the seismic range where it exceeds
  `seismicAccelerationThreshold` times g.
- Inside that range the regolith receives one dose of downslope diffusion, Δz = ∇·(K ∇z) with
  K = `Cs` v^`Ki_a` D^`Ki_b` g^`Ki_c` / l^`Ki_d` (Eq. 32; D is the impactor diameter, as fitted in
  Richardson et al. 2005, Fig. 14B; l is floored at the crater radius), applied in flux form within the
  explicit stability limit, followed by the collapse of slopes above `angleOfRepose`. The settling is
  local: the wrapped box of the seismic range, or at least 1.5 crater radii for the wall collapse. There
  is no global relaxation at output steps any more (one final check remains at the end of the run).
- Distant primaries (secondary-crater model) shake the domain from outside when it lies within
  their range. Ghost craters do not shake again: the primary's settling wraps around the domain.
- `isSeismicShaking 0` keeps only the local wall collapse after each crater.
- Cost control: the relaxation sweeps only the neighbourhood of the cells that moved in the previous
  iteration; elevation changes smaller than `minimumLayerThickness` are accumulated per cell and
  applied to the columns once they add up (and flushed before every output), so a whole-domain shake
  does not rewrite every column for micrometres of movement. The log reports the time spent in the
  settling phases.

Defaults (η = 10⁻⁵, Cs = 10⁻², f = 15 Hz, Q = 1500, l_s = 1 km, v_s = 3 km/s, threshold 1 g) sit in
the middle of Richardson's ranges; the overall degradation rate is uncertain by an order of magnitude.

### Debugging
`make debug` builds `build/apps/regolit_main_debug.run` with AddressSanitizer and
UndefinedBehaviorSanitizer enabled. Run it exactly like the release binary; it is a few times slower
and aborts with a stack trace on memory errors and undefined behavior.

### Python tools
The `python/regolit` package reads every output file and can drive runs. Install it in editable
mode with `pip install -e python/` (needs numpy and matplotlib), or set `PYTHONPATH=python`.

```python
import regolit

# Run with the repository config, overriding a few parameters. The binary is built if needed.
out = regolit.run({"regionWidth": 500, "endTime": 50, "randomSeed": 7}, workdir="runs/test")

out.x, out.y                    # cell-center coordinates of the (downsampled) output grid, m
out.times                       # time of every printed step, Ma
z = out.elevation()             # final elevation, indexed [y, x]
soot = out.surface_fraction("soot", step=3)        # surface composition at a step
soot_20cm = out.integrated_fraction("soot")        # composition integrated to depthToIntegrate
bins, counts = out.histogram("craters")            # craters, impactors, depth, existing_craters
craters = out.craters()                            # visible craters: x, y, diameter, depth, initial_depth
sub = out.subsurface()                             # full-resolution layer stacks (isPrintSubsurface 1 or 2)
ice_at_1m = sub.composition_at_depth(1.0)[..., 1]  # regolith, ice, soot fractions 1 m below the surface

regolit.quicklook(out, save="runs/test/quicklook.png")
```

Existing output can be read without running: `out = regolit.RegolitOutput("output")`.
From the shell: `python -m regolit.driver --set endTime=50 --workdir runs/test --quicklook` and
`python -m regolit.plot output --save quicklook.png`. The MATLAB readers in `vis/` remain available.

### Web UI
`web/` holds a small FastAPI application that runs the model with parameters chosen in the browser
and shows shaded-relief and composition maps, layered subsurface cross-sections along a line the
user places on the map, histograms and an animation of every run. `web/README.md` explains how to
run it locally and how to deploy it (on the compute machine behind a tunnel, or on the droplet).

### An example topography evolution simulation:
![Surface evolution](https://github.com/liorruba/surface_evolution/blob/master/craters.gif)
