# REGOLIT: REworking and Gardening of Lunar Impacted Terrains
## A surface evolution model for impacted surfaces to study mechanical gardening of ice and regolith
The model is based on [Richardson et al. 2009 CTEM model](https://www.sciencedirect.com/science/article/pii/S0019103509003194?casa_token=j6uDz1cdmAkAAAAA:UsgleZ2OBuARNfT8Gj0a2jaye59Fh9o4tzBj2rApSYEn_61GKxn3XCTfej-JPxHY2O2Un595JA)
and adds an efficient 3-D description of the subsurface using layers.

### Installation
1. Requirements: a C++17 compiler (g++ or clang++) and GNU make. There are no external dependencies.
2. To compile, simply run `make` in the `REGOLIT` directory (`make debug` builds without optimization and with debug symbols).
3. To run, execute `./build/apps/regolit_main.run` from the `REGOLIT` directory (it reads `config/` and writes `output/` and `log/`).
4. Secondary craters (`isEmplaceSecondaries 1`) follow `N(>r) = (r_max / r)^slope_secondaries` between
   one pixel and `secondaryLargestFraction` of the primary radius, with their own `secondaryDepthToDiameter`.
5. The output will be saved in the `output` directory: binary files containing the surface elevation, surface composition and subsurface layers, and `existing_craters.txt`, the list of craters still visible at the end of the run.

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
