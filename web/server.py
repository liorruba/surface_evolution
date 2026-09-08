"""REGOLIT web UI.

A small FastAPI service: the browser posts a parameter set, the server runs the C++ model in a
run directory, and the page shows maps (shaded relief, elevation, composition), layered subsurface
cross-sections, histograms and animations rendered on the server.
Run locally with ``uvicorn web.server:app --reload`` from the repository root, or see deploy/.

Environment variables
---------------------
REGOLIT_BINARY            path of the model executable (default build/apps/regolit_main.run)
REGOLIT_WEB_RUNS          directory holding the runs (default runs/web)
REGOLIT_WEB_USER/PASSWORD if set, HTTP basic authentication is required for every page
REGOLIT_WEB_CONCURRENCY   simultaneous model runs (default 2); REGOLIT_WEB_QUEUE waiting runs (default 8)
REGOLIT_WEB_TIMEOUT       seconds allowed per run (default 120)
REGOLIT_WEB_MAX_RUNS      runs kept on disk, oldest deleted first (default 500)
REGOLIT_WEB_MAX_DISK_GB   disk budget of the runs directory, oldest deleted first (default 200)
OMP_NUM_THREADS           threads per model run (default: physical cores / concurrency)
"""
from __future__ import annotations

import asyncio
import json
import math
import os
import re
import secrets
import shutil
import sys
import threading
import time
import zipfile
from collections import OrderedDict
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LightSource  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from fastapi import Depends, FastAPI, HTTPException, Query  # noqa: E402
from fastapi.responses import FileResponse, HTMLResponse, PlainTextResponse, Response  # noqa: E402
from fastapi.security import HTTPBasic, HTTPBasicCredentials  # noqa: E402
from fastapi.staticfiles import StaticFiles  # noqa: E402
from pydantic import BaseModel, Field  # noqa: E402

import regolit  # noqa: E402
from regolit import scaling  # noqa: E402
from regolit.io import RegolitOutput, Subsurface, read_config, read_layers  # noqa: E402

# ----------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------
STATIC_DIR = Path(__file__).resolve().parent / "static"
BINARY = Path(os.environ.get("REGOLIT_BINARY", REPO_ROOT / "build" / "apps" / "regolit_main.run"))
RUNS_DIR = Path(os.environ.get("REGOLIT_WEB_RUNS", REPO_ROOT / "runs" / "web"))
SETTINGS_DIR = Path(os.environ.get("REGOLIT_WEB_SETTINGS", RUNS_DIR.parent / "settings"))  # saved run setups
MAX_SAVED_SETTINGS = int(os.environ.get("REGOLIT_WEB_MAX_SETTINGS", "200"))
CONFIG_TEMPLATE = REPO_ROOT / "config" / "config.cfg"
LAYERS_TEMPLATE = REPO_ROOT / "config" / "layers.cfg"
MAX_CONCURRENT = int(os.environ.get("REGOLIT_WEB_CONCURRENCY", "4"))
MAX_QUEUE = int(os.environ.get("REGOLIT_WEB_QUEUE", "12"))
RUN_TIMEOUT = float(os.environ.get("REGOLIT_WEB_TIMEOUT", "43200"))
MAX_RUNS_KEPT = int(os.environ.get("REGOLIT_WEB_MAX_RUNS", "500"))
MAX_RUNS_DISK_GB = float(os.environ.get("REGOLIT_WEB_MAX_DISK_GB", "200"))
# Threads per model run (the slope relaxation is OpenMP-parallel): share the physical cores among
# the concurrent runs unless the environment says otherwise.
os.environ.setdefault("OMP_NUM_THREADS", str(max(1, ((os.cpu_count() or 2) // 2) // MAX_CONCURRENT)))
AUTH_USER = os.environ.get("REGOLIT_WEB_USER", "")
AUTH_PASSWORD = os.environ.get("REGOLIT_WEB_PASSWORD", "")

# Hard limits protecting the server, sized for the compute machine (16 cores, 125 GB, 20 TB data disk).
# Measured with 8 threads over 100 Ma: 2000 x 2000 cells run in 18 s and need 0.8 GB of memory,
# 4000 x 4000 cells run in 3 min and need 3.3 GB (memory grows roughly linearly with the cell count,
# about 0.2 GB per million cells). 8000 x 8000 cells therefore fit in about 15 GB, and four of them
# at a time in 60 GB; the wall-clock is bounded by the run timeout. One million small craters take
# about a minute. Every saved step writes seven maps of the output grid, so the number of steps is
# bounded through the total output size rather than by a fixed count.
MAX_GRID_SIDE = 8000
MAX_GRID_CELLS = MAX_GRID_SIDE * MAX_GRID_SIDE
MAX_STEPS = 500
MAX_OUTPUT_BYTES = 500 * 1024 * 1024 * 1024
MAPS_PER_STEP = 7
MAX_IMPACTS = 20_000_000
MAX_LAYER_ROWS = 40

RUN_ID_PATTERN = re.compile(r"^[0-9]{8}-[0-9]{6}-[0-9a-f]{6}$")

# Figure theme (matches the page's dark palette).
THEME = {"figure": "#141a21", "axes": "#0d1116", "ink": "#e4e8ee", "muted": "#8d98a8", "grid": "#2a3340", "accent": "#e0a458"}
# Composition colors of the layered cross-section: regolith, ice, soot.
COMPOSITION_COLORS = np.array([[0.62, 0.56, 0.48], [0.60, 0.84, 1.00], [0.10, 0.10, 0.12]])

# Map figure geometry, fixed so the page can overlay the cross-section line on the image. The map
# axes are exactly square in pixels (the domain is square), so the image fills them.
MAP_DPI = 110
MAP_LAYOUTS = {
    "with_colorbar": {"size": (5.7, 4.8), "axes": (0.125, 0.11, 0.81 * 4.8 / 5.7, 0.81), "cbar": (0.845, 0.11, 0.03, 0.81)},
    "without_colorbar": {"size": (4.8, 4.8), "axes": (0.14, 0.11, 0.81, 0.81), "cbar": None},
}
WIDE_FIG_SIZE = (8.6, 2.9)   # cross-sections and histograms

# Parameters the UI exposes: name, label, unit, min, max, kind, description. Values not listed
# here stay at the repository defaults (config/config.cfg).
PARAMETERS: List[Dict] = [
    dict(group="Domain", name="regionWidth", label="Region width", unit="m", min=50, max=1_000_000, kind="number",
         description="Side of the square, periodic domain; bounded only through the cell count."),
    dict(group="Domain", name="resolution", label="Resolution", unit="m/pixel", min=0.5, max=1000, kind="number",
         description="Cell size. Cells = (width / resolution)^2, at most {0} x {0}.".format(MAX_GRID_SIDE)),
    dict(group="Domain", name="downsamplingResolution", label="Map output resolution", unit="m/pixel", min=0.5, max=5000, kind="number",
         description="Maps are averaged to this resolution before saving (>= resolution). Cross-sections use the full resolution."),
    dict(group="Time", name="endTime", label="Duration", unit="Ma", min=0.1, max=4500, kind="number",
         description="Simulated time."),
    dict(group="Time", name="printTimeStep", label="Output interval", unit="Ma", min=0.01, max=4500, kind="number",
         description="Time between saved steps; up to 500 steps, bounded by the total output size."),
    dict(group="Time", name="histogramBinsPerDecade", label="Histogram bins per decade", unit="", min=4, max=200, kind="int",
         description="Log bins per decade of the size and depth histograms behind the figures."),
    dict(group="Time", name="randomSeed", label="Random seed", unit="", min=0, max=2**31 - 1, kind="int",
         description="Seed of the impactor sequence."),
    dict(group="Impactors", name="minimumImpactorDiameter", label="Minimum impactor diameter", unit="m", min=0.02, max=50, kind="number",
         description="Smallest impactor simulated."),
    dict(group="Impactors", name="fluxConstant_c", label="Flux constant c", unit="m^-2 Ma^-1", min=1e-14, max=1e-2, kind="number", auto="production_function",
         description="N(>1 m) impactors per m^2 per Ma. Set by the production function unless 'Power law' is chosen."),
    dict(group="Impactors", name="slope_b", label="Slope b", unit="", min=1.2, max=4.5, kind="number", auto="production_function",
         description="Cumulative slope: N(>d) ~ d^-b. Set by the production function unless manual."),
    dict(group="Impactors", name="earthFluxRatioCoefficient", label="Flux multiplier", unit="", min=0.001, max=100, kind="number",
         description="Scales the production function, e.g. to another body or epoch (1 = as published)."),
    dict(group="Impactors", name="impactorDensity", label="Impactor density", unit="kg/m^3", min=300, max=9000, kind="number"),
    dict(group="Target", name="g", label="Gravity", unit="m/s^2", min=0.01, max=30, kind="number", body=True),
    dict(group="Target", name="meanImpactVelocity", label="Impact velocity", unit="m/s", min=500, max=80000, kind="number", body=True),
    dict(group="Target", name="targetDensity", label="Target density", unit="kg/m^3", min=300, max=6000, kind="number", body=True),
    dict(group="Target", name="k1", label="Scaling constant K1", unit="", min=0.01, max=2, kind="number",
         description="Holsapple (1993) crater-volume scaling constant (0.132 for soils and regolith)."),
    dict(group="Target", name="k2", label="Strength constant K2", unit="", min=0, max=2, kind="number", body=True,
         description="Holsapple (1993) strength-regime constant (0.26 for dry soil and lunar regolith; 1 in the original model)."),
    dict(group="Target", name="mu", label="Scaling exponent mu", unit="", min=0.3, max=0.7, kind="number"),
    dict(group="Target", name="Ybar", label="Effective strength", unit="Pa", min=0, max=1e9, kind="number", body=True,
         description="Target strength; matters for the smallest craters. Set by the body preset from Williams et al. (2014)."),
    dict(group="Target", name="angleOfRepose", label="Angle of repose", unit="deg", min=5, max=80, kind="number",
         description="Slopes steeper than this fail at every output step. Not changed by the body preset."),
    dict(group="Craters", name="craterProfileType", label="Cavity shape", unit="", min=1, max=2, kind="choice",
         choices={"1": "parabolic", "2": "bowl (spherical cap)"}),
    dict(group="Craters", name="depthToDiameter", label="Depth / diameter", unit="", min=0.02, max=0.5, kind="number"),
    dict(group="Craters", name="rimToDiameter", label="Rim height / diameter", unit="", min=0, max=0.2, kind="number"),
    dict(group="Craters", name="rimDropoffExponent", label="Rim dropoff exponent", unit="", min=1, max=8, kind="number"),
    dict(group="Craters", name="isEmplaceEjecta", label="Ejecta blanket", unit="", min=0, max=1, kind="bool"),
    dict(group="Craters", name="ejectaSpread", label="Ejecta extent", unit="radii", min=2, max=16, kind="int"),
    dict(group="Craters", name="numberOfZModelShells", label="Z-model shells", unit="", min=10, max=500, kind="int"),
    dict(group="Craters", name="minimumLayerThickness", label="Minimum layer thickness", unit="m", min=0, max=1, kind="number",
         description="Thinner deposits are mixed into the surface layer."),
    dict(group="Secondary craters", name="isEmplaceSecondaries", label="Secondary craters", unit="", min=0, max=1, kind="bool", toggle=True,
         description="Fragments of the fast ejecta of every primary form craters where they land (Z-model launch speeds and ranges)."),
    dict(group="Secondary craters", name="secondaryMinimumVelocity", label="Minimum landing speed", unit="m/s", min=1, max=3000, kind="number",
         description="Ejecta landing slower than this only builds the blanket. Faster fragments form a crater where it would be deeper than the blanket there."),
    dict(group="Secondary craters", name="secondaryLargestFraction", label="Largest secondary / primary diameter", unit="", min=0.005, max=0.3, kind="number",
         description="Diameter of the largest secondary, made by the largest fragment landing about three radii out (about 0.05 on the Moon)."),
    dict(group="Secondary craters", name="secondaryVelocityExponent", label="Fragment size–speed exponent", unit="", min=0, max=3, kind="number",
         description="Faster launch annuli have smaller largest fragments, (speed / speed at three radii)^-exponent; 1 for spallation scaling."),
    dict(group="Secondary craters", name="slope_secondaries", label="Fragment size-distribution slope", unit="", min=1.5, max=8, kind="number",
         description="Cumulative slope: N(>L) = (L_max / L)^slope, one fragment at the largest size."),
    dict(group="Secondary craters", name="secondaryDepthToDiameter", label="Depth / diameter", unit="", min=0.02, max=0.5, kind="number",
         description=""),
    dict(group="Secondary craters", name="maximumSecondariesPerPrimary", label="Secondaries per primary (budget)", unit="", min=100, max=1_000_000, kind="int",
         description="Only the largest N secondaries of a primary are formed; the smaller fragments count as ejecta."),
    dict(group="Secondary craters", name="isEmplaceDistantSecondaries", label="Distant primaries", unit="", min=0, max=1, kind="bool", toggle=True,
         description="Also sample the primaries that form outside the domain and form the fragments they send in."),
    dict(group="Secondary craters", name="secondaryMaximumRange", label="Distant primaries out to", unit="m", min=1000, max=3_000_000, kind="number", depends="isEmplaceDistantSecondaries",
         description="Distance beyond the domain edge out to which distant primaries are sampled."),
    dict(group="Seismic shaking", name="isSeismicShaking", label="Seismic shaking", unit="", min=0, max=1, kind="bool", toggle=True,
         description="Every impact shakes its surroundings: downslope diffusion where the peak acceleration exceeds the threshold, plus slope collapse (Richardson 2005, 2009)."),
    dict(group="Seismic shaking", name="seismicEfficiency", label="Seismic efficiency", unit="", min=1e-8, max=1e-2, kind="number",
         description="Fraction of the impact energy radiated as seismic energy (10^-6 to 10^-4 in the literature)."),
    dict(group="Seismic shaking", name="Cs", label="Diffusion constant Cs", unit="", min=1e-6, max=1, kind="number",
         description="Dose per impact K = Cs v^0.5 D^1 g^0.5 / l^0.5 (m^2), with D the impactor diameter and l the distance (Richardson 2009, Eq. 32)."),
    dict(group="Seismic shaking", name="seismicAccelerationThreshold", label="Acceleration threshold", unit="g", min=0.05, max=10, kind="number",
         description="Regolith moves where the peak seismic acceleration exceeds this multiple of the surface gravity."),
    dict(group="Seismic shaking", name="prim_seis_freq", label="Seismic frequency", unit="Hz", min=1, max=100, kind="number",
         description="Primary frequency of the impact seismic signal (10 to 20 Hz)."),
    dict(group="Seismic shaking", name="Q_factor", label="Quality factor Q", unit="", min=10, max=10000, kind="number",
         description="Seismic attenuation; 1000 to 2000 for the fractured lunar crust."),
    dict(group="Seismic shaking", name="seis_mean_free", label="Scattering mean free path", unit="m", min=10, max=20000, kind="number",
         description="Sets the seismic diffusivity together with the wave velocity (125 to 2000 m on the Moon)."),
    dict(group="Seismic shaking", name="seis_wave_vel", label="Seismic velocity", unit="m/s", min=100, max=10000, kind="number",
         description="P-wave velocity in competent rock."),
    dict(group="Ice and soot", name="isVolatiles", label="Ice and soot processes", unit="", min=0, max=1, kind="bool", toggle=True,
         description="Periodic ice deposition, sublimation and the loss of ice and soot from ejecta. Off: ice and soot in the initial layers are passive tracers."),
    dict(group="Ice and soot", name="iceEmplacementInterval", label="Ice deposition interval", unit="Ma", min=0.01, max=4500, kind="number"),
    dict(group="Ice and soot", name="iceEmplacementThickness", label="Ice deposition thickness", unit="m", min=0, max=10, kind="number",
         description="Ice added on the surface at every interval (0 disables)."),
    dict(group="Ice and soot", name="sublimationInterval", label="Sublimation interval", unit="Ma", min=0.01, max=4500, kind="number",
         description="Reserved: the sublimation model is not implemented yet, so this has no effect."),
    dict(group="Ice and soot", name="sublimationThickness", label="Sublimation thickness", unit="m", min=0, max=10, kind="number",
         description="Reserved: the sublimation model is not implemented yet, so this has no effect."),
    dict(group="Ice and soot", name="ejectaVolatileRetention", label="Ice retained in ejecta", unit="fraction", min=0, max=1, kind="number",
         description="Fraction of the excavated ice that lands in the ejecta blanket."),
    dict(group="Ice and soot", name="ejectaSootRetention", label="Soot retained in ejecta", unit="fraction", min=0, max=1, kind="number",
         description="Fraction of the excavated soot that lands in the ejecta blanket."),
    # Test mode (not shown in the setup form; set by the Tests tab and kept by "Edit & run again"):
    dict(group="Test", name="testCraterDiameter", label="Test crater diameter", unit="m", min=0, max=100_000, kind="number", hidden=True,
         description="Form one crater of this final diameter instead of the random population (0 = off)."),
    dict(group="Test", name="testCraterX", label="Test crater x", unit="m", min=-500_000, max=500_000, kind="number", hidden=True, description=""),
    dict(group="Test", name="testCraterY", label="Test crater y", unit="m", min=-500_000, max=500_000, kind="number", hidden=True, description=""),
    dict(group="Subsurface", name="initialThickness", label="Basement thickness", unit="m", min=1, max=10000, kind="number", auto="basement",
         description="Automatic: the initial layers plus three times the depth of the largest expected crater. Untick to override."),
    dict(group="Subsurface", name="depthToIntegrate", label="Integration depth", unit="m", min=0.005, max=100, kind="number",
         description="Depth of the integrated-composition maps."),
]
PARAMETER_INDEX = {p["name"]: p for p in PARAMETERS}
MAP_KINDS = {
    "shaded_relief": ("Shaded relief", "", "gray"),
    "elevation": ("Elevation", "m", "terrain"),
    "surface_soot": ("Surface soot fraction", "", "viridis"),
    "surface_ice": ("Surface ice fraction", "", "viridis"),
    "surface_regolith": ("Surface regolith fraction", "", "viridis"),
    "integrated_soot": ("Soot fraction, integrated", "", "viridis"),
    "integrated_ice": ("Ice fraction, integrated", "", "viridis"),
    "integrated_regolith": ("Regolith fraction, integrated", "", "viridis"),
}

# ----------------------------------------------------------------------------------------------
# App, auth, concurrency
# ----------------------------------------------------------------------------------------------
security = HTTPBasic(auto_error=False)


def require_auth(credentials: Optional[HTTPBasicCredentials] = Depends(security)) -> None:
    if not AUTH_USER:
        return
    ok = credentials is not None and secrets.compare_digest(credentials.username, AUTH_USER) \
        and secrets.compare_digest(credentials.password, AUTH_PASSWORD)
    if not ok:
        raise HTTPException(status_code=401, detail="Authentication required", headers={"WWW-Authenticate": "Basic realm=REGOLIT"})


@asynccontextmanager
async def lifespan(_: FastAPI):
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    requeue = recover_runs()
    prune_runs()
    if not BINARY.exists():
        print("WARNING: model binary not found at {}; run `make` first.".format(BINARY), file=sys.stderr)
    for run_id, request in requeue:
        start_background(run_id, request["overrides"], request["layers"], request["presets"])
    watchdog = asyncio.create_task(watch_orphans())
    yield
    watchdog.cancel()


app = FastAPI(title="REGOLIT", docs_url=None, redoc_url=None, dependencies=[Depends(require_auth)], lifespan=lifespan)
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

RUN_SEMAPHORE = asyncio.Semaphore(MAX_CONCURRENT)
WAITING = {"count": 0}
PLOT_LOCK = threading.Lock()
SUBSURFACE_CACHE: "OrderedDict[str, Subsurface]" = OrderedDict()
SUBSURFACE_CACHE_SIZE = 4


class RunRequest(BaseModel):
    parameters: Dict[str, object] = Field(default_factory=dict)   # coerced and range-checked in validate() / estimate_for()
    layers: Optional[List[List[float]]] = None   # rows of (class, thickness, regolith, ice, soot), bottom-up
    presets: Optional[Dict[str, object]] = None  # UI choices kept with the run: body, production_function, basement_auto


class SettingsRequest(RunRequest):
    name: str = Field(min_length=1, max_length=80)


PRESET_KEYS = {"body": str, "production_function": str, "basement_auto": bool, "test": str}

# Predefined test scenarios (the Tests tab of the setup page). Parameters not listed keep the defaults.
TESTS: List[Dict] = [
    dict(id="single_crater_secondaries", title="A 1 km crater with secondaries",
         summary="10 km domain at 4 m/pixel. One 1 km crater at the centre, Moon parameters, with its ejecta blanket "
                 "and the secondary craters formed by its fragments. Two steps: before and after.",
         parameters=dict(regionWidth=10000, resolution=4, downsamplingResolution=4, endTime=1, printTimeStep=1,
                         initialThickness=1000, isEmplaceSecondaries=1, isEmplaceDistantSecondaries=0,
                         testCraterDiameter=1000, testCraterX=0, testCraterY=0, angleOfRepose=35,
                         **{key: scaling.BODIES["moon"][key] for key in scaling.BODY_FIELDS}),
         layers=[[0, 20, 1, 0, 0]],   # regolith only
         presets=dict(body="moon", production_function="power_law", basement_auto=False, test="single_crater_secondaries")),
]


def clean_presets(presets: Optional[Dict[str, object]]) -> Dict[str, object]:
    out: Dict[str, object] = {}
    for key, kind in PRESET_KEYS.items():
        if presets and key in presets:
            value = presets[key]
            if kind is bool:
                out[key] = bool(value)
            elif isinstance(value, str) and len(value) <= 40:
                out[key] = value
    return out


def layers_total_thickness(rows: List[List[float]]) -> float:
    """Thickness of the initial layer stack the model will use (the second class when there are several, else the first)."""
    blocks: List[List[List[float]]] = []
    previous = None
    for row in rows:
        if row[0] != previous:
            blocks.append([])
            previous = row[0]
        blocks[-1].append(row)
    block = blocks[1] if len(blocks) > 1 else blocks[0]
    return float(sum(r[1] for r in block))


def estimate_for(request: RunRequest) -> Dict:
    """Derived quantities for a (possibly half-edited) parameter set: values are coerced into range, never rejected."""
    effective = default_parameters()
    for name, value in request.parameters.items():
        spec = PARAMETER_INDEX.get(name)
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if spec is None or not math.isfinite(value):
            continue
        effective[name] = min(max(value, spec["min"]), spec["max"])
    rows = request.layers if request.layers else default_layers()
    try:
        total = layers_total_thickness([[float(v) for v in r] for r in rows if len(r) == 5])
    except (TypeError, ValueError, IndexError):
        total = layers_total_thickness(default_layers())
    presets = clean_presets(request.presets)
    key = presets.get("production_function", "power_law")
    if key not in scaling.PRODUCTION_FUNCTIONS or not scaling.PRODUCTION_FUNCTIONS[key].get("available"):
        key = "power_law"
    try:
        return scaling.estimate(effective, total, key, presets.get("body"))
    except (ValueError, ZeroDivisionError, OverflowError) as error:
        raise HTTPException(400, "cannot estimate: {}".format(error))


# ----------------------------------------------------------------------------------------------
# Helpers: configuration and validation
# ----------------------------------------------------------------------------------------------
def default_parameters() -> Dict[str, float]:
    config = read_config(CONFIG_TEMPLATE)
    return {p["name"]: config.get(p["name"], 0.0) for p in PARAMETERS}


def default_layers() -> List[List[float]]:
    return [list(map(float, row)) for row in read_layers(LAYERS_TEMPLATE)]


def validate(request: RunRequest) -> Tuple[Dict[str, float], Optional[List[List[float]]]]:
    """Check names, ranges and derived limits; returns the config overrides and the layer rows."""
    overrides: Dict[str, float] = {}
    for name, value in request.parameters.items():
        spec = PARAMETER_INDEX.get(name)
        if spec is None:
            raise HTTPException(400, "unknown parameter {!r}".format(name))
        try:
            value = float(value)
        except (TypeError, ValueError):
            raise HTTPException(400, "{}: not a number".format(name))
        if not math.isfinite(value) or value < spec["min"] or value > spec["max"]:
            raise HTTPException(400, "{}: must be between {:g} and {:g}".format(spec["label"], spec["min"], spec["max"]))
        if spec["kind"] in ("int", "bool", "choice"):
            value = int(round(value))
        overrides[name] = value

    effective = default_parameters()
    effective.update(overrides)
    width, res = effective["regionWidth"], effective["resolution"]
    side = round(width / res)
    cells = side ** 2
    if cells > MAX_GRID_CELLS:
        raise HTTPException(400, "the grid would have {0:,} x {0:,} cells; the limit is {1:,} x {1:,}. Increase the resolution or shrink the region.".format(int(side), MAX_GRID_SIDE))
    if effective["downsamplingResolution"] < res:
        overrides["downsamplingResolution"] = res
        effective["downsamplingResolution"] = res
    steps = math.ceil(effective["endTime"] / effective["printTimeStep"]) + 1
    if steps > MAX_STEPS:
        raise HTTPException(400, "{} output steps requested; the limit is {}. Increase the output interval.".format(steps, MAX_STEPS))
    output_cells = round(width / effective["downsamplingResolution"]) ** 2
    output_bytes = steps * output_cells * MAPS_PER_STEP * 8
    if output_bytes > MAX_OUTPUT_BYTES:
        raise HTTPException(400, "the run would write {:.1f} GB of maps ({} steps of {:,} x {:,} cells); the limit is {:.0f} GB. Increase the output interval or the map output resolution.".format(
            output_bytes / 2**30, steps, int(round(width / effective["downsamplingResolution"])), int(round(width / effective["downsamplingResolution"])), MAX_OUTPUT_BYTES / 2**30))
    impacts = effective["fluxConstant_c"] * effective["minimumImpactorDiameter"] ** (-effective["slope_b"]) * width ** 2 \
        * effective["endTime"] * effective["earthFluxRatioCoefficient"]
    if impacts > MAX_IMPACTS:
        raise HTTPException(400, "about {:,} impacts would be simulated; the limit is {:,}. Shorten the run, shrink the region or raise the minimum impactor diameter.".format(int(impacts), MAX_IMPACTS))
    overrides["isPrintSubsurface"] = 2   # layer stacks of the final state, for the cross-sections

    layers = None
    if request.layers is not None:
        if not 1 <= len(request.layers) <= MAX_LAYER_ROWS:
            raise HTTPException(400, "between 1 and {} layer rows are allowed".format(MAX_LAYER_ROWS))
        layers = []
        for row in request.layers:
            if len(row) != 5:
                raise HTTPException(400, "each layer row needs 5 values: class, thickness, regolith, ice, soot")
            cls, thickness, reg, ice, soot = [float(v) for v in row]
            if not (0 <= cls <= 100 and cls == int(cls)):
                raise HTTPException(400, "layer class must be an integer between 0 and 100")
            if not (0 < thickness <= 10000) or min(reg, ice, soot) < 0 or reg + ice + soot <= 0:
                raise HTTPException(400, "layer thickness must be positive and the composition non-negative with a positive sum")
            layers.append([int(cls), thickness, reg, ice, soot])
    return overrides, layers


# ----------------------------------------------------------------------------------------------
# Helpers: runs on disk
# ----------------------------------------------------------------------------------------------
def new_run_id() -> str:
    return time.strftime("%Y%m%d-%H%M%S") + "-" + secrets.token_hex(3)


def run_dir(run_id: str) -> Path:
    if not RUN_ID_PATTERN.match(run_id) or not (RUNS_DIR / run_id / "summary.json").exists():
        raise HTTPException(404, "unknown run")
    return RUNS_DIR / run_id


def directory_size(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def prune_runs() -> None:
    """Delete the oldest runs beyond the count limit and until the runs directory fits the disk budget."""
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    runs = sorted([p for p in RUNS_DIR.iterdir() if p.is_dir() and read_status(p) == "done"], key=lambda p: p.stat().st_mtime)
    for old in runs[: max(0, len(runs) - MAX_RUNS_KEPT)]:
        shutil.rmtree(old, ignore_errors=True)
    runs = runs[max(0, len(runs) - MAX_RUNS_KEPT):]
    sizes = {p: directory_size(p) for p in runs}
    total = sum(sizes.values())
    budget = MAX_RUNS_DISK_GB * 1024 ** 3
    for old in runs:
        if total <= budget:
            break
        shutil.rmtree(old, ignore_errors=True)
        total -= sizes[old]


def read_status(path: Path) -> str:
    try:
        return json.loads((path / "summary.json").read_text()).get("status", "done")
    except (OSError, ValueError):
        return "unknown"


def write_summary(run_id: str, summary: Dict) -> None:
    target = RUNS_DIR / run_id / "summary.json"
    tmp = target.with_name("summary.json.tmp")
    tmp.write_text(json.dumps(summary, indent=1))
    os.replace(tmp, target)


def queue_run(run_id: str, overrides: Dict[str, float], layers: Optional[List[List[float]]], presets: Optional[Dict[str, object]]) -> Dict:
    """Create the run directory with a provisional summary (status 'queued') so it shows up at once."""
    effective = default_parameters()
    effective.update(overrides)
    workdir = RUNS_DIR / run_id
    workdir.mkdir(parents=True, exist_ok=True)
    summary = {
        "id": run_id,
        "status": "queued",
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "queued_at": time.time(),
        "parameters": {p["name"]: effective.get(p["name"]) for p in PARAMETERS},
        "layers": layers if layers is not None else default_layers(),
        "presets": clean_presets(presets),
        "steps": [], "times": [],
    }
    (workdir / "request.json").write_text(json.dumps({"overrides": overrides, "layers": layers, "presets": presets}))
    write_summary(run_id, summary)
    return summary


def pid_alive(pid: Optional[int]) -> bool:
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError):
        return False
    return True


PROGRESS_PATTERN = re.compile(r"Progress: ([0-9.]+)%")


def read_tail(path: Path, max_bytes: int) -> str:
    """The last max_bytes of a text file (the model log can reach hundreds of MB)."""
    with open(path, "rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - max_bytes))
        data = handle.read()
    if size > max_bytes:
        data = data.split(b"\n", 1)[-1]   # drop the partial first line
    return data.decode(errors="replace")


def progress_of(run_id: str, summary: Dict) -> Dict:
    """Progress fields for an unfinished run, read from the model's log."""
    info = {"progress": None, "phase": "waiting for a free slot" if summary.get("status") == "queued" else "starting", "elapsed_s": None}
    if summary.get("started_at"):
        info["elapsed_s"] = round(time.time() - summary["started_at"], 1)
    log_path = RUNS_DIR / run_id / "log" / "log.txt"
    if summary.get("status") == "running" and log_path.exists():
        try:
            lines = read_tail(log_path, 64 * 1024).splitlines()
        except OSError:
            lines = []
        for line in reversed(lines):
            match = PROGRESS_PATTERN.search(line)
            if match:
                info["progress"] = float(match.group(1))
                break
        for line in reversed(lines):   # the latest non-warning line describes the current phase
            if "WARNING" not in line:
                info["phase"] = line.split("\t", 1)[-1].strip()
                break
    return info


def execute(run_id: str, overrides: Dict[str, float], layers: Optional[List[List[float]]], presets: Optional[Dict[str, object]] = None) -> Dict:
    """Run the model (blocking) in the run directory prepared by queue_run, then finalize it."""
    workdir = RUNS_DIR / run_id
    started = time.time()
    summary = load_summary(run_id)
    summary.update(status="running", started_at=started, pid=os.getpid())
    write_summary(run_id, summary)
    regolit.run(overrides, workdir, binary=BINARY, layers=layers if layers is not None else LAYERS_TEMPLATE,
                build_if_missing=False, quiet=True, timeout=RUN_TIMEOUT)
    return finalize(run_id, summary, started)


def finalize(run_id: str, summary: Dict, started: float) -> Dict:
    """Compact the layer stacks of a finished model run and write the final summary.json."""
    workdir = RUNS_DIR / run_id
    out = RegolitOutput(workdir)
    series = out.elevation_series()
    log_text = out.log()
    match = re.search(r"Number of craters in simulation: (\d+)", log_text)
    secondaries = re.search(r"Secondary craters: (\d+) from \d+ primaries inside the domain(?:, (\d+) from \d+ distant primaries)?", log_text)
    secondary_count = (int(secondaries.group(1)) + int(secondaries.group(2) or 0)) if secondaries else 0
    test_crater = 1 if float(out.config.get("testCraterDiameter", 0) or 0) > 0 else 0

    # Store the final layer stacks compactly (the raw file is large) and drop the raw file.
    has_subsurface = False
    for suffix in out.subsurface_steps():
        raw = out.output_dir / "subsurface_{}.out".format(suffix)
        if raw.exists():
            sub = regolit.read_subsurface(raw)
            sub.save(out.output_dir / "subsurface_{}.npz".format(suffix))
            raw.unlink()
        has_subsurface = True

    x_full, _ = out.full_resolution_coordinates()
    summary.update({
        "status": "done",
        "duration_s": round(time.time() - started, 2),
        "parameters": {p["name"]: out.config.get(p["name"]) for p in PARAMETERS},
        "steps": out.steps,
        "times": [float(t) for t in out.times],
        "grid": out.n,
        "resolution": out.resolution,
        "full_grid": int(len(x_full)),
        "full_resolution": float(out.config.get("resolution", out.resolution)),
        "extent": out.extent,
        "elevation_min": float(series.min()),
        "elevation_max": float(series.max()),
        "total_craters": (int(match.group(1)) + test_crater + secondary_count) if match else None,
        "primary_craters": (int(match.group(1)) + test_crater) if match else None,
        "secondary_craters": secondary_count,
        "visible_craters": int(len(out.craters()["x"])),
        "has_subsurface": has_subsurface,
    })
    summary.pop("started_at", None); summary.pop("queued_at", None); summary.pop("pid", None)
    (workdir / "cache").mkdir(exist_ok=True)
    write_summary(run_id, summary)
    (workdir / "request.json").unlink(missing_ok=True)
    return summary


def fail_run(run_id: str, error: str) -> None:
    try:
        summary = load_summary(run_id)
    except (OSError, ValueError, HTTPException):
        return
    summary.update(status="failed", error=error[:2000], duration_s=round(time.time() - summary.get("started_at", time.time()), 2))
    summary.pop("started_at", None); summary.pop("queued_at", None); summary.pop("pid", None)
    write_summary(run_id, summary)


def recover_runs() -> List[Tuple[str, Dict]]:
    """After a (re)start: running runs whose worker is alive continue on their own; running runs whose
    worker is gone are marked failed; queued runs are returned so they can be started again."""
    requeue = []
    for path in RUNS_DIR.iterdir() if RUNS_DIR.exists() else []:
        if not (path.is_dir() and RUN_ID_PATTERN.match(path.name)):
            continue
        status = read_status(path)
        if status == "running":
            summary = json.loads((path / "summary.json").read_text())
            if not pid_alive(summary.get("pid")):
                fail_run(path.name, "the worker process disappeared while this run was in progress")
        elif status == "queued":
            try:
                requeue.append((path.name, json.loads((path / "request.json").read_text())))
            except (OSError, ValueError):
                fail_run(path.name, "the server was restarted before this run started")
    return requeue


def list_runs() -> List[Dict]:
    """Short descriptions of the stored runs, newest first."""
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    runs = []
    for path in RUNS_DIR.iterdir():
        summary_path = path / "summary.json"
        if not (path.is_dir() and RUN_ID_PATTERN.match(path.name) and summary_path.exists()):
            continue
        try:
            summary = json.loads(summary_path.read_text())
        except ValueError:
            continue
        item = {key: summary.get(key) for key in ("id", "status", "error", "created", "duration_s", "grid", "resolution", "total_craters", "visible_craters")}
        item["status"] = item["status"] or "done"
        item["test"] = (summary.get("presets") or {}).get("test")
        item["parameters"] = {key: summary.get("parameters", {}).get(key) for key in ("regionWidth", "resolution", "endTime", "randomSeed")}
        runs.append(item)
    runs.sort(key=lambda item: item["id"], reverse=True)
    return runs


def load_summary(run_id: str) -> Dict:
    return json.loads((run_dir(run_id) / "summary.json").read_text())


def require_done(run_id: str) -> Dict:
    summary = load_summary(run_id)
    if summary.get("status", "done") != "done":
        raise HTTPException(409, "run {} is {}".format(run_id, summary.get("status")))
    return summary


def output_of(run_id: str) -> RegolitOutput:
    require_done(run_id)
    return RegolitOutput(run_dir(run_id))


def subsurface_of(run_id: str, out: RegolitOutput) -> Subsurface:
    """Final layer stacks of a run, cached in memory for the interactive cross-sections."""
    with PLOT_LOCK:
        if run_id in SUBSURFACE_CACHE:
            SUBSURFACE_CACHE.move_to_end(run_id)
            return SUBSURFACE_CACHE[run_id]
    sub = out.subsurface(-1)
    with PLOT_LOCK:
        SUBSURFACE_CACHE[run_id] = sub
        while len(SUBSURFACE_CACHE) > SUBSURFACE_CACHE_SIZE:
            SUBSURFACE_CACHE.popitem(last=False)
    return sub


def step_index(summary: Dict, step: int) -> int:
    if not 0 <= step < len(summary["steps"]):
        raise HTTPException(404, "step out of range")
    return step


# ----------------------------------------------------------------------------------------------
# Helpers: figures
# ----------------------------------------------------------------------------------------------
def themed_figure(size: Tuple[float, float], dpi: int = MAP_DPI):  # noqa: D103
    fig = plt.figure(figsize=size, dpi=dpi, facecolor=THEME["figure"])
    return fig


def style_axes(ax) -> None:
    ax.set_facecolor(THEME["axes"])
    for spine in ax.spines.values():
        spine.set_color(THEME["grid"])
    ax.tick_params(colors=THEME["muted"], labelsize=8.5)
    ax.xaxis.label.set_color(THEME["ink"])
    ax.yaxis.label.set_color(THEME["ink"])
    ax.title.set_color(THEME["ink"])


def map_geometry() -> Dict:
    """Pixel position of the map axes inside the PNG (origin top-left) for both layouts, for the page overlay."""
    geometry = {}
    for name, layout in MAP_LAYOUTS.items():
        width_px = layout["size"][0] * MAP_DPI
        height_px = layout["size"][1] * MAP_DPI
        left, bottom, width, height = layout["axes"]
        geometry[name] = {
            "width": width_px, "height": height_px,
            "x0": left * width_px, "x1": (left + width) * width_px,
            "y0": (1 - bottom - height) * height_px, "y1": (1 - bottom) * height_px,
        }
    return geometry


def hillshade(z: np.ndarray, resolution: float, azimuth: float, altitude: float, reference: Optional[np.ndarray] = None) -> np.ndarray:
    """Hillshade in [0, 1]; the grey stretch is taken from the reference elevation (the final state) so
    that all steps of a run share the same contrast."""
    source = LightSource(azdeg=azimuth, altdeg=altitude)
    shade = source.hillshade(z, vert_exag=1.0, dx=resolution, dy=resolution)
    ref = shade if reference is None else source.hillshade(reference, vert_exag=1.0, dx=resolution, dy=resolution)
    lo, hi = np.percentile(ref, [0.5, 99.5])
    if hi - lo < 1e-6:
        return shade
    return np.clip((shade - lo) / (hi - lo), 0.0, 1.0)


def atomic_savefig(fig, target: Path) -> None:
    """Save into a temporary file and move it into place, so a concurrent reader never sees a partial file."""
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    fig.savefig(str(temporary), format="png", facecolor=fig.get_facecolor())
    os.replace(str(temporary), str(target))


def not_found_if_deleted(function, *args):
    """Run a renderer; a run deleted while rendering surfaces as 404 instead of a server error."""
    try:
        return function(*args)
    except FileNotFoundError:
        raise HTTPException(404, "the run is no longer available")


def render_map(run_id: str, kind: str, step: int, azimuth: float = 315.0, altitude: float = 25.0) -> Path:
    summary = require_done(run_id)
    step = step_index(summary, step)
    if kind not in MAP_KINDS:
        raise HTTPException(404, "unknown map kind")
    suffix = "_{:g}_{:g}".format(azimuth, altitude) if kind == "shaded_relief" else ""
    target = run_dir(run_id) / "cache" / "map_{}_{}{}.png".format(kind, step, suffix)
    if target.exists():
        return target
    out = output_of(run_id)
    label, unit, cmap = MAP_KINDS[kind]
    show_colorbar = True
    if kind == "elevation":
        data, vmin, vmax = out.elevation(step), summary["elevation_min"], summary["elevation_max"]
    elif kind == "shaded_relief":
        reference = out.elevation(-1) if step != len(summary["steps"]) - 1 else None
        data, vmin, vmax = hillshade(out.elevation(step), summary["resolution"], azimuth, altitude, reference), 0.0, 1.0
        show_colorbar = False
        label = "Shaded relief (sun from {:g} deg, {:g} deg high)".format(azimuth, altitude)
    else:
        where, species = kind.split("_")
        data = out.surface_fraction(species, step) if where == "surface" else out.integrated_fraction(species, step)
        vmin, vmax = 0.0, 1.0
        if where == "integrated":
            label = "{} (top {:g} m)".format(label, summary["parameters"].get("depthToIntegrate", float("nan")))
    layout = MAP_LAYOUTS["with_colorbar" if show_colorbar else "without_colorbar"]
    with PLOT_LOCK:
        fig = themed_figure(layout["size"])
        ax = fig.add_axes(layout["axes"])
        style_axes(ax)
        image = ax.imshow(data, origin="lower", extent=summary["extent"], cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title("{}, t = {:g} Ma".format(label, summary["times"][step]), fontsize=10)
        if show_colorbar:
            cax = fig.add_axes(layout["cbar"])
            colorbar = fig.colorbar(image, cax=cax)
            colorbar.set_label(unit, color=THEME["ink"])
            colorbar.ax.tick_params(colors=THEME["muted"], labelsize=8.5)
            colorbar.outline.set_edgecolor(THEME["grid"])
        atomic_savefig(fig, target)
        plt.close(fig)
    return target


def render_section(run_id: str, axis: str, at: float, depth: float, step: int) -> bytes:
    """Layered subsurface along a west-east (axis x, at y = at) or south-north (axis y, at x = at) line."""
    summary = require_done(run_id)
    step = step_index(summary, step)
    if axis not in ("x", "y"):
        raise HTTPException(400, "axis must be x or y")
    if not summary.get("has_subsurface"):
        raise HTTPException(404, "this run has no layer stacks")
    out = output_of(run_id)
    sub = subsurface_of(run_id, out)
    x_full, y_full = out.full_resolution_coordinates()
    half = summary["parameters"]["regionWidth"] / 2
    at = float(min(max(at, -half), half))
    positions, stacks, surface = sub.section(x_full, y_full, axis, at)
    cell = summary["full_resolution"]

    floor = float(surface.min() - depth)
    ceiling = float(surface.max() + 0.08 * (surface.max() - floor + 1e-9))
    polygons, colors = [], []
    for position, stack, top in zip(positions, stacks, surface):
        z_top = top
        for k in range(stack.shape[0] - 1, -1, -1):
            z_bottom = floor if k == 0 else z_top - stack[k, 0]
            z_bottom = max(z_bottom, floor)
            if z_top > floor:
                polygons.append([(position - cell / 2, z_bottom), (position + cell / 2, z_bottom), (position + cell / 2, z_top), (position - cell / 2, z_top)])
                colors.append(np.clip(stack[k, 1:4] @ COMPOSITION_COLORS, 0, 1))
            z_top = z_bottom
            if z_top <= floor:
                break

    # Surface of the selected step (from the output maps) for comparison with the final state:
    other = None
    if step != len(summary["steps"]) - 1:
        z = out.elevation(step)
        if axis == "x":
            j = int(np.argmin(np.abs(out.y - at)))
            other = (out.x, z[j, :])
        else:
            i = int(np.argmin(np.abs(out.x - at)))
            other = (out.y, z[:, i])

    with PLOT_LOCK:
        fig = themed_figure(WIDE_FIG_SIZE)
        ax = fig.add_axes((0.08, 0.19, 0.895, 0.69))
        style_axes(ax)
        ax.add_collection(PolyCollection(polygons, facecolors=colors, edgecolors="none"))
        ax.plot(positions, surface, color=THEME["ink"], lw=0.9)
        handles = [Patch(facecolor=COMPOSITION_COLORS[k], edgecolor="none", label=name) for k, name in enumerate(["regolith", "ice", "soot"])]
        if other is not None:
            line, = ax.plot(other[0], other[1], color=THEME["accent"], lw=1.0, ls="--", label="surface at t = {:g} Ma".format(summary["times"][step]))
            handles.append(line)
        ax.set_xlim(positions[0] - cell / 2, positions[-1] + cell / 2)
        ax.set_ylim(floor, ceiling)
        ax.set_xlabel("{} [m]".format("x" if axis == "x" else "y"))
        ax.set_ylabel("elevation [m]")
        line_name = "y = {:.1f} m (west-east)".format(at) if axis == "x" else "x = {:.1f} m (south-north)".format(at)
        ax.set_title("Subsurface layering along {}, final state (t = {:g} Ma)".format(line_name, summary["times"][-1]), fontsize=9.5)
        legend = ax.legend(handles=handles, loc="lower right", fontsize=8, facecolor=THEME["figure"], edgecolor=THEME["grid"], labelcolor=THEME["ink"], ncol=len(handles))
        legend.get_frame().set_alpha(0.9)
        ax.grid(alpha=0.15, color=THEME["muted"])
        buffer = __import__("io").BytesIO()
        fig.savefig(buffer, format="png", facecolor=fig.get_facecolor())
        plt.close(fig)
    return buffer.getvalue()


FIGURE_NAMES = ("impactors", "craters", "rplot", "turnover")
FIGURE_ALIASES = {"sizes": "craters", "depths": "turnover"}


def cumulative_from_histogram(bins: np.ndarray, counts: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """N(>= bin edge) from per-bin counts; only the edges from the first to the last populated bin."""
    populated = np.nonzero(counts)[0]
    if not len(populated):
        return np.array([]), np.array([])
    cumulative = np.cumsum(counts[::-1])[::-1]
    lo, hi = populated[0], populated[-1] + 1
    return bins[lo:hi + 1], np.append(cumulative[lo:hi], 0)[: hi + 1 - lo]


def bins_per_decade(bins: np.ndarray) -> float:
    return 1.0 / np.log10(bins[1] / bins[0]) if len(bins) > 1 else 20.0


def rebin(bins: np.ndarray, counts: np.ndarray, factor: int) -> Tuple[np.ndarray, np.ndarray]:
    """Aggregate `factor` consecutive bins (counts[k] belongs to [bins[k], bins[k+1]))."""
    factor = max(1, factor)
    n = (len(bins) - 1) // factor
    edges = bins[: n * factor + 1 : factor]
    merged = counts[: n * factor].reshape(n, factor).sum(axis=1)
    return edges, merged


_TURNOVER_TABLE: Dict[float, Tuple[np.ndarray, np.ndarray]] = {}


def poisson_lambda_for_count(n: int, probability: float) -> float:
    """The Poisson mean at which n or more events occur with the given probability (Gault et al. 1974, Table 1)."""
    from math import exp, lgamma, log, sqrt
    if n > 200:   # normal approximation with continuity correction; z is the (1 - p) quantile of the standard normal
        z = {0.5: 0.0, 0.9: -1.2816, 0.99: -2.3263, 0.1: 1.2816}.get(probability)
        if z is not None:
            root = (-z + sqrt(z * z + 4 * (n - 0.5))) / 2
            return root * root

    def survival(mu: float) -> float:   # P(N >= n)
        cdf = sum(exp(k * log(mu) - mu - lgamma(k + 1)) for k in range(n))
        return 1 - cdf

    lo, hi = 0.0, max(10.0, 4.0 * n)
    while survival(hi) < probability:
        hi *= 2
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if survival(mid) < probability:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def turnover_count(lam: np.ndarray, probability: float) -> np.ndarray:
    """Number of turnovers n reached with the given probability for Poisson means lam, interpolated smoothly
    between integer n as in Gault et al. (1974, Fig. 2); below the lambda of n = 1 the result is < 1."""
    if probability not in _TURNOVER_TABLE:
        ns = sorted(set(list(range(1, 11)) + [int(round(10 ** e)) for e in np.arange(1.0, 6.01, 0.125)]))
        lams = np.array([poisson_lambda_for_count(n, probability) for n in ns])
        _TURNOVER_TABLE[probability] = (np.log(lams), np.log(np.array(ns, dtype=float)))
    log_lams, log_ns = _TURNOVER_TABLE[probability]
    lam = np.asarray(lam, dtype=float)
    out = np.zeros_like(lam)
    ok = lam > 0
    log_lam = np.log(np.where(ok, lam, 1.0))
    out[ok] = np.exp(np.interp(log_lam[ok], log_lams, log_ns))
    # below n = 1 continue with the slope of the first segment so the curve leaves the plot smoothly
    small = ok & (log_lam < log_lams[0])
    slope = (log_ns[1] - log_ns[0]) / (log_lams[1] - log_lams[0])
    out[small] = np.exp(log_ns[0] + slope * (log_lam[small] - log_lams[0]))
    large = ok & (log_lam > log_lams[-1])
    out[large] = lam[large]
    return out
    return out


def plot_size_series(axes, series, area_km2: float, exact: Optional[np.ndarray] = None, exact_label: str = "") -> None:
    """Cumulative N(>D) per km² on axes[0] and counts per bin on axes[1] for (label, color, style, width, bins, counts) series."""
    lo, hi = np.inf, 0.0
    for label, color, style, width, bins, counts in series:
        edges, cumulative = cumulative_from_histogram(bins, counts)
        if not len(edges):
            continue
        axes[0].plot(edges[:-1], cumulative[:-1] / area_km2, color=color, ls=style, lw=width, label=label)
        centers = np.sqrt(bins[:-1] * bins[1:])
        mask = counts[:-1] > 0
        axes[1].plot(centers[mask], counts[:-1][mask], color=color, ls=style, lw=width)
        lo, hi = min(lo, edges[0]), max(hi, edges[-1])
    if exact is not None and len(exact):
        axes[0].plot(exact, np.arange(len(exact), 0, -1) / area_km2, color="#7ed0a8", ls=":", lw=1.0, label=exact_label)
    if np.isfinite(lo) and hi > 0:
        for ax in axes:
            ax.set_xlim(lo * 0.8, hi * 1.5)


def render_figure(run_id: str, name: str, scale: int = 1) -> Path:
    """One of the run's distribution figures as lines (scale 2 = the high-resolution version)."""
    name = FIGURE_ALIASES.get(name, name)
    if name not in FIGURE_NAMES:
        raise HTTPException(404, "unknown figure")
    target = run_dir(run_id) / "cache" / "{}_{}x.png".format(name, scale)
    if target.exists():
        return target
    out = output_of(run_id)
    summary = load_summary(run_id)
    p = summary["parameters"]
    area_km2 = p["regionWidth"] ** 2 / 1e6
    with PLOT_LOCK:
        fig = themed_figure(WIDE_FIG_SIZE, dpi=MAP_DPI * scale)
        if name == "turnover":
            axes = [fig.add_axes((0.075, 0.19, 0.66, 0.60))]
        else:
            axes = [fig.add_axes((0.075, 0.19, 0.395, 0.69)), fig.add_axes((0.585, 0.19, 0.395, 0.69))]
        for ax in axes:
            style_axes(ax)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(alpha=0.15, color=THEME["muted"], which="both")
        legend_kwargs = dict(fontsize=7.5, facecolor=THEME["figure"], edgecolor=THEME["grid"], labelcolor=THEME["ink"])

        if name == "impactors":
            bins, counts = out.histogram("impactors")
            plot_size_series(axes, [("impactors", "#79a6ff", "-", 1.6, bins, counts.astype(float))], area_km2)
            axes[0].set_xlabel("impactor diameter [m]"); axes[0].set_ylabel("N(> d) per km²"); axes[0].set_title("Cumulative impactor distribution", fontsize=10)
            axes[1].set_xlabel("impactor diameter [m]"); axes[1].set_ylabel("count per bin ({:.0f} per decade)".format(bins_per_decade(bins))); axes[1].set_title("Differential impactor distribution", fontsize=10)

        elif name == "craters":
            bins, counts = out.histogram("craters")
            visible = np.sort(np.asarray(out.craters()["diameter"], dtype=float))
            visible_counts, _ = np.histogram(visible, bins=bins)
            series = [("craters formed", THEME["accent"], "-", 2.0, bins, counts.astype(float)),
                      ("craters visible at the end", "#7ed0a8", "--", 1.4, bins, np.append(visible_counts, 0).astype(float))]
            plot_size_series(axes, series, area_km2, visible, "visible, exact (one point per crater)")
            axes[0].set_xlabel("crater diameter [m]"); axes[0].set_ylabel("N(> D) per km²"); axes[0].set_title("Cumulative crater distributions", fontsize=10)
            axes[1].set_xlabel("crater diameter [m]"); axes[1].set_ylabel("count per bin ({:.0f} per decade)".format(bins_per_decade(bins))); axes[1].set_title("Differential crater distributions", fontsize=10)
            legend = axes[0].legend(loc="lower left", **legend_kwargs); legend.get_frame().set_alpha(0.9)

        elif name == "rplot":
            # Relative size-frequency distribution (Arvidson et al. 1979): R = D_mean^3 dN/dD / A, in sqrt(2) bins.
            bins, counts = out.histogram("craters")
            factor = max(1, int(round(bins_per_decade(bins) * np.log10(np.sqrt(2)))))
            visible = np.asarray(out.craters()["diameter"], dtype=float)
            visible_counts, _ = np.histogram(visible, bins=bins)
            area_m2 = area_km2 * 1e6
            for label, color, style, width, c in [("craters formed", THEME["accent"], "-", 2.0, counts.astype(float)), ("craters visible at the end", "#7ed0a8", "--", 1.4, np.append(visible_counts, 0).astype(float))]:
                edges, merged = rebin(bins, c, factor)
                mean_d = np.sqrt(edges[:-1] * edges[1:])
                r = mean_d ** 3 * merged / (area_m2 * (edges[1:] - edges[:-1]))
                mask = merged > 0
                axes[0].plot(mean_d[mask], r[mask], color=color, ls=style, lw=width, marker="o", ms=3, label=label)
                if mask.any():
                    error = r[mask] / np.sqrt(merged[mask])
                    axes[0].errorbar(mean_d[mask], r[mask], yerr=error, fmt="none", ecolor=color, elinewidth=0.6, alpha=0.6, capsize=0)
                # the same craters as a differential per sqrt(2) bin, for reference
                axes[1].plot(mean_d[mask], merged[mask], color=color, ls=style, lw=width, marker="o", ms=3)
            from matplotlib.ticker import NullFormatter
            axes[0].yaxis.set_minor_formatter(NullFormatter())
            low, high = axes[0].get_ylim()
            axes[0].set_ylim(10 ** np.floor(np.log10(max(low, 1e-6))), 10 ** np.ceil(np.log10(max(high, 0.5))))
            # Geometric saturation (Gault 1970: N = 1.54 D^-2 per unit area, i.e. R = 3.08) and the empirical
            # saturation band at 10% and 1% of it (Hartmann; Richardson 2009).
            geometric = 3.08
            low, high = axes[0].get_ylim()
            axes[0].set_ylim(low, max(high, geometric * 2))
            axes[0].axhline(geometric, color=THEME["muted"], lw=0.9, ls="--")
            axes[0].axhline(0.1 * geometric, color=THEME["muted"], lw=0.8, ls=":"); axes[0].axhline(0.01 * geometric, color=THEME["muted"], lw=0.8, ls=":")
            # Each label sits just above its own line: the geometric one at the left end (the legend is at the
            # upper right), the band labels at the right end.
            x_left, x_right = axes[0].get_xlim()[0] * 1.08, axes[0].get_xlim()[1] / 1.08
            axes[0].text(x_left, geometric * 1.12, "geometric saturation", fontsize=6.5, color=THEME["muted"], ha="left", va="bottom")
            for level, text in [(0.1 * geometric, "10% of geometric (empirical saturation)"), (0.01 * geometric, "1% of geometric")]:
                axes[0].text(x_right, level * 1.12, text, fontsize=6.5, color=THEME["muted"], ha="right", va="bottom")
            axes[0].set_xlabel("crater diameter [m]"); axes[0].set_ylabel("R"); axes[0].set_title("R-plot (√2 bins)", fontsize=10)
            axes[1].set_xlabel("crater diameter [m]"); axes[1].set_ylabel("count per √2 bin"); axes[1].set_title("Counts per √2 bin", fontsize=10)
            legend = axes[0].legend(loc="upper right", **legend_kwargs); legend.get_frame().set_alpha(0.9)

        else:   # turnover: Gault et al. (1974), Eqs. 12-15 and Fig. 9
            # A point is turned over to depth d = D_a / 8 when a crater of diameter >= D_a forms within
            # c D_a / 2 of it (c = 0.75, craters of depth D_a / 4 excavating only their upper half):
            # lambda(d, t) = N(>= 8 d) pi (c 4 d)^2 t, with N the cumulative number of craters per unit
            # area and time. The flux is constant, so N at any age follows from the run's population by
            # scaling with t / T. The curves are the number of turnovers reached with 50% and 99%
            # probability, i.e. Poisson quantiles of lambda (Gault's Table 1 and Fig. 2).
            bins, counts = out.histogram("craters")
            end_time = float(p["endTime"])
            area_m2 = area_km2 * 1e6
            weights = counts[:-1].astype(float)
            cumulative = np.cumsum(weights[::-1])[::-1]   # craters with D >= bins[i]
            edges = bins[:-1]
            nonzero = weights > 0
            if nonzero.any():
                first, last = int(np.argmax(nonzero)), len(weights) - 1 - int(np.argmax(nonzero[::-1]))
                d_min, d_max = edges[first] / 8, bins[last + 1] / 8
                d = np.logspace(np.log10(d_min), np.log10(d_max), 300)
                keep = cumulative > 0
                n_final = np.exp(np.interp(np.log(8 * d), np.log(edges[keep]), np.log(cumulative[keep]))) / area_m2
                n_final[8 * d > bins[last + 1]] = 0.0
                ages_ma = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 4600.0]   # Gault's 10^4 ... 10^9 and 4.6 x 10^9 yr
                if not any(abs(np.log10(a / end_time)) < 0.15 for a in ages_ma):
                    ages_ma = sorted(ages_ma + [end_time])
                cmap = plt.get_cmap("viridis")
                n_max, lam_min = 1.0, np.inf
                for k, t in enumerate(ages_ma):
                    lam = n_final * t / end_time * np.pi * (0.75 * 4 * d) ** 2
                    shown = lam > 1e-3
                    if not shown.any():
                        continue
                    is_run = abs(t - end_time) < 1e-9 * end_time
                    color = cmap(0.25 + 0.75 * k / max(len(ages_ma) - 1, 1))
                    years = t * 1e6
                    exponent = int(np.floor(np.log10(years)))
                    mantissa = years / 10 ** exponent
                    label = ("{:g}×10$^{{{}}}$ yr".format(round(mantissa, 2), exponent) if abs(mantissa - 1) > 0.01 else "10$^{{{}}}$ yr".format(exponent)) + (" (run)" if is_run else "")
                    lw = 2.0 if is_run else 1.4
                    axes[0].plot(d[shown], lam[shown], color=color, lw=lw * 0.6, ls=":", alpha=0.9, label=label)
                    lam_min = min(lam_min, lam[shown].min())
                    n50, n99 = turnover_count(lam, 0.5), turnover_count(lam, 0.99)
                    if (n50 >= 1).any():
                        axes[0].plot(d[n50 >= 1], n50[n50 >= 1], color=color, lw=lw)
                        n_max = max(n_max, n50.max())
                    if (n99 >= 1).any():
                        axes[0].plot(d[n99 >= 1], n99[n99 >= 1], color=color, lw=lw * 0.75, ls="--")
                axes[0].plot([], [], color=THEME["ink"], lw=1.4, label="n at 50% probability"); axes[0].plot([], [], color=THEME["ink"], lw=1.0, ls="--", label="n at 99% probability")
                axes[0].plot([], [], color=THEME["ink"], lw=0.9, ls=":", label="expected number λ")
                axes[0].axhline(1, color=THEME["muted"], lw=0.7, ls="-", alpha=0.6)
                axes[0].set_xlim(d_min, d_max)
                low = 10 ** np.floor(np.log10(max(lam_min, 1e-3))) if np.isfinite(lam_min) else 1e-2
                axes[0].set_ylim(low, 10 ** np.ceil(np.log10(n_max * 1.2)))
                smallest = edges[first]
                axes[0].text(0.01, 0.02, "smallest crater in the run {:.3g} m, so no turnover below {:.3g} m is modelled".format(smallest, smallest / 8),
                             transform=axes[0].transAxes, fontsize=6.5, color=THEME["muted"])
                secondary = axes[0].secondary_xaxis("top", functions=(lambda x: 8 * x, lambda x: x / 8))
                secondary.set_xlabel("crater diameter $D_a$ = 8 × depth [m]", fontsize=8, color=THEME["muted"], labelpad=4)
                secondary.tick_params(colors=THEME["muted"], labelsize=7, which="both")
                for spine in secondary.spines.values():
                    spine.set_color(THEME["grid"])
            else:
                axes[0].text(0.5, 0.5, "no craters", transform=axes[0].transAxes, ha="center", color=THEME["muted"])
            axes[0].set_xlabel("depth of turnover [m]"); axes[0].set_ylabel("number of times turned over")
            axes[0].set_title("Regolith turnover, Gault et al. (1974) method, constant flux (compare their Fig. 9)", fontsize=10)
            legend = axes[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0, **legend_kwargs); legend.get_frame().set_alpha(0.9)
        atomic_savefig(fig, target)
        plt.close(fig)
    return target


def render_histograms(run_id: str) -> Path:   # kept for older links
    return render_figure(run_id, "craters", 1)


def render_animation(run_id: str, kind: str) -> Path:
    from PIL import Image

    summary = require_done(run_id)
    if kind not in MAP_KINDS:
        raise HTTPException(404, "unknown map kind")
    target = run_dir(run_id) / "cache" / "animation_{}.gif".format(kind)
    if target.exists():
        return target
    frames = [Image.open(render_map(run_id, kind, step)).convert("P", palette=Image.Palette.ADAPTIVE) for step in range(len(summary["steps"]))]
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    frames[0].save(str(temporary), format="GIF", save_all=True, append_images=frames[1:], duration=500, loop=0)
    os.replace(str(temporary), str(target))
    return target


def make_zip(run_id: str) -> Path:
    require_done(run_id)
    directory = run_dir(run_id)
    target = directory / "cache" / "regolit_{}.zip".format(run_id)
    if target.exists():
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    with zipfile.ZipFile(str(temporary), "w", zipfile.ZIP_DEFLATED) as archive:
        for sub in ("config", "output", "log"):
            for path in sorted((directory / sub).glob("*")):
                archive.write(str(path), "regolit_{}/{}/{}".format(run_id, sub, path.name))
        archive.write(str(directory / "summary.json"), "regolit_{}/summary.json".format(run_id))
    os.replace(str(temporary), str(target))
    return target


# ----------------------------------------------------------------------------------------------
# Routes
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# Saved settings (stored on the server, one JSON file per name)
# ----------------------------------------------------------------------------------------------
SETTINGS_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9-]{0,79}$")


def settings_id(name: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", name.strip().lower()).strip("-")[:80]
    if not slug:
        raise HTTPException(422, "the settings name needs at least one letter or digit")
    return slug


def settings_path(settings_id_: str) -> Path:
    if not SETTINGS_ID_PATTERN.match(settings_id_):
        raise HTTPException(404, "unknown settings")
    return SETTINGS_DIR / (settings_id_ + ".json")


def settings_entry(doc: Dict) -> Dict:
    presets = doc.get("presets") or {}
    return dict(id=doc["id"], name=doc["name"], saved=doc["saved"], body=presets.get("body"),
                production_function=presets.get("production_function"))


def list_settings() -> List[Dict]:
    SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
    out = []
    for path in SETTINGS_DIR.glob("*.json"):
        try:
            out.append(settings_entry(json.loads(path.read_text())))
        except (OSError, ValueError, KeyError):
            continue
    return sorted(out, key=lambda d: d["name"].lower())


def clean_settings(request: SettingsRequest) -> Dict:
    """Store the numeric parameters, layers and presets; nothing else from the client."""
    parameters: Dict[str, float] = {}
    for p in PARAMETERS:
        if p["name"] in request.parameters:
            try:
                parameters[p["name"]] = float(request.parameters[p["name"]])
            except (TypeError, ValueError):
                raise HTTPException(422, "{} is not a number".format(p["label"]))
    layers = None
    if request.layers:
        if len(request.layers) > MAX_LAYER_ROWS:
            raise HTTPException(422, "at most {} layer rows".format(MAX_LAYER_ROWS))
        layers = [[float(v) for v in row[:5]] for row in request.layers if len(row) == 5]
    return dict(id=settings_id(request.name), name=request.name.strip(), saved=time.strftime("%Y-%m-%dT%H:%M:%S"),
                parameters=parameters, layers=layers, presets=clean_presets(request.presets))


def save_settings(request: SettingsRequest) -> Dict:
    doc = clean_settings(request)
    path = settings_path(doc["id"])
    SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
    if not path.exists() and len(list(SETTINGS_DIR.glob("*.json"))) >= MAX_SAVED_SETTINGS:
        raise HTTPException(429, "too many saved settings ({}); delete some first".format(MAX_SAVED_SETTINGS))
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(doc, indent=1))
    os.replace(tmp, path)
    return doc


@app.get("/", response_class=HTMLResponse)
async def index() -> HTMLResponse:
    return HTMLResponse((STATIC_DIR / "index.html").read_text())


@app.get("/health")
async def health() -> Dict:
    return {"status": "ok", "binary": BINARY.exists(), "waiting": WAITING["count"]}


@app.get("/api/meta")
async def meta() -> Dict:
    return {
        "parameters": PARAMETERS,
        "defaults": default_parameters(),
        "layers": default_layers(),
        "bodies": scaling.BODIES,
        "body_fields": list(scaling.BODY_FIELDS),
        "production_functions": {k: {"label": v["label"], "available": v["available"], "reference": v["reference"]} for k, v in scaling.PRODUCTION_FUNCTIONS.items()},
        "default_presets": {"body": "moon", "production_function": "williams", "basement_auto": True},
        "map_kinds": {k: v[0] for k, v in MAP_KINDS.items()},
        "default_kind": "shaded_relief",
        "map_geometry": map_geometry(),
        "tests": TESTS,
        "limits": {"max_grid_cells": MAX_GRID_CELLS, "max_steps": MAX_STEPS, "max_output_mb": MAX_OUTPUT_BYTES // 2**20, "max_impacts": MAX_IMPACTS, "max_layer_rows": MAX_LAYER_ROWS},
    }


@app.post("/api/estimate")
async def post_estimate(request: RunRequest) -> Dict:
    """Fitted flux law, expected impacts, largest crater and suggested basement for a parameter set (no run)."""
    return estimate_for(request)


@app.get("/api/settings")
async def get_settings_list() -> List[Dict]:
    return list_settings()


@app.post("/api/settings")
async def create_settings(request: SettingsRequest) -> Dict:
    return save_settings(request)


@app.get("/api/settings/{settings_id_}")
async def get_settings(settings_id_: str) -> Dict:
    path = settings_path(settings_id_)
    if not path.exists():
        raise HTTPException(404, "unknown settings")
    return json.loads(path.read_text())


@app.post("/api/settings/{settings_id_}/delete")
async def delete_settings(settings_id_: str) -> Dict:
    path = settings_path(settings_id_)
    if not path.exists():
        raise HTTPException(404, "unknown settings")
    path.unlink()
    return {"deleted": settings_id_}


@app.get("/api/runs")
async def get_runs() -> List[Dict]:
    return list_runs()


BACKGROUND_TASKS: set = set()


async def run_in_background(run_id: str, overrides: Dict[str, float], layers: Optional[List[List[float]]], presets: Optional[Dict[str, object]]) -> None:
    """Wait for a slot, then run the model in a detached worker process (python -m web.worker).

    The worker lives in its own session, so it survives a reload or restart of this server, a
    dropped client connection and a proxy timeout; it writes summary.json itself when it is done."""
    WAITING["count"] += 1
    acquired = False
    try:
        async with RUN_SEMAPHORE:
            WAITING["count"] -= 1
            acquired = True
            if not (RUNS_DIR / run_id / "summary.json").exists():   # deleted while queued
                return
            try:
                returncode, stderr = await asyncio.to_thread(spawn_worker, run_id)
                if returncode != 0 and read_status(RUNS_DIR / run_id) in ("queued", "running"):
                    fail_run(run_id, "the worker exited with status {}: {}".format(returncode, stderr[-1500:]))
            except Exception as error:
                fail_run(run_id, "could not start the worker: {}".format(error))
    finally:
        if not acquired:
            WAITING["count"] -= 1
        await asyncio.to_thread(prune_runs)


def spawn_worker(run_id: str) -> Tuple[int, str]:
    """Run the detached worker to completion (blocking; called in a thread).

    Uses the standard library rather than the event loop's subprocess support on purpose: under
    uvloop the children inherit every inheritable descriptor, including uvicorn's listening socket,
    and a worker holding the socket keeps the port busy long after the server restarts."""
    import subprocess
    process = subprocess.Popen([sys.executable, "-m", "web.worker", run_id], cwd=str(REPO_ROOT), start_new_session=True,
                               close_fds=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    _, stderr = process.communicate()
    return process.returncode, stderr or ""


def start_background(run_id: str, overrides: Dict[str, float], layers: Optional[List[List[float]]], presets: Optional[Dict[str, object]]) -> None:
    task = asyncio.create_task(run_in_background(run_id, overrides, layers, presets))
    BACKGROUND_TASKS.add(task)
    task.add_done_callback(BACKGROUND_TASKS.discard)


async def watch_orphans() -> None:
    """Runs adopted from a previous server process: mark them failed if their worker vanished."""
    while True:
        await asyncio.sleep(15)
        try:
            for path in RUNS_DIR.iterdir():
                if path.is_dir() and RUN_ID_PATTERN.match(path.name) and read_status(path) == "running":
                    summary = json.loads((path / "summary.json").read_text())
                    if not pid_alive(summary.get("pid")):
                        fail_run(path.name, "the worker process disappeared while this run was in progress")
        except Exception:
            pass


@app.post("/api/runs", status_code=202)
async def create_run(request: RunRequest) -> Dict:
    """Validate, register the run and start it in the background; poll GET /api/runs/{id} for its status."""
    overrides, layers = validate(request)
    if WAITING["count"] >= MAX_QUEUE:
        raise HTTPException(429, "the server is busy; please try again in a minute")
    run_id = new_run_id()
    summary = queue_run(run_id, overrides, layers, request.presets)
    start_background(run_id, overrides, layers, request.presets)
    summary["waiting"] = WAITING["count"]
    return summary


@app.get("/api/runs/{run_id}")
async def get_run(run_id: str) -> Dict:
    summary = load_summary(run_id)
    if summary.get("status", "done") in ("queued", "running"):
        summary.update(progress_of(run_id, summary))
        summary["waiting"] = WAITING["count"]
    return summary


@app.post("/api/runs/{run_id}/delete")
async def delete_run(run_id: str) -> Dict:
    if load_summary(run_id).get("status", "done") == "running":
        raise HTTPException(409, "this run is still running; wait for it to finish before deleting it")
    shutil.rmtree(run_dir(run_id), ignore_errors=True)
    with PLOT_LOCK:
        SUBSURFACE_CACHE.pop(run_id, None)
    return {"deleted": run_id}


@app.get("/api/runs/{run_id}/log", response_class=PlainTextResponse)
async def get_log(run_id: str, tail: Optional[int] = Query(None, ge=1, le=64 * 1024 * 1024)) -> str:
    """The model log; ?tail=N returns only its last N bytes."""
    path = run_dir(run_id) / "log" / "log.txt"
    if not path.exists():
        return ""
    if tail:
        return await asyncio.to_thread(read_tail, path, tail)
    return FileResponse(str(path), media_type="text/plain")


@app.get("/api/runs/{run_id}/map/{kind}/{step}.png")
async def get_map(run_id: str, kind: str, step: int,
                  az: float = Query(315.0, ge=0, le=360), alt: float = Query(25.0, ge=1, le=89)) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(not_found_if_deleted, render_map, run_id, kind, step, az, alt)), media_type="image/png")


@app.get("/api/runs/{run_id}/section.png")
async def get_section(run_id: str, axis: str = Query("x", pattern="^[xy]$"), at: float = 0.0,
                      depth: float = Query(2.0, gt=0, le=1000), step: int = -1) -> Response:
    summary = load_summary(run_id)
    if step < 0:
        step = len(summary["steps"]) + step
    png = await asyncio.to_thread(not_found_if_deleted, render_section, run_id, axis, at, depth, step)
    return Response(content=png, media_type="image/png")


@app.get("/api/runs/{run_id}/histograms.png")
async def get_histograms(run_id: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(not_found_if_deleted, render_histograms, run_id)), media_type="image/png")


@app.get("/api/runs/{run_id}/figure/{name}.png")
async def get_figure(run_id: str, name: str, scale: int = Query(1, ge=1, le=3)) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(not_found_if_deleted, render_figure, run_id, name, scale)), media_type="image/png")


@app.get("/api/runs/{run_id}/animation/{kind}.gif")
async def get_animation(run_id: str, kind: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(not_found_if_deleted, render_animation, run_id, kind)), media_type="image/gif")


@app.get("/api/runs/{run_id}/craters.csv")
async def get_craters(run_id: str) -> FileResponse:
    return FileResponse(str(run_dir(run_id) / "output" / "existing_craters.txt"), media_type="text/csv", filename="craters_{}.csv".format(run_id))


@app.get("/api/runs/{run_id}/download.zip")
async def get_download(run_id: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(not_found_if_deleted, make_zip, run_id)), media_type="application/zip", filename="regolit_{}.zip".format(run_id))
