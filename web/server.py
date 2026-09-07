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
REGOLIT_WEB_MAX_RUNS      runs kept on disk, oldest deleted first (default 200)
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
from regolit.io import RegolitOutput, Subsurface, read_config, read_layers  # noqa: E402

# ----------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------
STATIC_DIR = Path(__file__).resolve().parent / "static"
BINARY = Path(os.environ.get("REGOLIT_BINARY", REPO_ROOT / "build" / "apps" / "regolit_main.run"))
RUNS_DIR = Path(os.environ.get("REGOLIT_WEB_RUNS", REPO_ROOT / "runs" / "web"))
CONFIG_TEMPLATE = REPO_ROOT / "config" / "config.cfg"
LAYERS_TEMPLATE = REPO_ROOT / "config" / "layers.cfg"
MAX_CONCURRENT = int(os.environ.get("REGOLIT_WEB_CONCURRENCY", "2"))
MAX_QUEUE = int(os.environ.get("REGOLIT_WEB_QUEUE", "8"))
RUN_TIMEOUT = float(os.environ.get("REGOLIT_WEB_TIMEOUT", "120"))
MAX_RUNS_KEPT = int(os.environ.get("REGOLIT_WEB_MAX_RUNS", "200"))
AUTH_USER = os.environ.get("REGOLIT_WEB_USER", "")
AUTH_PASSWORD = os.environ.get("REGOLIT_WEB_PASSWORD", "")

# Hard limits protecting the server. Every saved step writes seven maps of the output grid, so the
# number of steps is bounded through the total output size rather than by a fixed count.
MAX_GRID_CELLS = 500 * 500
MAX_STEPS = 500
MAX_OUTPUT_BYTES = 800 * 1024 * 1024
MAPS_PER_STEP = 7
MAX_IMPACTS = 200_000
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
    dict(group="Domain", name="regionWidth", label="Region width", unit="m", min=50, max=5000, kind="number",
         description="Side of the square, periodic domain."),
    dict(group="Domain", name="resolution", label="Resolution", unit="m/pixel", min=0.5, max=50, kind="number",
         description="Cell size. Cells = (width / resolution)^2, at most 500 x 500."),
    dict(group="Domain", name="downsamplingResolution", label="Map output resolution", unit="m/pixel", min=0.5, max=200, kind="number",
         description="Maps are averaged to this resolution before saving (>= resolution). Cross-sections use the full resolution."),
    dict(group="Time", name="endTime", label="Duration", unit="Ma", min=0.1, max=4500, kind="number",
         description="Simulated time."),
    dict(group="Time", name="printTimeStep", label="Output interval", unit="Ma", min=0.01, max=4500, kind="number",
         description="Time between saved steps; up to 500 steps, bounded by the total output size."),
    dict(group="Time", name="randomSeed", label="Random seed", unit="", min=0, max=2**31 - 1, kind="int",
         description="Seed of the impactor sequence."),
    dict(group="Impactors", name="minimumImpactorDiameter", label="Minimum impactor diameter", unit="m", min=0.02, max=50, kind="number",
         description="Smallest impactor drawn from the power-law distribution."),
    dict(group="Impactors", name="slope_b", label="Size-distribution slope b", unit="", min=1.2, max=4.5, kind="number",
         description="Cumulative slope: N(>D) ~ D^-b."),
    dict(group="Impactors", name="fluxConstant_c", label="Flux constant c", unit="m^-2 Ma^-1", min=1e-12, max=1e-3, kind="number",
         description="N(>1 m) per m^2 per Ma at Earth."),
    dict(group="Impactors", name="earthFluxRatioCoefficient", label="Body/Earth flux ratio", unit="", min=0.01, max=50, kind="number",
         description="Scales the flux to the target body."),
    dict(group="Impactors", name="impactorDensity", label="Impactor density", unit="kg/m^3", min=300, max=9000, kind="number"),
    dict(group="Impactors", name="meanImpactVelocity", label="Impact velocity", unit="m/s", min=500, max=80000, kind="number"),
    dict(group="Impactors", name="isEmplaceSecondaries", label="Secondary craters", unit="", min=0, max=1, kind="bool",
         description="Form secondary craters around each primary."),
    dict(group="Impactors", name="slope_secondaries", label="Secondaries slope", unit="", min=1.5, max=8, kind="number"),
    dict(group="Target", name="g", label="Gravity", unit="m/s^2", min=0.01, max=30, kind="number"),
    dict(group="Target", name="targetDensity", label="Target density", unit="kg/m^3", min=300, max=6000, kind="number"),
    dict(group="Target", name="k1", label="Scaling constant k1", unit="", min=0.01, max=2, kind="number",
         description="Holsapple (1993) crater-volume scaling constant."),
    dict(group="Target", name="mu", label="Scaling exponent mu", unit="", min=0.3, max=0.7, kind="number"),
    dict(group="Target", name="Ybar", label="Effective strength", unit="Pa", min=0, max=1e9, kind="number"),
    dict(group="Target", name="angleOfRepose", label="Angle of repose", unit="deg", min=5, max=80, kind="number",
         description="Slopes steeper than this fail at every output step."),
    dict(group="Craters", name="craterProfileType", label="Cavity shape", unit="", min=1, max=2, kind="choice",
         choices={"1": "parabolic", "2": "bowl (spherical cap)"}),
    dict(group="Craters", name="depthToDiameter", label="Depth / diameter", unit="", min=0.02, max=0.5, kind="number"),
    dict(group="Craters", name="rimToDiameter", label="Rim height / diameter", unit="", min=0, max=0.2, kind="number"),
    dict(group="Craters", name="rimDropoffExponent", label="Rim dropoff exponent", unit="", min=1, max=8, kind="number"),
    dict(group="Craters", name="isEmplaceEjecta", label="Ejecta blanket", unit="", min=0, max=1, kind="bool"),
    dict(group="Craters", name="ejectaSpread", label="Ejecta extent", unit="radii", min=2, max=16, kind="int"),
    dict(group="Craters", name="numberOfZModelShells", label="Z-model shells", unit="", min=10, max=500, kind="int"),
    dict(group="Craters", name="ejectaVolatileRetention", label="Ice retained in ejecta", unit="fraction", min=0, max=1, kind="number"),
    dict(group="Craters", name="ejectaSootRetention", label="Soot retained in ejecta", unit="fraction", min=0, max=1, kind="number"),
    dict(group="Craters", name="minimumLayerThickness", label="Minimum layer thickness", unit="m", min=0, max=1, kind="number",
         description="Thinner deposits are mixed into the surface layer."),
    dict(group="Subsurface", name="initialThickness", label="Basement thickness", unit="m", min=1, max=10000, kind="number"),
    dict(group="Subsurface", name="depthToIntegrate", label="Integration depth", unit="m", min=0.005, max=100, kind="number",
         description="Depth of the integrated-composition maps."),
    dict(group="Subsurface", name="iceEmplacementInterval", label="Ice deposition interval", unit="Ma", min=0.01, max=4500, kind="number"),
    dict(group="Subsurface", name="iceEmplacementThickness", label="Ice deposition thickness", unit="m", min=0, max=10, kind="number",
         description="Ice added on icy surfaces at every interval (0 disables)."),
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
    prune_runs()
    if not BINARY.exists():
        print("WARNING: model binary not found at {}; run `make` first.".format(BINARY), file=sys.stderr)
    yield


app = FastAPI(title="REGOLIT", docs_url=None, redoc_url=None, dependencies=[Depends(require_auth)], lifespan=lifespan)
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

RUN_SEMAPHORE = asyncio.Semaphore(MAX_CONCURRENT)
WAITING = {"count": 0}
PLOT_LOCK = threading.Lock()
SUBSURFACE_CACHE: "OrderedDict[str, Subsurface]" = OrderedDict()
SUBSURFACE_CACHE_SIZE = 4


class RunRequest(BaseModel):
    parameters: Dict[str, float] = Field(default_factory=dict)
    layers: Optional[List[List[float]]] = None   # rows of (class, thickness, regolith, ice, soot), bottom-up


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
    cells = round(width / res) ** 2
    if cells > MAX_GRID_CELLS:
        raise HTTPException(400, "the grid would have {:,} cells; the limit is 500 x 500. Increase the resolution or shrink the region.".format(int(cells)))
    if effective["downsamplingResolution"] < res:
        overrides["downsamplingResolution"] = res
        effective["downsamplingResolution"] = res
    steps = math.ceil(effective["endTime"] / effective["printTimeStep"]) + 1
    if steps > MAX_STEPS:
        raise HTTPException(400, "{} output steps requested; the limit is {}. Increase the output interval.".format(steps, MAX_STEPS))
    output_cells = round(width / effective["downsamplingResolution"]) ** 2
    output_bytes = steps * output_cells * MAPS_PER_STEP * 8
    if output_bytes > MAX_OUTPUT_BYTES:
        raise HTTPException(400, "the run would write {:.0f} MB of maps ({} steps of {:,} cells); the limit is {:.0f} MB. Increase the output interval or the map output resolution.".format(
            output_bytes / 2**20, steps, int(output_cells), MAX_OUTPUT_BYTES / 2**20))
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


def prune_runs() -> None:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    runs = sorted([p for p in RUNS_DIR.iterdir() if p.is_dir()], key=lambda p: p.stat().st_mtime)
    for old in runs[: max(0, len(runs) - MAX_RUNS_KEPT)]:
        shutil.rmtree(old, ignore_errors=True)


def execute(run_id: str, overrides: Dict[str, float], layers: Optional[List[List[float]]]) -> Dict:
    """Run the model (blocking), compact the layer stacks and write summary.json."""
    workdir = RUNS_DIR / run_id
    started = time.time()
    out = regolit.run(overrides, workdir, binary=BINARY, layers=layers if layers is not None else LAYERS_TEMPLATE,
                      build_if_missing=False, quiet=True, timeout=RUN_TIMEOUT)
    series = out.elevation_series()
    match = re.search(r"Number of craters in simulation: (\d+)", out.log())

    # Store the final layer stacks compactly (the raw file is large) and drop the raw file.
    has_subsurface = False
    for suffix in out.subsurface_steps():
        raw = out.output_dir / "subsurface_{}.out".format(suffix)
        if raw.exists():
            Subsurface.load  # noqa: B018 (keep the name imported for readers)
            sub = regolit.read_subsurface(raw)
            sub.save(out.output_dir / "subsurface_{}.npz".format(suffix))
            raw.unlink()
        has_subsurface = True

    x_full, _ = out.full_resolution_coordinates()
    summary = {
        "id": run_id,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "duration_s": round(time.time() - started, 2),
        "parameters": {p["name"]: out.config.get(p["name"]) for p in PARAMETERS},
        "layers": layers if layers is not None else default_layers(),
        "steps": out.steps,
        "times": [float(t) for t in out.times],
        "grid": out.n,
        "resolution": out.resolution,
        "full_grid": int(len(x_full)),
        "full_resolution": float(out.config.get("resolution", out.resolution)),
        "extent": out.extent,
        "elevation_min": float(series.min()),
        "elevation_max": float(series.max()),
        "total_craters": int(match.group(1)) if match else None,
        "visible_craters": int(len(out.craters()["x"])),
        "has_subsurface": has_subsurface,
    }
    (workdir / "summary.json").write_text(json.dumps(summary, indent=1))
    (workdir / "cache").mkdir(exist_ok=True)
    return summary


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
        item = {key: summary.get(key) for key in ("id", "created", "duration_s", "grid", "resolution", "total_craters", "visible_craters")}
        item["parameters"] = {key: summary.get("parameters", {}).get(key) for key in ("regionWidth", "resolution", "endTime", "randomSeed")}
        runs.append(item)
    runs.sort(key=lambda item: item["id"], reverse=True)
    return runs


def load_summary(run_id: str) -> Dict:
    return json.loads((run_dir(run_id) / "summary.json").read_text())


def output_of(run_id: str) -> RegolitOutput:
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
    summary = load_summary(run_id)
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
    summary = load_summary(run_id)
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


def render_histograms(run_id: str) -> Path:
    target = run_dir(run_id) / "cache" / "histograms.png"
    if target.exists():
        return target
    out = output_of(run_id)
    with PLOT_LOCK:
        fig = themed_figure(WIDE_FIG_SIZE)
        axes = [fig.add_axes((0.07, 0.19, 0.40, 0.69)), fig.add_axes((0.585, 0.19, 0.40, 0.69))]
        for ax in axes:
            style_axes(ax)
        styles = [("craters", THEME["accent"], "-", 2.2, "craters formed"), ("existing_craters", "#7ed0a8", "--", 1.4, "craters visible at the end"), ("impactors", "#79a6ff", "-", 1.4, "impactors")]
        upper = 1.0
        for name, color, style, width, label in styles:
            bins, counts = out.histogram(name)
            axes[0].step(bins, np.maximum(counts, 0.5), where="post", color=color, ls=style, lw=width, label=label)
            if counts.any():
                top = np.nonzero(counts)[0].max() + 1
                upper = max(upper, bins[top] if top < len(bins) else bins[-1])
        axes[0].set_xlim(bins[0] * 0.8, upper * 3)
        axes[0].set_xscale("log")
        axes[0].set_yscale("log")
        axes[0].set_xlabel("diameter [m]")
        axes[0].set_ylabel("count per bin")
        axes[0].set_title("Size distributions", fontsize=10)
        legend = axes[0].legend(fontsize=8, facecolor=THEME["figure"], edgecolor=THEME["grid"], labelcolor=THEME["ink"])
        legend.get_frame().set_alpha(0.9)
        bins, counts = out.histogram("depth")
        axes[1].step(bins, np.maximum(counts, 0.5), where="post", color="#f08a7a")
        if counts.any():
            top = np.nonzero(counts)[0].max() + 1
            axes[1].set_xlim(bins[0] * 0.8, (bins[top] if top < len(bins) else bins[-1]) * 3)
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
        axes[1].set_xlabel("current depth [m]")
        axes[1].set_ylabel("count per bin")
        axes[1].set_title("Depths of visible craters", fontsize=10)
        for ax in axes:
            ax.grid(alpha=0.15, color=THEME["muted"], which="both")
        atomic_savefig(fig, target)
        plt.close(fig)
    return target


def render_animation(run_id: str, kind: str) -> Path:
    from PIL import Image

    summary = load_summary(run_id)
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
        "map_kinds": {k: v[0] for k, v in MAP_KINDS.items()},
        "default_kind": "shaded_relief",
        "map_geometry": map_geometry(),
        "limits": {"max_grid_cells": MAX_GRID_CELLS, "max_steps": MAX_STEPS, "max_output_mb": MAX_OUTPUT_BYTES // 2**20, "max_impacts": MAX_IMPACTS, "max_layer_rows": MAX_LAYER_ROWS},
    }


@app.get("/api/runs")
async def get_runs() -> List[Dict]:
    return list_runs()


@app.post("/api/runs")
async def create_run(request: RunRequest) -> Dict:
    overrides, layers = validate(request)
    if WAITING["count"] >= MAX_QUEUE:
        raise HTTPException(429, "the server is busy; please try again in a minute")
    run_id = new_run_id()
    WAITING["count"] += 1
    acquired = False
    try:
        async with RUN_SEMAPHORE:
            WAITING["count"] -= 1
            acquired = True
            try:
                summary = await asyncio.to_thread(execute, run_id, overrides, layers)
            except Exception as error:  # model failure: report the message, drop the directory
                shutil.rmtree(RUNS_DIR / run_id, ignore_errors=True)
                raise HTTPException(500, "the model run failed: {}".format(str(error)[:800]))
    finally:
        if not acquired:
            WAITING["count"] -= 1
    prune_runs()
    return summary


@app.get("/api/runs/{run_id}")
async def get_run(run_id: str) -> Dict:
    return load_summary(run_id)


@app.post("/api/runs/{run_id}/delete")
async def delete_run(run_id: str) -> Dict:
    shutil.rmtree(run_dir(run_id), ignore_errors=True)
    with PLOT_LOCK:
        SUBSURFACE_CACHE.pop(run_id, None)
    return {"deleted": run_id}


@app.get("/api/runs/{run_id}/log", response_class=PlainTextResponse)
async def get_log(run_id: str) -> str:
    path = run_dir(run_id) / "log" / "log.txt"
    return path.read_text() if path.exists() else ""


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


@app.get("/api/runs/{run_id}/animation/{kind}.gif")
async def get_animation(run_id: str, kind: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(not_found_if_deleted, render_animation, run_id, kind)), media_type="image/gif")


@app.get("/api/runs/{run_id}/craters.csv")
async def get_craters(run_id: str) -> FileResponse:
    return FileResponse(str(run_dir(run_id) / "output" / "existing_craters.txt"), media_type="text/csv", filename="craters_{}.csv".format(run_id))


@app.get("/api/runs/{run_id}/download.zip")
async def get_download(run_id: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(not_found_if_deleted, make_zip, run_id)), media_type="application/zip", filename="regolit_{}.zip".format(run_id))
