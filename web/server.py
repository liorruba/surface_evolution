"""REGOLIT web UI.

A small FastAPI service: the browser posts a parameter set, the server runs the C++ model in a
run directory, and the page shows maps, cross-sections and histograms rendered on the server.
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
import io
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
from pathlib import Path
from contextlib import asynccontextmanager
from typing import Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from fastapi import Depends, FastAPI, HTTPException, Query  # noqa: E402
from fastapi.responses import FileResponse, HTMLResponse, PlainTextResponse  # noqa: E402
from fastapi.security import HTTPBasic, HTTPBasicCredentials  # noqa: E402
from fastapi.staticfiles import StaticFiles  # noqa: E402
from pydantic import BaseModel, Field  # noqa: E402

import regolit  # noqa: E402
from regolit.io import RegolitOutput, read_config, read_layers  # noqa: E402

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

# Hard limits protecting the server.
MAX_GRID_CELLS = 500 * 500
MAX_OUTPUT_CELLS = 500 * 500
MAX_STEPS = 50
MAX_IMPACTS = 200_000
MAX_LAYER_ROWS = 40

RUN_ID_PATTERN = re.compile(r"^[0-9]{8}-[0-9]{6}-[0-9a-f]{6}$")

# Parameters the UI exposes: name, label, unit, min, max, kind, description. Values not listed
# here stay at the repository defaults (config/config.cfg).
PARAMETERS: List[Dict] = [
    # Domain
    dict(group="Domain", name="regionWidth", label="Region width", unit="m", min=50, max=5000, kind="number",
         description="Side of the square, periodic domain."),
    dict(group="Domain", name="resolution", label="Resolution", unit="m/pixel", min=0.5, max=50, kind="number",
         description="Cell size. Cells = (width / resolution)^2, at most 500 x 500."),
    dict(group="Domain", name="downsamplingResolution", label="Output resolution", unit="m/pixel", min=0.5, max=200, kind="number",
         description="Maps are averaged to this resolution before saving (>= resolution)."),
    # Time
    dict(group="Time", name="endTime", label="Duration", unit="Ma", min=0.1, max=4500, kind="number",
         description="Simulated time."),
    dict(group="Time", name="printTimeStep", label="Output interval", unit="Ma", min=0.01, max=4500, kind="number",
         description="Time between saved steps; at most 50 steps per run."),
    dict(group="Time", name="randomSeed", label="Random seed", unit="", min=0, max=2**31 - 1, kind="int",
         description="Seed of the impactor sequence."),
    # Impactors
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
    # Target
    dict(group="Target", name="g", label="Gravity", unit="m/s^2", min=0.01, max=30, kind="number"),
    dict(group="Target", name="targetDensity", label="Target density", unit="kg/m^3", min=300, max=6000, kind="number"),
    dict(group="Target", name="k1", label="Scaling constant k1", unit="", min=0.01, max=2, kind="number",
         description="Holsapple (1993) crater-volume scaling constant."),
    dict(group="Target", name="mu", label="Scaling exponent mu", unit="", min=0.3, max=0.7, kind="number"),
    dict(group="Target", name="Ybar", label="Effective strength", unit="Pa", min=0, max=1e9, kind="number"),
    dict(group="Target", name="angleOfRepose", label="Angle of repose", unit="deg", min=5, max=80, kind="number",
         description="Slopes steeper than this fail at every output step."),
    # Craters
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
    # Subsurface
    dict(group="Subsurface", name="initialThickness", label="Basement thickness", unit="m", min=1, max=10000, kind="number"),
    dict(group="Subsurface", name="depthToIntegrate", label="Integration depth", unit="m", min=0.005, max=100, kind="number",
         description="Depth of the integrated-composition maps."),
    dict(group="Subsurface", name="iceEmplacementInterval", label="Ice deposition interval", unit="Ma", min=0.01, max=4500, kind="number"),
    dict(group="Subsurface", name="iceEmplacementThickness", label="Ice deposition thickness", unit="m", min=0, max=10, kind="number",
         description="Ice added on icy surfaces at every interval (0 disables)."),
]
PARAMETER_INDEX = {p["name"]: p for p in PARAMETERS}
MAP_KINDS = {
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


class RunRequest(BaseModel):
    parameters: Dict[str, float] = Field(default_factory=dict)
    layers: Optional[List[List[float]]] = None   # rows of (class, thickness, regolith, ice, soot), bottom-up


# ----------------------------------------------------------------------------------------------
# Helpers
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
    steps = math.ceil(effective["endTime"] / effective["printTimeStep"])
    if steps > MAX_STEPS:
        raise HTTPException(400, "{} output steps requested; the limit is {}. Increase the output interval.".format(steps, MAX_STEPS))
    impacts = effective["fluxConstant_c"] * effective["minimumImpactorDiameter"] ** (-effective["slope_b"]) * width ** 2 \
        * effective["endTime"] * effective["earthFluxRatioCoefficient"]
    if impacts > MAX_IMPACTS:
        raise HTTPException(400, "about {:,} impacts would be simulated; the limit is {:,}. Shorten the run, shrink the region or raise the minimum impactor diameter.".format(int(impacts), MAX_IMPACTS))
    overrides["isPrintSubsurface"] = 0   # full layer stacks are too large to store for web runs

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
    """Run the model (blocking) and write summary.json."""
    workdir = RUNS_DIR / run_id
    started = time.time()
    out = regolit.run(overrides, workdir, binary=BINARY, layers=layers if layers is not None else LAYERS_TEMPLATE,
                      build_if_missing=False, quiet=True, timeout=RUN_TIMEOUT)
    series = out.elevation_series()
    match = re.search(r"Number of craters in simulation: (\d+)", out.log())
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
        "extent": out.extent,
        "elevation_min": float(series.min()),
        "elevation_max": float(series.max()),
        "total_craters": int(match.group(1)) if match else None,
        "visible_craters": int(len(out.craters()["x"])),
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


def step_index(summary: Dict, step: int) -> int:
    if not 0 <= step < len(summary["steps"]):
        raise HTTPException(404, "step out of range")
    return step


def render_map(run_id: str, kind: str, step: int) -> Path:
    summary = load_summary(run_id)
    step = step_index(summary, step)
    if kind not in MAP_KINDS:
        raise HTTPException(404, "unknown map kind")
    target = run_dir(run_id) / "cache" / "map_{}_{}.png".format(kind, step)
    if target.exists():
        return target
    out = output_of(run_id)
    label, unit, cmap = MAP_KINDS[kind]
    if kind == "elevation":
        data, vmin, vmax = out.elevation(step), summary["elevation_min"], summary["elevation_max"]
    else:
        where, species = kind.split("_")
        data = out.surface_fraction(species, step) if where == "surface" else out.integrated_fraction(species, step)
        vmin, vmax = 0.0, 1.0
        if where == "integrated":
            label = "{} (top {:g} m)".format(label, summary["parameters"].get("depthToIntegrate", float("nan")))
    with PLOT_LOCK:
        fig, ax = plt.subplots(figsize=(6.6, 5.6), dpi=110)
        image = ax.imshow(data, origin="lower", extent=summary["extent"], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title("{} at t = {:g} Ma".format(label, summary["times"][step]))
        fig.colorbar(image, ax=ax, fraction=0.046, label=unit)
        fig.tight_layout()
        fig.savefig(str(target))
        plt.close(fig)
    return target


def render_profile(run_id: str, step: int) -> Path:
    summary = load_summary(run_id)
    step = step_index(summary, step)
    target = run_dir(run_id) / "cache" / "profile_{}.png".format(step)
    if target.exists():
        return target
    out = output_of(run_id)
    final = out.elevation(-1)
    row = int(np.unravel_index(np.argmin(final), final.shape)[0])   # through the lowest point of the final surface
    z = out.elevation(step)
    with PLOT_LOCK:
        fig, ax = plt.subplots(figsize=(6.6, 3.2), dpi=110)
        ax.plot(out.x, z[row, :], color="C0", lw=1.2)
        ax.set_ylim(summary["elevation_min"], summary["elevation_max"])
        ax.set_xlabel("x [m]")
        ax.set_ylabel("elevation [m]")
        ax.set_title("Cross-section at y = {:.1f} m, t = {:g} Ma".format(out.y[row], summary["times"][step]))
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(str(target))
        plt.close(fig)
    return target


def render_histograms(run_id: str) -> Path:
    target = run_dir(run_id) / "cache" / "histograms.png"
    if target.exists():
        return target
    out = output_of(run_id)
    with PLOT_LOCK:
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), dpi=110)
        styles = [("craters", "C0", "-", 2.2, "craters formed"), ("existing_craters", "C2", "--", 1.4, "craters visible at the end"), ("impactors", "C1", "-", 1.4, "impactors")]
        upper = 1.0
        for name, color, style, width, label in styles:
            bins, counts = out.histogram(name)
            axes[0].step(bins, np.maximum(counts, 0.5), where="post", color=color, ls=style, lw=width, label=label)
            if counts.any():
                upper = max(upper, bins[np.nonzero(counts)[0].max() + 1] if np.nonzero(counts)[0].max() + 1 < len(bins) else bins[-1])
        axes[0].set_xlim(bins[0] * 0.8, upper * 3)
        axes[0].set_xscale("log")
        axes[0].set_yscale("log")
        axes[0].set_xlabel("diameter [m]")
        axes[0].set_ylabel("count per bin")
        axes[0].set_title("Size distributions")
        axes[0].legend(fontsize=8)
        bins, counts = out.histogram("depth")
        axes[1].step(bins, np.maximum(counts, 0.5), where="post", color="C3")
        if counts.any():
            top = np.nonzero(counts)[0].max() + 1
            axes[1].set_xlim(bins[0] * 0.8, (bins[top] if top < len(bins) else bins[-1]) * 3)
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
        axes[1].set_xlabel("current depth [m]")
        axes[1].set_ylabel("count per bin")
        axes[1].set_title("Depths of visible craters")
        for ax in axes:
            ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(str(target))
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
    frames[0].save(str(target), save_all=True, append_images=frames[1:], duration=500, loop=0)
    return target


def make_zip(run_id: str) -> Path:
    directory = run_dir(run_id)
    target = directory / "cache" / "regolit_{}.zip".format(run_id)
    if target.exists():
        return target
    with zipfile.ZipFile(str(target), "w", zipfile.ZIP_DEFLATED) as archive:
        for sub in ("config", "output", "log"):
            for path in sorted((directory / sub).glob("*")):
                archive.write(str(path), "regolit_{}/{}/{}".format(run_id, sub, path.name))
        archive.write(str(directory / "summary.json"), "regolit_{}/summary.json".format(run_id))
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
        "limits": {"max_grid_cells": MAX_GRID_CELLS, "max_steps": MAX_STEPS, "max_impacts": MAX_IMPACTS, "max_layer_rows": MAX_LAYER_ROWS},
    }


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


@app.get("/api/runs")
async def get_runs() -> List[Dict]:
    return list_runs()


@app.get("/api/runs/{run_id}")
async def get_run(run_id: str) -> Dict:
    return load_summary(run_id)


@app.post("/api/runs/{run_id}/delete")
async def delete_run(run_id: str) -> Dict:
    shutil.rmtree(run_dir(run_id), ignore_errors=True)
    return {"deleted": run_id}


@app.get("/api/runs/{run_id}/log", response_class=PlainTextResponse)
async def get_log(run_id: str) -> str:
    path = run_dir(run_id) / "log" / "log.txt"
    return path.read_text() if path.exists() else ""


@app.get("/api/runs/{run_id}/map/{kind}/{step}.png")
async def get_map(run_id: str, kind: str, step: int) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(render_map, run_id, kind, step)), media_type="image/png")


@app.get("/api/runs/{run_id}/profile/{step}.png")
async def get_profile(run_id: str, step: int) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(render_profile, run_id, step)), media_type="image/png")


@app.get("/api/runs/{run_id}/histograms.png")
async def get_histograms(run_id: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(render_histograms, run_id)), media_type="image/png")


@app.get("/api/runs/{run_id}/animation/{kind}.gif")
async def get_animation(run_id: str, kind: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(render_animation, run_id, kind)), media_type="image/gif")


@app.get("/api/runs/{run_id}/craters.csv")
async def get_craters(run_id: str) -> FileResponse:
    return FileResponse(str(run_dir(run_id) / "output" / "existing_craters.txt"), media_type="text/csv", filename="craters_{}.csv".format(run_id))


@app.get("/api/runs/{run_id}/download.zip")
async def get_download(run_id: str) -> FileResponse:
    return FileResponse(str(await asyncio.to_thread(make_zip, run_id)), media_type="application/zip", filename="regolit_{}.zip".format(run_id))
