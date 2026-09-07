"""Readers (and config writers) for REGOLIT.

Conventions
-----------
* Every 2-D map is returned indexed ``[y, x]`` (row = y index, column = x index). This matches the
  MATLAB readers in ``vis/`` and ``matplotlib.pyplot.imshow(z, origin="lower", extent=[x0, x1, y0, y1])``.
* The C++ code writes maps x-major (all y for the first x, then the next x, ...), so a file is
  reshaped to ``(n, n)`` and transposed.
* Surface and integrated-subsurface maps are written at ``downsamplingResolution``; the full layer
  stacks (``subsurface_XX.out``) are written at the grid resolution.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

SPECIES = ("regolith", "ice", "soot")
SURFACE_FILES = {"regolith": "regolithFraction", "ice": "iceFraction", "soot": "sootFraction"}
INTEGRATED_FILES = {"regolith": "depthRegolithFraction", "ice": "depthIceFraction", "soot": "depthSootFraction"}
HISTOGRAM_FILES = {
    "craters": "craters_histogram.txt",
    "impactors": "impactor_histogram.txt",
    "depth": "depth_histogram.txt",
    "existing_craters": "existing_craters_histogram.txt",
}

PathLike = Union[str, Path]


# ----------------------------------------------------------------------------------------------
# Configuration files
# ----------------------------------------------------------------------------------------------
def read_config(path: PathLike) -> Dict[str, float]:
    """Parse ``config.cfg`` (``key value`` lines, ``//`` comments) into a dict of floats."""
    params: Dict[str, float] = {}
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            params[parts[0]] = float(parts[1])
        except ValueError:
            continue
    return params


def _format_value(value) -> str:
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        text = "{:.12g}".format(float(value))
        return text
    return str(value)


def write_config(path: PathLike, params: Dict[str, object], template: Optional[PathLike] = None,
                 allow_new_keys: bool = False) -> None:
    """Write ``config.cfg``.

    With a template, its comments and ordering are preserved and only the values of the keys in
    ``params`` are replaced (a key missing from the template raises, unless ``allow_new_keys``, in
    which case it is appended). Without a template a bare ``key value`` file is written.
    """
    params = dict(params)
    lines: List[str] = []
    if template is not None:
        for raw in Path(template).read_text().splitlines():
            stripped = raw.strip()
            if stripped and not stripped.startswith("//"):
                key = stripped.split()[0]
                if key in params:
                    indent = raw[: len(raw) - len(raw.lstrip())]
                    lines.append("{}{:<24s} {}".format(indent, key, _format_value(params.pop(key))))
                    continue
            lines.append(raw)
        if params and not allow_new_keys:
            raise KeyError("keys not present in the template config: {}".format(sorted(params)))
        if params:
            lines.append("")
            lines.append("// Added by regolit.driver:")
    for key, value in params.items():
        lines.append("{:<24s} {}".format(key, _format_value(value)))
    Path(path).write_text("\n".join(lines) + "\n")


def read_layers(path: PathLike) -> List[Tuple[int, float, float, float, float]]:
    """Parse ``layers.cfg`` into ``(class, thickness, regolith, ice, soot)`` rows, file order."""
    rows = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        rows.append((int(float(parts[0])), *[float(v) for v in parts[1:5]]))
    return rows


def write_layers(path: PathLike, rows: Iterable[Sequence[float]]) -> None:
    """Write ``layers.cfg`` from ``(class, thickness, regolith, ice, soot)`` rows.

    Rows are bottom-up within a class (the first row of a class is the deepest layer), exactly as in
    the file format. Consecutive rows with the same class form one column class.
    """
    text = ["// Written by regolit.io.write_layers: class thickness regolith ice soot (bottom-up)"]
    for row in rows:
        text.append("{:d} {:.12g} {:.12g} {:.12g} {:.12g}".format(int(row[0]), *[float(v) for v in row[1:5]]))
    Path(path).write_text("\n".join(text) + "\n")


# ----------------------------------------------------------------------------------------------
# Binary maps
# ----------------------------------------------------------------------------------------------
def read_coordinates(output_dir: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Cell-center coordinates (m) of the output maps, from ``x.out`` and ``y.out``."""
    output_dir = Path(output_dir)
    x = np.fromfile(output_dir / "x.out", dtype="<f8")
    y = np.fromfile(output_dir / "y.out", dtype="<f8")
    return x, y


def read_matrix(path: PathLike, n: int) -> np.ndarray:
    """Read one ``n x n`` map of doubles and return it indexed ``[y, x]``."""
    data = np.fromfile(path, dtype="<f8")
    if data.size != n * n:
        raise ValueError("{} holds {} values, expected {}x{}".format(path, data.size, n, n))
    return data.reshape(n, n).T.copy()


def list_steps(output_dir: PathLike, prefix: str = "elevation") -> List[str]:
    """Sorted index strings of the printed steps, e.g. ``['01', ..., '11']``."""
    pattern = re.compile(r"^{}_(\d+)\.out$".format(re.escape(prefix)))
    steps = []
    for path in Path(output_dir).iterdir():
        match = pattern.match(path.name)
        if match:
            steps.append(match.group(1))
    return sorted(steps, key=int)


# ----------------------------------------------------------------------------------------------
# Layer stacks
# ----------------------------------------------------------------------------------------------
class Subsurface:
    """Full-resolution layer stacks of one time step.

    ``layers[j][i]`` is an ``(n_layers, 4)`` array of ``thickness, regolith, ice, soot`` for the
    column at ``y[j], x[i]``, ordered bottom-up (the last row is the surface layer). The bottom
    layer is the basement and is treated as extending indefinitely downward.
    """

    def __init__(self, elevation: np.ndarray, layers: List[List[np.ndarray]]):
        self.elevation = elevation
        self.layers = layers
        self.n = elevation.shape[0]

    def column(self, i: int, j: int) -> np.ndarray:
        """Layers of the column at x index ``i``, y index ``j`` (bottom-up)."""
        return self.layers[j][i]

    def interface_depths(self, i: int, j: int) -> np.ndarray:
        """Depths below the surface of the layer interfaces, from the surface downward."""
        thickness = self.layers[j][i][:, 0][::-1]
        return np.concatenate([[0.0], np.cumsum(thickness)])

    def composition_at_depth(self, depth: float) -> np.ndarray:
        """Composition ``[y, x, (regolith, ice, soot)]`` of the layer found ``depth`` m below the surface."""
        out = np.empty((self.n, self.n, 3))
        for j in range(self.n):
            for i in range(self.n):
                stack = self.layers[j][i]
                remaining = depth
                k = stack.shape[0] - 1
                while k > 0 and remaining >= stack[k, 0]:
                    remaining -= stack[k, 0]
                    k -= 1
                out[j, i] = stack[k, 1:4]
        return out

    def integrated_composition(self, depth: float) -> np.ndarray:
        """Thickness-weighted composition ``[y, x, 3]`` of the top ``depth`` m (as ``depthXFraction`` files)."""
        out = np.empty((self.n, self.n, 3))
        for j in range(self.n):
            for i in range(self.n):
                stack = self.layers[j][i]
                remaining = depth
                acc = np.zeros(3)
                for k in range(stack.shape[0] - 1, -1, -1):
                    take = remaining if k == 0 else min(remaining, stack[k, 0])
                    if take > 0:
                        acc += take * stack[k, 1:4]
                        remaining -= take
                    if remaining <= 0:
                        break
                total = acc.sum()
                out[j, i] = acc / total if total > 0 else acc
        return out

    def number_of_layers(self) -> np.ndarray:
        """Number of layers in every column, ``[y, x]``."""
        return np.array([[self.layers[j][i].shape[0] for i in range(self.n)] for j in range(self.n)])


def read_subsurface(path: PathLike) -> Subsurface:
    """Read ``subsurface_XX.out``: per column a header (number of layers, elevation, -1, -1) then the layers."""
    data = np.fromfile(path, dtype="<f8")
    columns = []
    offset = 0
    while offset < data.size:
        n_layers = int(data[offset])
        elevation = data[offset + 1]
        offset += 4
        columns.append((elevation, data[offset: offset + 4 * n_layers].reshape(n_layers, 4).copy()))
        offset += 4 * n_layers
    n = int(round(math.sqrt(len(columns))))
    if n * n != len(columns):
        raise ValueError("{} holds {} columns, not a square grid".format(path, len(columns)))
    elevation = np.empty((n, n))
    layers: List[List[np.ndarray]] = [[None] * n for _ in range(n)]  # type: ignore[list-item]
    for k, (z, stack) in enumerate(columns):
        i, j = divmod(k, n)  # written with the x index outer and the y index inner
        elevation[j, i] = z
        layers[j][i] = stack
    return Subsurface(elevation, layers)


# ----------------------------------------------------------------------------------------------
# Text outputs
# ----------------------------------------------------------------------------------------------
def read_histogram(path: PathLike) -> Tuple[np.ndarray, np.ndarray]:
    """Read a histogram file: ``bins`` (log-spaced edges) and ``counts`` (``counts[k]`` is for ``bins[k] <= v < bins[k+1]``)."""
    data = np.loadtxt(path, ndmin=2)
    return data[:, 0], data[:, 1].astype(int)


def read_craters(path: PathLike) -> Dict[str, np.ndarray]:
    """Read ``existing_craters.txt`` into arrays: x, y, diameter, depth, initial_depth (m)."""
    with open(path) as handle:
        header = [name.strip() for name in handle.readline().split(",")]
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    if data.size == 0:
        return {name: np.array([]) for name in header}
    return {name: data[:, k] for k, name in enumerate(header)}


# ----------------------------------------------------------------------------------------------
# One run
# ----------------------------------------------------------------------------------------------
class RegolitOutput:
    """The output of one REGOLIT run.

    ``path`` may be the ``output`` directory or a run directory containing ``output/`` (and, if
    available, ``config/config.cfg`` and ``log/log.txt``).
    """

    def __init__(self, path: PathLike, config: Optional[Dict[str, float]] = None):
        path = Path(path)
        if (path / "output" / "x.out").exists():
            self.run_dir: Optional[Path] = path
            self.output_dir = path / "output"
        elif (path / "x.out").exists():
            self.output_dir = path
            self.run_dir = path.parent if (path.parent / "config" / "config.cfg").exists() else None
        else:
            raise FileNotFoundError("no REGOLIT output (x.out) found in {}".format(path))

        if config is None and self.run_dir is not None and (self.run_dir / "config" / "config.cfg").exists():
            config = read_config(self.run_dir / "config" / "config.cfg")
        self.config: Dict[str, float] = dict(config or {})

        self.x, self.y = read_coordinates(self.output_dir)
        self.n = len(self.x)
        self.steps = list_steps(self.output_dir)
        if not self.steps:
            raise FileNotFoundError("no elevation_XX.out files in {}".format(self.output_dir))

    # -- bookkeeping ---------------------------------------------------------------------------
    def __repr__(self) -> str:
        return "RegolitOutput({!r}: {} steps, {}x{} cells)".format(str(self.output_dir), len(self.steps), self.n, self.n)

    @property
    def resolution(self) -> float:
        """Cell size of the output maps (m)."""
        return float(self.x[1] - self.x[0]) if self.n > 1 else float(self.config.get("resolution", float("nan")))

    @property
    def extent(self) -> List[float]:
        """``[x0, x1, y0, y1]`` for ``imshow``."""
        half = self.resolution / 2
        return [self.x[0] - half, self.x[-1] + half, self.y[0] - half, self.y[-1] + half]

    @property
    def times(self) -> np.ndarray:
        """Model time (Ma) of every step: the loop prints every ``printTimeStep`` from 0, the last print is at ``endTime``."""
        dt = self.config.get("printTimeStep")
        end = self.config.get("endTime")
        if dt is None or end is None:
            return np.arange(len(self.steps), dtype=float)
        times = np.array([k * dt for k in range(len(self.steps) - 1)] + [end], dtype=float)
        return np.minimum(times, end)

    def _step(self, step: Union[int, str]) -> str:
        """Resolve a step: an ``int`` is a position in ``steps`` (negative counts from the end); a
        ``str`` is the file index as written in the file names, matched by value (``"3"`` == ``"03"``)."""
        if isinstance(step, str):
            if step in self.steps:
                return step
            matches = [s for s in self.steps if int(s) == int(step)]
            if not matches:
                raise KeyError("step {!r} not in {}".format(step, self.steps))
            return matches[0]
        return self.steps[step]

    def _map(self, prefix: str, step: Union[int, str]) -> np.ndarray:
        return read_matrix(self.output_dir / "{}_{}.out".format(prefix, self._step(step)), self.n)

    # -- maps ----------------------------------------------------------------------------------
    def elevation(self, step: Union[int, str] = -1) -> np.ndarray:
        """Surface elevation (m), ``[y, x]``."""
        return self._map("elevation", step)

    def surface_fraction(self, species: str, step: Union[int, str] = -1) -> np.ndarray:
        """Fraction of ``species`` (regolith, ice or soot) in the surface layer, ``[y, x]``."""
        return self._map(SURFACE_FILES[species], step)

    def integrated_fraction(self, species: str, step: Union[int, str] = -1) -> np.ndarray:
        """Fraction of ``species`` integrated over the top ``depthToIntegrate`` m, ``[y, x]``."""
        return self._map(INTEGRATED_FILES[species], step)

    def elevation_series(self) -> np.ndarray:
        """All steps stacked: ``[step, y, x]``."""
        return np.stack([self.elevation(s) for s in self.steps])

    def fraction_series(self, species: str, kind: str = "surface") -> np.ndarray:
        """All steps of a composition map stacked: ``[step, y, x]``; ``kind`` is ``surface`` or ``integrated``."""
        reader = self.surface_fraction if kind == "surface" else self.integrated_fraction
        return np.stack([reader(species, s) for s in self.steps])

    # -- layer stacks ---------------------------------------------------------------------------
    def has_subsurface(self) -> bool:
        return (self.output_dir / "subsurface_{}.out".format(self.steps[-1])).exists()

    def subsurface(self, step: Union[int, str] = -1) -> Subsurface:
        """Full-resolution layer stacks (requires ``isPrintSubsurface 1``)."""
        return read_subsurface(self.output_dir / "subsurface_{}.out".format(self._step(step)))

    # -- text outputs --------------------------------------------------------------------------
    def histogram(self, name: str) -> Tuple[np.ndarray, np.ndarray]:
        """``craters``, ``impactors``, ``depth`` or ``existing_craters`` histogram: ``(bins, counts)``."""
        return read_histogram(self.output_dir / HISTOGRAM_FILES[name])

    def craters(self) -> Dict[str, np.ndarray]:
        """Craters still visible at the end of the run: x, y, diameter, depth, initial_depth (m)."""
        return read_craters(self.output_dir / "existing_craters.txt")

    def log(self) -> str:
        """Contents of ``log/log.txt`` (empty string if not available)."""
        if self.run_dir is not None and (self.run_dir / "log" / "log.txt").exists():
            return (self.run_dir / "log" / "log.txt").read_text()
        return ""
