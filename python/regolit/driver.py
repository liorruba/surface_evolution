"""Run REGOLIT from Python.

The C++ binary reads ``./config/config.cfg`` and ``./config/layers.cfg`` (optionally
``./config/pixelIndex.cfg``) relative to its working directory and writes ``./output`` and
``./log`` there. :func:`run` prepares such a directory, executes the binary in it and returns a
:class:`~regolit.io.RegolitOutput`.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Union

import numpy as np

from .io import RegolitOutput, read_config, write_config, write_layers

PathLike = Union[str, Path]

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BINARY = REPO_ROOT / "build" / "apps" / "regolit_main.run"
DEFAULT_CONFIG = REPO_ROOT / "config" / "config.cfg"
DEFAULT_LAYERS = REPO_ROOT / "config" / "layers.cfg"


def build(repo_root: PathLike = REPO_ROOT, debug: bool = False, quiet: bool = True) -> Path:
    """Run ``make`` (or ``make debug``) in the repository and return the path of the binary."""
    target = "debug" if debug else "main"
    subprocess.run(["make", target], cwd=str(repo_root), check=True,
                   stdout=subprocess.DEVNULL if quiet else None)
    name = "regolit_main_debug.run" if debug else "regolit_main.run"
    return Path(repo_root) / "build" / "apps" / name


def prepare_run_directory(workdir: PathLike, overrides: Optional[Dict[str, object]] = None, *,
                          config_template: PathLike = DEFAULT_CONFIG,
                          layers: Union[PathLike, Iterable[Sequence[float]]] = DEFAULT_LAYERS,
                          pixel_index: Optional[np.ndarray] = None) -> Dict[str, float]:
    """Create ``workdir/config`` with the configuration of a run; returns the effective config."""
    workdir = Path(workdir)
    config_dir = workdir / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    write_config(config_dir / "config.cfg", overrides or {}, template=config_template)
    config = read_config(config_dir / "config.cfg")

    if isinstance(layers, (str, Path)):
        shutil.copyfile(str(layers), str(config_dir / "layers.cfg"))
    else:
        write_layers(config_dir / "layers.cfg", layers)

    mask_path = config_dir / "pixelIndex.cfg"
    if pixel_index is not None:
        grid_size = int(round(config["regionWidth"] / config["resolution"]))
        mask = np.asarray(pixel_index)
        if mask.shape != (grid_size, grid_size):
            raise ValueError("pixel_index must have shape ({0}, {0}) = (regionWidth/resolution)^2, got {1}".format(grid_size, mask.shape))
        # The C++ code reads the class of cell (x index i, y index j) at linear index j * gridSize + i,
        # i.e. the array is indexed [y, x] and written row-major.
        mask.astype(np.int8).tofile(str(mask_path))
    elif mask_path.exists():
        mask_path.unlink()
    return config


def run(overrides: Optional[Dict[str, object]] = None, workdir: Optional[PathLike] = None, *,
        binary: Optional[PathLike] = None, config_template: PathLike = DEFAULT_CONFIG,
        layers: Union[PathLike, Iterable[Sequence[float]]] = DEFAULT_LAYERS,
        pixel_index: Optional[np.ndarray] = None, build_if_missing: bool = True,
        quiet: bool = False, timeout: Optional[float] = None) -> RegolitOutput:
    """Run the model and return its output.

    Parameters
    ----------
    overrides
        ``{"parameter": value}`` replacing values of the template ``config.cfg``.
    workdir
        Directory for the run (``config/``, ``output/``, ``log/``). Default: ``runs/<timestamp>``
        under the repository root. An existing ``output/`` in it is removed by the binary.
    binary
        Path of the executable; default ``build/apps/regolit_main.run``, built if missing.
    config_template, layers
        Template ``config.cfg`` and ``layers.cfg`` (a path, or rows for :func:`write_layers`).
    pixel_index
        Optional ``int8`` array ``[y, x]`` of layer classes (see ``layers.cfg``).
    quiet
        Suppress the binary's progress output.
    """
    if workdir is None:
        workdir = REPO_ROOT / "runs" / _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    workdir = Path(workdir)

    binary = Path(binary) if binary is not None else DEFAULT_BINARY
    if not binary.exists():
        if not build_if_missing:
            raise FileNotFoundError("REGOLIT binary not found: {}".format(binary))
        binary = build(REPO_ROOT, quiet=quiet)

    config = prepare_run_directory(workdir, overrides, config_template=config_template, layers=layers, pixel_index=pixel_index)

    result = subprocess.run([str(binary.resolve())], cwd=str(workdir), timeout=timeout,
                            stdout=subprocess.PIPE if quiet else None, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        log_path = workdir / "log" / "log.txt"
        tail = "\n".join(log_path.read_text().splitlines()[-10:]) if log_path.exists() else ""
        raise RuntimeError("REGOLIT exited with status {}\n{}\n{}".format(result.returncode, result.stderr.strip(), tail))
    return RegolitOutput(workdir, config=config)


def _parse_value(text: str):
    try:
        value = float(text)
    except ValueError:
        return text
    return int(value) if value.is_integer() and "." not in text and "e" not in text.lower() else value


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run REGOLIT with modified parameters.")
    parser.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override a config parameter (repeatable)")
    parser.add_argument("--workdir", default=None, help="run directory (default: runs/<timestamp>)")
    parser.add_argument("--binary", default=None, help="path of the executable")
    parser.add_argument("--layers", default=str(DEFAULT_LAYERS), help="layers.cfg to use")
    parser.add_argument("--quiet", action="store_true", help="hide the progress output")
    parser.add_argument("--quicklook", action="store_true", help="save quicklook.png in the run directory")
    args = parser.parse_args(argv)

    overrides = {}
    for item in args.set:
        if "=" not in item:
            parser.error("--set expects KEY=VALUE, got {!r}".format(item))
        key, value = item.split("=", 1)
        overrides[key.strip()] = _parse_value(value.strip())

    out = run(overrides, args.workdir, binary=args.binary, layers=args.layers, quiet=args.quiet)
    print(out)
    print("times [Ma]:", np.array2string(out.times, precision=2))
    craters = out.craters()
    print("visible craters at the end: {}".format(len(craters["x"])))
    if args.quicklook:
        from .plot import quicklook
        target = Path(out.run_dir or out.output_dir.parent) / "quicklook.png"
        quicklook(out, save=target)
        print("saved", target)
    return 0


if __name__ == "__main__":
    sys.exit(main())
