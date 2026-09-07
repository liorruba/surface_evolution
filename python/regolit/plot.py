"""Quick-look figures for REGOLIT output."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np

from .io import RegolitOutput


def quicklook(out: RegolitOutput, step: Union[int, str] = -1, species: str = "soot", save: Optional[Union[str, Path]] = None,
              cross_section_y: Optional[float] = None):
    """Elevation, surface and integrated composition maps, and an elevation cross-section.

    The cross-section runs along x through ``cross_section_y`` (m); by default through the lowest
    point of the map. Returns the matplotlib figure.
    """
    import matplotlib.pyplot as plt

    z = out.elevation(step)
    surface = out.surface_fraction(species, step)
    integrated = out.integrated_fraction(species, step)
    extent = out.extent
    times = out.times
    index = out.steps.index(out._step(step))
    title_time = "t = {:g} Ma".format(times[index])

    if cross_section_y is None:
        j = int(np.unravel_index(np.argmin(z), z.shape)[0])
    else:
        j = int(np.argmin(np.abs(out.y - cross_section_y)))

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    panels = [
        (axes[0, 0], z, "terrain", "Elevation [m]"),
        (axes[0, 1], surface, "viridis", "Surface {} fraction".format(species)),
        (axes[1, 0], integrated, "viridis", "{} fraction, top {:g} m".format(species.capitalize(), out.config.get("depthToIntegrate", float("nan")))),
    ]
    for ax, data, cmap, label in panels:
        image = ax.imshow(data, origin="lower", extent=extent, cmap=cmap)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title("{} ({})".format(label, title_time))
        fig.colorbar(image, ax=ax, fraction=0.046)
    axes[0, 0].axhline(out.y[j], color="k", ls="--", lw=0.8)

    ax = axes[1, 1]
    ax.plot(out.x, z[j, :], color="C0")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("elevation [m]")
    ax.set_title("Cross-section at y = {:.1f} m".format(out.y[j]))
    ax.grid(alpha=0.3)

    fig.tight_layout()
    if save is not None:
        fig.savefig(str(save), dpi=110)
    return fig


def plot_histogram(out: RegolitOutput, name: str = "craters", ax=None, **kwargs):
    """Log-log step plot of one of the histograms (``craters``, ``impactors``, ``depth``, ``existing_craters``)."""
    import matplotlib.pyplot as plt

    bins, counts = out.histogram(name)
    if ax is None:
        _, ax = plt.subplots()
    ax.step(bins, counts, where="post", label=name, **kwargs)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("diameter [m]" if name != "depth" else "depth [m]")
    ax.set_ylabel("count")
    ax.legend()
    return ax


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Quick-look figure of a REGOLIT run.")
    parser.add_argument("path", help="output directory or run directory")
    parser.add_argument("--step", default="-1", help="file index of the step (1-based, as in elevation_XX.out), or a negative position counted from the end (default: -1, the last step)")
    parser.add_argument("--species", default="soot", choices=["regolith", "ice", "soot"])
    parser.add_argument("--save", default=None, help="write the figure to this file instead of showing it")
    args = parser.parse_args(argv)

    if args.save:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        number = int(args.step)
    except ValueError:
        parser.error("--step must be an integer")
    step: Union[int, str] = number if number < 0 else str(number)
    out = RegolitOutput(args.path)
    quicklook(out, step=step, species=args.species, save=args.save)
    if args.save:
        print("saved", args.save)
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
