"""Python tools for the REGOLIT surface evolution model.

- :mod:`regolit.io`     readers for the binary and text output files
- :mod:`regolit.driver` write a configuration, run the C++ binary, collect the results
- :mod:`regolit.plot`   quick-look figures

``run``, ``build``, ``quicklook`` and ``plot_histogram`` are imported lazily so that
``python -m regolit.driver`` and ``python -m regolit.plot`` work without warnings and matplotlib is
only loaded when a figure is requested.
"""
from .io import (RegolitOutput, Subsurface, read_config, write_config, read_layers, write_layers,
                 read_matrix, read_subsurface, read_histogram, read_craters)
from . import scaling

__all__ = ["RegolitOutput", "Subsurface", "read_config", "write_config", "read_layers", "write_layers",
           "read_matrix", "read_subsurface", "read_histogram", "read_craters", "scaling",
           "run", "build", "REPO_ROOT", "DEFAULT_BINARY", "quicklook", "plot_histogram"]

_LAZY = {"run": "driver", "build": "driver", "REPO_ROOT": "driver", "DEFAULT_BINARY": "driver",
         "quicklook": "plot", "plot_histogram": "plot"}


def __getattr__(name):
    if name in _LAZY:
        import importlib
        module = importlib.import_module("." + _LAZY[name], __name__)
        return getattr(module, name)
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))
