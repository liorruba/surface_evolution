"""Crater scaling and impactor production functions.

The crater scaling mirrors ``src/crater.cpp`` (Holsapple pi-scaling, final radius 1.18 x transient)
so that estimates made here agree with what the model forms. Production functions give the
cumulative number of *craters* larger than a diameter; they are converted to impactor size
distributions through that scaling and approximated by the single power law
``N(>d) = c d^-b`` that the model samples from, fitted over the impactor size range of a run.

Units: lengths in m, densities in kg/m^3, velocities in m/s, cumulative densities per m^2 per Ma.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

# ----------------------------------------------------------------------------------------------
# Target bodies (editable defaults; the angle of repose is deliberately not part of a preset)
# ----------------------------------------------------------------------------------------------
BODIES: Dict[str, Dict] = {
    "moon": {"label": "Moon", "g": 1.62, "meanImpactVelocity": 17500, "targetDensity": 1700,
             "note": "mean impact velocity of the present asteroid flux ~17-20 km/s"},
    "mercury": {"label": "Mercury", "g": 3.70, "meanImpactVelocity": 42000, "targetDensity": 1700,
                "note": "mean impact velocity ~42 km/s"},
    "mars": {"label": "Mars", "g": 3.71, "meanImpactVelocity": 9600, "targetDensity": 1700,
             "note": "mean impact velocity ~9.6 km/s"},
    "ceres": {"label": "Ceres", "g": 0.28, "meanImpactVelocity": 4800, "targetDensity": 1500,
              "note": "mean impact velocity ~4.8 km/s"},
    "vesta": {"label": "Vesta", "g": 0.25, "meanImpactVelocity": 4750, "targetDensity": 1700,
              "note": "mean impact velocity ~4.75 km/s"},
    "custom": {"label": "Custom", "note": "keep the values entered below"},
}
BODY_FIELDS = ("g", "meanImpactVelocity", "targetDensity")

# ----------------------------------------------------------------------------------------------
# Production functions: cumulative crater density N(>D) per m^2 per Ma at the present-day rate
# ----------------------------------------------------------------------------------------------
NEUKUM_COEFFICIENTS = [-3.0876, -3.557528, 0.781027, 1.021521, -0.156012, -0.444058, 0.019977,
                       0.086850, -0.005874, -0.006809, 8.25e-4, 5.54e-5]

PRODUCTION_FUNCTIONS: Dict[str, Dict] = {
    "neukum": {
        "label": "Neukum et al. (2001), lunar",
        "kind": "log10_polynomial",
        "coefficients": NEUKUM_COEFFICIENTS,        # log10 N(>D) [km^-2 Ga^-1] = sum a_j (log10 D_km)^j
        "valid_range_m": [10.0, 300000.0],
        "reference": "Neukum, Ivanov & Hartmann (2001), Space Sci. Rev. 96, 55. Lunar production function for "
                     "1 Ga at the present rate; extrapolated below 10 m craters.",
        "available": True,
    },
    "marchi": {
        "label": "Marchi et al. (2009)",
        "kind": "log10_polynomial",
        "coefficients": None,
        "valid_range_m": None,
        "reference": "Marchi et al. (2009), AJ 137, 4936. Coefficients not entered yet: add them here to enable.",
        "available": False,
    },
    "daubar": {
        "label": "Daubar et al. (2013), current Mars",
        "kind": "power_law",
        "reference_diameter_m": 3.9,
        "reference_rate": 1.65e-6 * 1e-6 * 1e6,     # 1.65e-6 craters km^-2 yr^-1 at D >= 3.9 m -> per m^2 per Ma
        "cumulative_slope": 2.5,
        "reference": "Daubar et al. (2013), Icarus 225, 506: present-day rate 1.65e-6 km^-2 yr^-1 for D >= 3.9 m. "
                     "The cumulative slope entered here (2.5) should be checked against the paper.",
        "available": True,
    },
    "power_law": {
        "label": "Power law (manual)",
        "kind": "manual",
        "reference": "N(>d) = c d^-b with the flux constant c and slope b entered by hand.",
        "available": True,
    },
}


def crater_density(key: str, diameter_m) -> np.ndarray:
    """Cumulative crater density N(>D) per m^2 per Ma for a production function."""
    spec = PRODUCTION_FUNCTIONS[key]
    D = np.asarray(diameter_m, dtype=float)
    if spec["kind"] == "log10_polynomial":
        if not spec.get("coefficients"):
            raise ValueError("production function {} has no coefficients".format(key))
        log_d_km = np.log10(D / 1000.0)
        log_n = sum(a * log_d_km ** j for j, a in enumerate(spec["coefficients"]))
        return 10.0 ** log_n * 1e-9            # km^-2 Ga^-1 -> m^-2 Ma^-1
    if spec["kind"] == "power_law":
        return spec["reference_rate"] * (D / spec["reference_diameter_m"]) ** (-spec["cumulative_slope"])
    raise ValueError("production function {} is not evaluable".format(key))


# ----------------------------------------------------------------------------------------------
# Crater scaling (mirrors crater.cpp)
# ----------------------------------------------------------------------------------------------
def transient_volume(impactor_diameter: float, p: Dict[str, float]) -> float:
    a = impactor_diameter / 2.0
    v = p["meanImpactVelocity"]
    mass = 4.0 / 3.0 * math.pi * a ** 3 * p["impactorDensity"]
    buff1 = (p["g"] * a / v ** 2) * (p["targetDensity"] / p["impactorDensity"]) ** (-1.0 / 3.0)
    buff2 = (p["Ybar"] / p["targetDensity"] / v ** 2) ** ((2.0 + p["mu"]) / 2.0)
    return p["k1"] * (mass / p["targetDensity"]) * (buff1 + buff2) ** (-3 * p["mu"] / (2 + p["mu"]))


def final_crater_radius(impactor_diameter: float, p: Dict[str, float]) -> float:
    """Final crater radius (m) formed by an impactor of the given diameter (m)."""
    return 1.18 * (3 * transient_volume(impactor_diameter, p) / math.pi) ** (1.0 / 3.0)


def impactor_diameter_for_crater(crater_diameter: float, p: Dict[str, float]) -> float:
    """Inverse of final_crater_radius (crater diameter in m), by bisection in log space."""
    lo, hi = 1e-4, 1e5
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if 2 * final_crater_radius(mid, p) < crater_diameter:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def impactor_density(key: str, impactor_diameters, p: Dict[str, float]) -> np.ndarray:
    """N(>d) per m^2 per Ma of impactors, from a crater production function and the scaling."""
    d = np.asarray(impactor_diameters, dtype=float)
    craters = np.array([2 * final_crater_radius(float(x), p) for x in d])
    return crater_density(key, craters)


# ----------------------------------------------------------------------------------------------
# Fits and estimates
# ----------------------------------------------------------------------------------------------
def fit_power_law(key: str, p: Dict[str, float], dmin: float, dmax: float, points: int = 40) -> Tuple[float, float]:
    """Least-squares power law N(>d) = c d^-b through the production function over [dmin, dmax]."""
    d = np.logspace(math.log10(dmin), math.log10(max(dmax, dmin * 1.5)), points)
    n = impactor_density(key, d, p)
    mask = n > 0
    slope, intercept = np.polyfit(np.log(d[mask]), np.log(n[mask]), 1)
    return math.exp(intercept), -slope


def largest_expected_impactor(key: str, p: Dict[str, float], area: float, time: float, multiplier: float, dmin: float) -> float:
    """Impactor diameter for which one impact is expected over the area and time (bisection on the production function)."""
    def expected(d):
        return float(impactor_density(key, [d], p)[0]) * area * time * multiplier
    if expected(dmin) <= 1:
        return dmin
    lo, hi = dmin, 1e5
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if expected(mid) > 1:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def estimate(p: Dict[str, float], layers_total: float, production_function: str = "power_law") -> Dict:
    """Derived quantities for a parameter set: fitted flux law, expected impacts, largest crater, basement thickness.

    ``p`` needs regionWidth, endTime, minimumImpactorDiameter, earthFluxRatioCoefficient, fluxConstant_c,
    slope_b, depthToDiameter and the scaling parameters (g, k1, mu, Ybar, targetDensity, impactorDensity,
    meanImpactVelocity).
    """
    area = p["regionWidth"] ** 2
    time = p["endTime"]
    multiplier = p["earthFluxRatioCoefficient"]
    dmin = p["minimumImpactorDiameter"]
    spec = PRODUCTION_FUNCTIONS.get(production_function, PRODUCTION_FUNCTIONS["power_law"])

    if spec["kind"] == "manual" or not spec.get("available"):
        c, b = p["fluxConstant_c"], p["slope_b"]
        fit_range = None
    else:
        dmax_pf = largest_expected_impactor(production_function, p, area, time, multiplier, dmin)
        c, b = fit_power_law(production_function, p, dmin, dmax_pf)
        fit_range = [dmin, max(dmax_pf, dmin * 1.5)]

    impacts = c * dmin ** (-b) * area * time * multiplier
    largest_impactor = (c * area * time * multiplier) ** (1.0 / b) if impacts > 1 else dmin
    largest_radius = final_crater_radius(largest_impactor, p)
    largest_depth = p["depthToDiameter"] * 2 * largest_radius
    smallest_crater = 2 * final_crater_radius(dmin, p)
    basement = max(10.0, math.ceil((layers_total + 3.0 * largest_depth) * 10) / 10)
    return {
        "production_function": production_function,
        "fluxConstant_c": c,
        "slope_b": b,
        "fit_range_m": fit_range,
        "expected_impacts": impacts,
        "largest_impactor_m": largest_impactor,
        "largest_crater_diameter_m": 2 * largest_radius,
        "largest_crater_depth_m": largest_depth,
        "smallest_crater_diameter_m": smallest_crater,
        "layers_total_m": layers_total,
        "suggested_basement_m": basement,
    }
