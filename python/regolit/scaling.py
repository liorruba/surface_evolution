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
# Target material properties follow Holsapple (1993) as tabulated by Williams, Pathare & Aharonson
# (2014, Icarus 235, 23, Table 2): lunar regolith K1 = 0.132, K2 = 0.26, mu = 0.41, Y = 0.01 MPa,
# rho = 1500 kg/m^3; their nominal Mars target is dry soil with Y = 65 kPa and rho = 2000 kg/m^3.
# flux_ratio_to_earth is the body/Earth impactor flux ratio used with the fireball flux of
# Williams et al. (2014): 0.725 for the Moon (Ivanov 2006) and 1.885 for Mars (Hartmann 2005's 2.6
# times 0.725). Bodies without a published ratio use 1 and rely on the flux multiplier.
BODIES: Dict[str, Dict] = {
    "moon": {"label": "Moon", "g": 1.62, "meanImpactVelocity": 17500, "targetDensity": 1500, "Ybar": 1.0e4, "k2": 0.26,
             "flux_ratio_to_earth": 0.725, "note": "lunar regolith target (Williams et al. 2014, Table 2); mean impact velocity ~17.5 km/s"},
    "mercury": {"label": "Mercury", "g": 3.70, "meanImpactVelocity": 42000, "targetDensity": 1500, "Ybar": 1.0e4, "k2": 0.26,
                "flux_ratio_to_earth": None, "note": "lunar-regolith target properties; mean impact velocity ~42 km/s; no published Earth flux ratio (1 assumed)"},
    "mars": {"label": "Mars", "g": 3.71, "meanImpactVelocity": 10200, "targetDensity": 2000, "Ybar": 6.5e4, "k2": 0.26,
             "flux_ratio_to_earth": 1.885, "note": "dry-soil target and 10.2 km/s mean entry velocity after Williams et al. (2014); no atmosphere in this model"},
    "ceres": {"label": "Ceres", "g": 0.28, "meanImpactVelocity": 4800, "targetDensity": 1500, "Ybar": 1.0e4, "k2": 0.26,
              "flux_ratio_to_earth": None, "note": "lunar-regolith target properties; mean impact velocity ~4.8 km/s"},
    "vesta": {"label": "Vesta", "g": 0.25, "meanImpactVelocity": 4750, "targetDensity": 1500, "Ybar": 1.0e4, "k2": 0.26,
              "flux_ratio_to_earth": None, "note": "lunar-regolith target properties; mean impact velocity ~4.75 km/s"},
    "custom": {"label": "Custom", "flux_ratio_to_earth": None, "note": "keep the values entered below"},
}
BODY_FIELDS = ("g", "meanImpactVelocity", "targetDensity", "Ybar", "k2")

# ----------------------------------------------------------------------------------------------
# Production functions
# ----------------------------------------------------------------------------------------------
# Crater production functions give N(>D) of craters per m^2 per Ma at the present-day rate and are
# converted to impactor sizes through the crater scaling; the fireball flux is already an impactor
# size distribution and is used directly.
NEUKUM_COEFFICIENTS = [-3.0876, -3.557528, 0.781027, 1.021521, -0.156012, -0.444058, 0.019977,
                       0.086850, -0.005874, -0.006809, 8.25e-4, 5.54e-5]
EARTH_SURFACE_AREA_M2 = 5.10e14
KILOTON_J = 4.184e12

PRODUCTION_FUNCTIONS: Dict[str, Dict] = {
    "williams": {
        "label": "Williams et al. (2014): fireball flux scaled to the body",
        "kind": "impactor_power_law",
        "a0": 0.5677, "b0": 0.90,       # log10 N(>E) = a0 - b0 log10 E, bolides per year on Earth, E in kt (Brown et al. 2002)
        "reference": "Williams, Pathare & Aharonson (2014), Icarus 235, 23: the annual flux of terrestrial fireballs "
                     "(Brown et al. 2002, a0 = 0.5677, b0 = 0.90) converted to impactor diameters with the impactor density and the "
                     "body's mean impact velocity, and scaled by the body/Earth flux ratio (Moon 0.725, Mars 1.885). "
                     "An impactor size distribution with cumulative slope 3 b0 = 2.7; no atmospheric filtering.",
        "available": True,
    },
    "neukum": {
        "label": "Neukum et al. (2001), lunar",
        "kind": "log10_polynomial",
        "coefficients": NEUKUM_COEFFICIENTS,        # log10 N(>D) [km^-2 Ga^-1] = sum a_j (log10 D_km)^j
        "valid_range_m": [10.0, 300000.0],
        "reference": "Neukum, Ivanov & Hartmann (2001), Space Sci. Rev. 96, 55. Lunar production function for "
                     "1 Ga at the present rate; extrapolated below 10 m craters.",
        "available": True,
    },
    "daubar": {
        "label": "Daubar et al. (2013), current Mars",
        "kind": "power_law",
        "reference_diameter_m": 3.9,
        "reference_rate": 1.65e-6 * 1e-6 * 1e6,     # 1.65e-6 craters km^-2 yr^-1 at D >= 3.9 m -> per m^2 per Ma
        "cumulative_slope": 1.45,                   # differential slope -2.45 +- 0.36 -> cumulative slope 1.45
        "reference": "Daubar et al. (2013), Icarus 225, 506: present-day rate of 1.65e-6 craters km^-2 yr^-1 with effective "
                     "D >= 3.9 m, differential slope -2.45 +- 0.36 (cumulative slope 1.45). Primary craters on Mars, "
                     "measured under the present atmosphere.",
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
    """Cumulative crater density N(>D) per m^2 per Ma for a crater production function."""
    spec = PRODUCTION_FUNCTIONS[key]
    D = np.asarray(diameter_m, dtype=float)
    if spec["kind"] == "log10_polynomial":
        log_d_km = np.log10(D / 1000.0)
        log_n = sum(a * log_d_km ** j for j, a in enumerate(spec["coefficients"]))
        return 10.0 ** log_n * 1e-9            # km^-2 Ga^-1 -> m^-2 Ma^-1
    if spec["kind"] == "power_law":
        return spec["reference_rate"] * (D / spec["reference_diameter_m"]) ** (-spec["cumulative_slope"])
    raise ValueError("production function {} is not a crater production function".format(key))


def body_flux_ratio(body: Optional[str]) -> float:
    ratio = BODIES.get(body or "custom", {}).get("flux_ratio_to_earth")
    return float(ratio) if ratio else 1.0


# ----------------------------------------------------------------------------------------------
# Crater scaling (mirrors crater.cpp)
# ----------------------------------------------------------------------------------------------
def transient_volume(impactor_diameter: float, p: Dict[str, float]) -> float:
    a = impactor_diameter / 2.0
    v = p["meanImpactVelocity"]
    mass = 4.0 / 3.0 * math.pi * a ** 3 * p["impactorDensity"]
    buff1 = (p["g"] * a / v ** 2) * (p["targetDensity"] / p["impactorDensity"]) ** (-1.0 / 3.0)
    buff2 = (p["Ybar"] / p["targetDensity"] / v ** 2) ** ((2.0 + p["mu"]) / 2.0)
    return p["k1"] * (mass / p["targetDensity"]) * (buff1 + p.get("k2", 1.0) * buff2) ** (-3 * p["mu"] / (2 + p["mu"]))


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


def impactor_density(key: str, impactor_diameters, p: Dict[str, float], body: Optional[str] = None) -> np.ndarray:
    """N(>d) per m^2 per Ma of impactors larger than d, for any production function."""
    d = np.asarray(impactor_diameters, dtype=float)
    spec = PRODUCTION_FUNCTIONS[key]
    if spec["kind"] == "impactor_power_law":
        mass = math.pi / 6.0 * d ** 3 * p["impactorDensity"]
        energy_kt = 0.5 * mass * p["meanImpactVelocity"] ** 2 / KILOTON_J
        per_earth_per_year = 10.0 ** spec["a0"] * energy_kt ** (-spec["b0"])
        return per_earth_per_year / EARTH_SURFACE_AREA_M2 * 1e6 * body_flux_ratio(body)
    craters = np.array([2 * final_crater_radius(float(x), p) for x in d])
    return crater_density(key, craters)


# ----------------------------------------------------------------------------------------------
# Fits and estimates
# ----------------------------------------------------------------------------------------------
def fit_power_law(key: str, p: Dict[str, float], dmin: float, dmax: float, points: int = 40, body: Optional[str] = None) -> Tuple[float, float]:
    """Least-squares power law N(>d) = c d^-b through the production function over [dmin, dmax]."""
    d = np.logspace(math.log10(dmin), math.log10(max(dmax, dmin * 1.5)), points)
    n = impactor_density(key, d, p, body)
    mask = n > 0
    slope, intercept = np.polyfit(np.log(d[mask]), np.log(n[mask]), 1)
    return math.exp(intercept), -slope


def largest_expected_impactor(key: str, p: Dict[str, float], area: float, time: float, multiplier: float, dmin: float, body: Optional[str] = None) -> float:
    """Impactor diameter for which one impact is expected over the area and time (bisection on the production function)."""
    def expected(d):
        return float(impactor_density(key, [d], p, body)[0]) * area * time * multiplier
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


def estimate(p: Dict[str, float], layers_total: float, production_function: str = "power_law", body: Optional[str] = None) -> Dict:
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
        dmax_pf = largest_expected_impactor(production_function, p, area, time, multiplier, dmin, body)
        c, b = fit_power_law(production_function, p, dmin, dmax_pf, body=body)
        fit_range = [dmin, max(dmax_pf, dmin * 1.5)]

    impacts = c * dmin ** (-b) * area * time * multiplier
    largest_impactor = (c * area * time * multiplier) ** (1.0 / b) if impacts > 1 else dmin
    largest_radius = final_crater_radius(largest_impactor, p)
    largest_depth = p["depthToDiameter"] * 2 * largest_radius
    smallest_crater = 2 * final_crater_radius(dmin, p)
    basement = max(10.0, math.ceil((layers_total + 3.0 * largest_depth) * 10) / 10)
    return {
        "production_function": production_function,
        "body": body,
        "body_flux_ratio": body_flux_ratio(body) if spec["kind"] == "impactor_power_law" else None,
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
