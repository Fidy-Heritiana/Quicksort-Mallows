"""Read-only numerical helpers shared by the application and notebook.

The original implementation lives in reproduce.py. This module calls those
functions and never replaces the manuscript data or figures.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import norm

from reproduce import finite_moments, sample_costs

ROOT = Path(__file__).resolve().parent
SAVED_PARAMETERS = ("0.25", "0.5", "0.8")
COLOURS = {"0.25": "#0072B2", "0.5": "#D55E00", "0.8": "#009E73"}


def load_manuscript(phi: str = "0.5") -> dict:
    """Load one saved experiment, including its original seed and sample size."""
    if phi not in SAVED_PARAMETERS:
        raise ValueError("Choose one of the three saved parameters.")
    metadata = json.loads((ROOT / "data/validation.json").read_text())
    with np.load(ROOT / f"data/moments_phi_{phi}.npz", allow_pickle=False) as data:
        mean = data["mean"].copy()
        moments = data["centred_moments"].copy()
    with np.load(ROOT / "data/simulated_costs.npz", allow_pickle=False) as data:
        costs = data[f"phi_{phi}"].copy()
    return {
        "phi": phi, "n": metadata["sample_n"],
        "repetitions": len(costs), "seed": metadata["parameters"][phi]["seed"],
        "mean": mean, "moments": moments, "costs": costs,
        "source": "Original manuscript data",
        "centring_residual": metadata["parameters"][phi]["max_relative_centring_residual"],
    }


def run_experiment(phi: str = "0.5", n: int = 200,
                   repetitions: int = 2000, seed: int = 20260922) -> dict:
    """Run a bounded interactive experiment without writing to disk.

    Bounds keep shared demos responsive and exclude zero-variance endpoints.
    Larger computations remain available through reproduce.py.
    """
    q = float(phi)
    if not math.isfinite(q) or not 0.01 <= q <= 0.99:
        raise ValueError("Interactive experiments require 0.01 <= phi <= 0.99.")
    for label, value, lower, upper in (
        ("n", n, 3, 1000),
        ("repetitions", repetitions, 500, 20000),
        ("seed", seed, 0, 2**32 - 1),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError(f"{label} must be an integer.")
        if not lower <= value <= upper:
            raise ValueError(f"{label} must be between {lower} and {upper}.")
    phi = str(phi)
    mean, moments, residual = finite_moments(phi, max_n=n, max_order=4)
    costs = sample_costs(q, n, repetitions, seed)
    if not np.isfinite(moments[:, n]).all() or moments[2, n] <= 0:
        raise ArithmeticError("Finite moments could not be normalised reliably.")
    return {
        "phi": phi, "n": n, "repetitions": repetitions, "seed": seed,
        "mean": mean, "moments": moments, "costs": costs,
        "source": "New simulation", "centring_residual": residual,
    }


def experiment_metadata(result: dict) -> dict:
    n = result["n"]
    return {
        "model": "First-pivot Quicksort with stable partition on Mallows input",
        "phi": float(result["phi"]), "n": n,
        "repetitions": result["repetitions"], "seed": result["seed"],
        "source": result["source"],
        "finite_mean": float(result["mean"][n]),
        "finite_variance": float(result["moments"][2, n]),
        "max_centred_moment_order": result["moments"].shape[0] - 1,
        "max_relative_centring_residual": result["centring_residual"],
        "standardisation": "(cost - finite_mean) / sqrt(finite_variance)",
    }


def standardised_costs(result: dict) -> np.ndarray:
    n = result["n"]
    return (result["costs"] - float(result["mean"][n])) / math.sqrt(
        float(result["moments"][2, n]))


def summary_table(result: dict) -> pd.DataFrame:
    n, costs = result["n"], result["costs"]
    return pd.DataFrame({
        "Quantity": ["Mean", "Variance"],
        "Finite recurrence": [float(result["mean"][n]), float(result["moments"][2, n])],
        "Monte Carlo": [float(np.mean(costs)), float(np.var(costs, ddof=1))],
    })


def sample_table(result: dict) -> pd.DataFrame:
    return pd.DataFrame({
        "phi": float(result["phi"]), "n": result["n"], "seed": result["seed"],
        "replicate": np.arange(1, result["repetitions"] + 1),
        "cost": result["costs"], "standardised_cost": standardised_costs(result),
    })


def moment_table(result: dict) -> pd.DataFrame:
    """Centred ordinary moments and normalised ratios; no factorial moments."""
    mean, m = result["mean"], result["moments"]
    n = np.arange(3, len(mean))
    q = float(result["phi"])
    variance = np.asarray(m[2, 3:], dtype=float)
    table = pd.DataFrame({
        "phi": q, "n": n, "mean": np.asarray(mean[3:], dtype=float),
        "variance": variance,
        "mean_ratio": np.asarray(mean[3:], dtype=float) / ((1 - q) * n**2 / 2),
        "variance_ratio": variance / (q * (1 - q) * n**3 / 3),
        "skewness": np.asarray(m[3, 3:], dtype=float) / variance**1.5,
        "kurtosis": np.asarray(m[4, 3:], dtype=float) / variance**2,
    })
    if m.shape[0] > 6:
        table["sixth_standardised_moment"] = np.asarray(m[6, 3:], dtype=float) / variance**3
    return table


def normal_figures(result: dict, bins: int = 45) -> tuple[go.Figure, go.Figure]:
    """Histogram and normal Q-Q plot, centred using finite theoretical moments."""
    z = standardised_costs(result)
    colour = COLOURS.get(result["phi"], "#0072B2")
    x = np.linspace(min(-4.0, float(z.min())), max(4.0, float(z.max())), 500)
    hist = go.Figure()
    hist.add_trace(go.Histogram(x=z, nbinsx=bins, histnorm="probability density",
                               name="Simulated costs", marker_color=colour, opacity=0.65))
    hist.add_trace(go.Scatter(x=x, y=norm.pdf(x), mode="lines", name="Standard normal",
                             line={"color": "#263544", "width": 2.5}))
    hist.update_layout(title="Normalised cost distribution", xaxis_title="(Cₙ − μₙ) / σₙ",
                       yaxis_title="Density", bargap=0.04, template="plotly_white",
                       legend={"orientation": "h", "y": -0.22}, margin={"t": 55, "b": 80})
    ordered = np.sort(z)
    indexes = np.unique(np.linspace(0, len(z) - 1, min(900, len(z)), dtype=int))
    theoretical = norm.ppf((indexes + 0.5) / len(z))
    qq = go.Figure()
    qq.add_trace(go.Scatter(x=theoretical, y=ordered[indexes], mode="markers",
                           name="Cost quantiles", marker={"color": colour, "size": 4}))
    low = min(float(theoretical.min()), float(ordered.min()))
    high = max(float(theoretical.max()), float(ordered.max()))
    qq.add_trace(go.Scatter(x=[low, high], y=[low, high], mode="lines", name="y = x",
                           line={"color": "#596977", "dash": "dash"}))
    qq.update_layout(title="Normal Q–Q plot", xaxis_title="Theoretical normal quantile",
                     yaxis_title="Standardised cost quantile", template="plotly_white",
                     legend={"orientation": "h", "y": -0.22}, margin={"t": 55, "b": 80})
    return hist, qq


METRICS = {
    "mean_ratio": ("Mean / leading quadratic term", 1.0),
    "variance_ratio": ("Variance / leading cubic term", 1.0),
    "skewness": ("Standardised third centred moment", 0.0),
    "kurtosis": ("Standardised fourth centred moment", 3.0),
    "sixth_standardised_moment": ("Standardised sixth centred moment", 15.0),
}


def convergence_figure(table: pd.DataFrame, metric: str) -> go.Figure:
    label, limit = METRICS[metric]
    fig = go.Figure()
    for q, group in table.groupby("phi"):
        key = str(q)
        fig.add_trace(go.Scatter(x=group["n"], y=group[metric], mode="lines",
                                name=f"φ = {q:g}", line={"color": COLOURS.get(key, "#0072B2")}))
    fig.add_hline(y=limit, line_dash="dash", line_color="#596977",
                  annotation_text=f"Asymptotic limit: {limit:g}")
    fig.update_layout(title=label, xaxis_title="Input size n", yaxis_title=label,
                       template="plotly_white", legend={"orientation": "h"})
    return fig
