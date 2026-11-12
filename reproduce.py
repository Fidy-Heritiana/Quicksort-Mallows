#!/usr/bin/env python3
"""Reproduce the numerical tables and figures of the detailed Quicksort-Mallows manuscript.

English translation of the original reproduction script.
Calculations, parameters, random seeds and saved data are unchanged.
All plot titles and axis labels are in English.

Run from any directory: python3 reproduce.py
The deterministic calculations use the full finite recurrence, not sampling.
The independent Monte Carlo experiment samples recursive Mallows splits.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from fractions import Fraction
from itertools import permutations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
PARAMETERS = ("0.25", "0.5", "0.8")
COLOURS = ("#0072B2", "#D55E00", "#009E73")
LINESTYLES = ("-", "--", "-.")
SEED = 20260921


def deterministic_cost(word):
    """Literal first-pivot Quicksort with stable sublists."""
    if len(word) < 2:
        return 0
    pivot = word[0]
    left = [x for x in word[1:] if x < pivot]
    right = [x for x in word[1:] if x > pivot]
    return len(word) - 1 + deterministic_cost(left) + deterministic_cost(right)


def triples(order):
    return [(a, b, order-a-b, math.comb(order, a)*math.comb(order-a, b))
            for a in range(order+1) for b in range(order-a+1)
            if a != 1 and b != 1]


def exhaustive_check(max_n=7):
    """Exact rational comparison with all permutations, phi=1/2."""
    q = Fraction(1, 2)
    p = 1-q
    order = 4
    mean = [Fraction(0)]*(max_n+1)
    moments = [[Fraction(0)]*(max_n+1) for _ in range(order+1)]
    moments[0] = [Fraction(1)]*(max_n+1)
    count = 0
    for n in range(1, max_n+1):
        weights = [p*q**j/(1-q**n) for j in range(n)]
        mean[n] = n-1+sum(weights[j]*(mean[j]+mean[n-1-j]) for j in range(n))
        exact_mean = sum((n-k)*p*(1+q**k)/(1-q**(k+1)) for k in range(1, n))
        assert mean[n] == exact_mean
        for r in range(2, order+1):
            value = Fraction(0)
            for j, weight in enumerate(weights):
                k = n-1-j
                d = n-1+mean[j]+mean[k]-mean[n]
                for a, b, c, coeff in triples(r):
                    value += weight*coeff*moments[a][j]*moments[b][k]*d**c
            moments[r][n] = value
        partition = Fraction(0)
        totals = [Fraction(0)]*(order+1)
        for word in permutations(range(n)):
            inv = sum(word[i] > word[j] for i in range(n) for j in range(i+1, n))
            weight = q**inv
            partition += weight
            centred = Fraction(deterministic_cost(word))-mean[n]
            for r in range(2, order+1):
                totals[r] += weight*centred**r
            count += 1
        assert all(totals[r]/partition == moments[r][n] for r in range(2, order+1))
    return {"phi": "1/2", "max_n": max_n, "permutations_checked": count,
            "centred_moment_orders": [2, 3, 4], "status": "passed"}


def finite_moments(q_string, max_n=2000, max_order=6):
    """Full finite sums in extended precision; no asymptotic substitutions."""
    q = np.longdouble(q_string)
    p = 1-q
    powers = q**np.arange(max_n+1, dtype=np.longdouble)
    mean = np.zeros(max_n+1, dtype=np.longdouble)
    moments = np.zeros((max_order+1, max_n+1), dtype=np.longdouble)
    moments[0, :] = 1
    terms = {r: triples(r) for r in range(2, max_order+1)}
    marginal = np.longdouble(0)
    max_relative_centring_residual = np.longdouble(0)
    for n in range(1, max_n+1):
        if n >= 2:
            marginal += p*(1+powers[n-1])/(1-powers[n])
        mean[n] = mean[n-1]+marginal
        weights = p*powers[:n]/(1-powers[n])
        d = n-1+mean[:n]+mean[n-1::-1]-mean[n]
        residual = abs(np.dot(weights, d))/max(1, n*n)
        max_relative_centring_residual = max(max_relative_centring_residual, residual)
        d_powers = [np.ones(n, dtype=np.longdouble)]
        for _ in range(max_order):
            d_powers.append(d_powers[-1]*d)
        for r in range(2, max_order+1):
            summand = np.zeros(n, dtype=np.longdouble)
            for a, b, c, coeff in terms[r]:
                summand += coeff*moments[a, :n]*moments[b, n-1::-1]*d_powers[c]
            moments[r, n] = np.dot(weights, summand)
    return mean, moments, float(max_relative_centring_residual)


def sample_costs(q, n, repetitions, seed):
    """Independent copies of the exact recursive cost law, batched by depth."""
    rng = np.random.default_rng(seed)
    costs = np.zeros(repetitions, dtype=np.int64)
    sizes = np.full(repetitions, n, dtype=np.int64)
    owners = np.arange(repetitions, dtype=np.int64)
    log_q = math.log(q)
    while sizes.size:
        np.add.at(costs, owners, sizes-1)
        uniforms = rng.random(sizes.size)
        left = np.floor(np.log1p(-uniforms*(1-np.power(q, sizes)))/log_q).astype(np.int64)
        left = np.minimum(left, sizes-1)  # Guard the floating-point endpoint.
        right = sizes-1-left
        left_mask = left >= 2
        right_mask = right >= 2
        sizes = np.concatenate((left[left_mask], right[right_mask]))
        owners = np.concatenate((owners[left_mask], owners[right_mask]))
    return costs


def make_figures(results, samples, sample_n, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelsize": 9, "legend.fontsize": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "savefig.bbox": "tight"})
    fig, axes = plt.subplots(1, 2, figsize=(6.25, 2.55), constrained_layout=True)
    for key, colour, style in zip(PARAMETERS, COLOURS, LINESTYLES):
        mean, moment, _ = results[key]
        q = float(key)
        n = np.arange(20, len(mean))
        v = q*(1-q)/3
        ratio = np.asarray(moment[2, 20:]/(v*n**3), dtype=float)
        axes[0].plot(n, ratio, color=colour, linestyle=style, label=rf"$\phi={key}$")
        z = np.sort((samples[key]-float(mean[sample_n]))/math.sqrt(float(moment[2, sample_n])))
        indexes = np.linspace(12, len(z)-13, 450, dtype=int)
        quantiles = norm.ppf((indexes+0.5)/len(z))
        axes[1].plot(quantiles, z[indexes], color=colour, linestyle=style, linewidth=1.05)
    axes[0].axhline(1, color="#555555", linewidth=0.9, linestyle=":")
    axes[0].set(xlabel="$n$", ylabel=r"$\mathrm{Var}(C_n)/(v_\phi n^3)$", title="Finite-size variance")
    axes[0].legend(loc="lower right", frameon=False)
    axes[1].plot([-3.2, 3.2], [-3.2, 3.2], ":", color="#555555", linewidth=0.9)
    axes[1].set(xlabel="Theoretical normal quantile", ylabel="Standardised cost quantile",
                title=rf"Normal Q-Q plot, $n={sample_n}$", xlim=(-3.2, 3.2), ylim=(-3.4, 3.4))
    for ax in axes:
        ax.grid(alpha=0.15)
    fig.savefig(out/"variance_and_qq.pdf")
    fig.savefig(out/"variance_and_qq.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(6.25, 2.05), constrained_layout=True)
    x = np.linspace(-4, 4, 400)
    for ax, key, colour in zip(axes, PARAMETERS, COLOURS):
        mean, moment, _ = results[key]
        z = (samples[key]-float(mean[sample_n]))/math.sqrt(float(moment[2, sample_n]))
        ax.hist(z, bins=np.linspace(-4.5, 4.5, 46), density=True,
                color=colour, alpha=0.48, edgecolor="white", linewidth=0.4)
        ax.plot(x, norm.pdf(x), color="#292929", linewidth=1.2)
        ax.set(xlabel="Standardised cost", title=rf"$\phi={key}$", xlim=(-4, 4), ylim=(0, 0.48))
    axes[0].set_ylabel("Density")
    fig.savefig(out/"normal_histograms.pdf")
    fig.savefig(out/"normal_histograms.png", dpi=220)
    plt.close(fig)


def write_latex_rows(table, path):
    with path.open("w") as stream:
        stream.write("\\newcommand{\\NumericalRows}{%\n")
        for row in table:
            stream.write(f"{row['phi']} & {int(row['n'])} & {float(row['variance_ratio']):.6f} & "
                         f"{float(row['skewness']):.6f} & {float(row['kurtosis']):.6f} & "
                         f"{float(row['sixth_standardised_moment']):.6f} \\\\\n")
        stream.write("}\n")


def redraw_saved(data, figures):
    """Regenerate the figures from the saved numerical data."""
    summary = json.loads((data/"validation.json").read_text())
    results, samples = {}, {}
    with np.load(data/"simulated_costs.npz") as archive:
        for key in PARAMETERS:
            samples[key] = archive[f"phi_{key}"].copy()
            with np.load(data/f"moments_phi_{key}.npz") as saved:
                results[key] = (saved["mean"].copy(), saved["centred_moments"].copy(), 0)
    make_figures(results, samples, summary["sample_n"], figures)
    print("Figures regenerated from the saved data, n =", summary["sample_n"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=2000)
    parser.add_argument("--sample-n", type=int, default=1000)
    parser.add_argument("--samples", type=int, default=12000)
    parser.add_argument("--figures-only", action="store_true",
                        help="Redraw the figures from the saved data without recomputing the moments or samples.")
    args = parser.parse_args()
    if not (3 <= args.sample_n <= args.max_n) or args.samples < 100:
        parser.error("Require 3 <= sample-n <= max-n and at least 100 samples.")
    data = ROOT/"data"
    figures = ROOT/"figures"
    data.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    if args.figures_only:
        redraw_saved(data, figures)
        return
    check = exhaustive_check()
    print("Exhaustive rational verification:", check, flush=True)
    results, samples, rows = {}, {}, []
    summary = {"verification": check, "parameters": {}, "sample_n": args.sample_n,
               "repetitions": args.samples, "base_seed": SEED,
               "python": platform.python_version(), "numpy": np.__version__,
               "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
               "longdouble_precision": int(np.finfo(np.longdouble).precision)}
    table = []
    for index, key in enumerate(PARAMETERS):
        mean, moment, residual = finite_moments(key, args.max_n)
        results[key] = (mean, moment, residual)
        q = float(key)
        v = q*(1-q)/3
        for n in range(3, args.max_n+1):
            var = moment[2, n]
            row = {"phi": key, "n": n, "mean": float(mean[n]), "variance": float(var),
                   "variance_ratio": float(var/(v*n**3)),
                   "skewness": float(moment[3, n]/var**np.longdouble("1.5")),
                   "kurtosis": float(moment[4, n]/var**2),
                   "sixth_standardised_moment": float(moment[6, n]/var**3)}
            rows.append(row)
            if n in {500, args.max_n}:
                table.append(row)
        samples[key] = sample_costs(q, args.sample_n, args.samples, SEED+index)
        summary["parameters"][key] = {
            "seed": SEED+index, "max_relative_centring_residual": residual,
            "finite_mean_at_sample_n": float(mean[args.sample_n]),
            "finite_variance_at_sample_n": float(moment[2, args.sample_n]),
            "sample_mean": float(np.mean(samples[key])),
            "sample_variance_unbiased": float(np.var(samples[key], ddof=1)),
        }
        np.savez_compressed(data/f"moments_phi_{key}.npz", mean=np.asarray(mean, dtype=float),
                            centred_moments=np.asarray(moment, dtype=float))
        print("Completed phi="+key, summary["parameters"][key], flush=True)
    with (data/"finite_moments.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_latex_rows(table, data/"numerical_rows.tex")
    np.savez_compressed(data/"simulated_costs.npz", **{f"phi_{k}": v for k, v in samples.items()})
    (data/"validation.json").write_text(json.dumps(summary, indent=2)+"\n")
    make_figures(results, samples, args.sample_n, figures)
    print("All tables, data and figures regenerated in", ROOT, flush=True)


if __name__ == "__main__":
    main()
