# Quicksort under the Mallows distribution

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://quicksort-mallows-ak6ykvfg6rf4taeffvvgbc.streamlit.app/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Fidy-Heritiana/Quicksort-Mallows/blob/main/demo.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Fidy-Heritiana/Quicksort-Mallows/main?labpath=demo.ipynb)
[![Checks](https://github.com/Fidy-Heritiana/Quicksort-Mallows/actions/workflows/ci.yml/badge.svg)](https://github.com/Fidy-Heritiana/Quicksort-Mallows/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24186089.svg)](https://doi.org/10.5281/zenodo.24186089).

**An analytic study of Quicksort comparisons: finite recurrences, moment
asymptotics, normal fluctuations, and reproducible experiments.**

This repository accompanies research on first-pivot Quicksort with a
Mallows-distributed input. The analytic approach uses generating functions and
moment recurrences. The Python code computes finite-size quantities, checks
identities using exact rational arithmetic, and illustrates the asymptotic
results with Monte Carlo experiments.

![Normalised manuscript samples](figures/normal_histograms.png)

## Explore interactively

- **Colab / Binder:** click a badge above and run `demo.ipynb` from top to bottom.
  The notebook provides parameter controls and figures. GitHub displays the
  notebook, while Colab and Binder execute its code.
- **Streamlit:** [open the live application](https://quicksort-mallows-ak6ykvfg6rf4taeffvvgbc.streamlit.app/)
  to explore the manuscript data, run simulations, inspect moment convergence,
  and download results directly in your browser. No installation is required.
  To run the application locally, use
  `python -m streamlit run streamlit_app.py`; see
  [the deployment instructions](DEPLOYMENT.md) for details.
  
| Application page | What you can explore |
| --- | --- |
| Manuscript data | Original samples for φ = 0.25, 0.5 and 0.8, histograms, normal Q–Q plots and CSV downloads. |
| New simulation | Choose φ, n, the number of simulations and the seed; compare empirical statistics with finite-recurrence calculations. |
| Moment convergence | Finite means, variances and standardised centred moments through n = 2,000. |
| Model and verification | The model, moment conventions and bundled exact-check reports. |

Plots support hovering and zooming. The app and notebook do not overwrite the
original data. Interactive runs use 0.01 ≤ φ ≤ 0.99, 3 ≤ n ≤ 1,000 and
500–20,000 simulations to keep the shared interface responsive. Notebook sliders
use smaller upper bounds; their computational helper accepts the full range.
The command-line script supports larger computations.

## Mathematical setting

For a permutation π of size n,

$$\mathbb P_\phi(\pi)=\frac{\phi^{\operatorname{inv}(\pi)}}{Z_n(\phi)},
\qquad 0<\phi<1.$$

The algorithm chooses the **first element as pivot** and uses **stable
partition**. Let Cₙ count key comparisons, with C₀ = C₁ = 0. Conditional on the
left-subproblem size J, the two subproblem costs are independent and

$$C_n\overset d=n-1+C_J+C'_{n-1-J},\qquad
\mathbb P(J=j)=\frac{(1-\phi)\phi^j}{1-\phi^n},\quad 0\le j<n.$$

The simulator samples these recursive splits directly. It has the same cost
law as running the specified algorithm on a Mallows permutation.
The finite expectation is

$$\mu_n=\sum_{k=1}^{n-1}(n-k)
\frac{(1-\phi)(1+\phi^k)}{1-\phi^{k+1}}.$$

The asymptotic results studied in the manuscript, for **fixed** φ ∈ (0,1), are

$$\mu_n=\frac{1-\phi}{2}n^2+O(n),\qquad
\operatorname{Var}(C_n)=\frac{\phi(1-\phi)}{3}n^3+O(n^2),$$

$$\frac{C_n-\mu_n}{\sqrt{\operatorname{Var}(C_n)}}
\xrightarrow{d}\mathcal N(0,1).$$

The uniform case φ = 1 belongs to a different asymptotic regime. At φ = 0 the
input is increasing and the cost is deterministic, so variance normalisation
is undefined. Neither endpoint is included in this demonstration.

### Moment conventions

| Quantity | Definition |
| --- | --- |
| Factorial moment fᵣ,ₙ | E[Cₙ(Cₙ − 1)⋯(Cₙ − r + 1)] |
| Ordinary moment mᵣ,ₙ | E[Cₙʳ] |
| Centred moment Mᵣ(n) | E[(Cₙ − μₙ)ʳ] |

Differentiating the probability generating function with respect to u,
evaluated at u = 1, gives **factorial moments**:

$$\left.\partial_u^r\sum_{n\ge0}\mathbb E[u^{C_n}]z^n\right|_{u=1}
=\sum_{n\ge0}f_{r,n}z^n.$$

In particular, Var(Cₙ) = f₂,ₙ + f₁,ₙ − f₁,ₙ². The numerical arrays contain
**centred ordinary moments** Mᵣ(n). Kurtosis in the tables tends to 3; it is
ordinary kurtosis. The sixth standardised centred moment tends to 15.

## Install and run locally

Python **3.12** is the reference environment.

```bash
git clone https://github.com/Fidy-Heritiana/Quicksort-Mallows.git
cd Quicksort-Mallows
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on Linux/macOS or
`.venv\Scripts\Activate.ps1` in Windows PowerShell. Then run:

```bash
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py
```

For the numerical scripts alone, `requirements-core.txt` is sufficient.
For a local notebook, install `requirements-notebook.txt` and Jupyter Notebook,
then open `demo.ipynb`. Its setup cell installs the notebook requirements.

## Reproduce the manuscript figures

```bash
python reproduce.py --figures-only
```

This reads the saved data and redraws `figures/variance_and_qq.{pdf,png}` and
`figures/normal_histograms.{pdf,png}`. The samples and normalisation match the
manuscript; the labels are English.

![Finite variance and normal Q–Q plot](figures/variance_and_qq.png)

To recompute all tables, samples and figures, run `python reproduce.py`.

| Setting | Original value |
| --- | --- |
| Mallows parameters | 0.25, 0.5, 0.8 |
| Maximum deterministic size | 2,000 |
| Maximum centred-moment order | 6 |
| Simulated input size | 1,000 |
| Samples per parameter | 12,000 |
| Seeds, in parameter order | 20260921, 20260922, 20260923 |

The full command replaces its generated files in `data/` and `figures/`.
For a smaller trial, use a separate copy of the repository:

```bash
python reproduce.py --max-n 60 --sample-n 30 --samples 500
```

The data-integrity check records the original manuscript dataset in
`data/SHA256SUMS`. Deliberately replacing that dataset requires updating the
checksums and the recorded experiment parameters.

## Verify the calculations

```bash
python verify_identities.py
python verify_details.py
python -m unittest discover -s tests -v
```

The first two commands use only the standard library. They check distributions
against all 5,913 permutations of sizes 1 through 7 and verify recurrences
through size 10 and moment order 6 for φ = 1/3, 1/2 and 2/3. They also check
moment conversions, generating-function identities, the fourth-moment source,
sampling intervals, exact centring and reciprocal-parameter symmetry.

The application tests check navigation, submitted parameters, data integrity,
and finite-size values independently known for n = 3. GitHub Actions runs the
checks and redraws the manuscript figures on pushes and pull requests.

Finite calculations and simulations complement the analytic proofs; they do
not establish the asymptotic statements.

## Repository contents

| Path | Purpose |
| --- | --- |
| `reproduce.py` | Original English reproduction script: full finite recurrences and recursive-law simulation. |
| `verify_identities.py`, `verify_details.py` | Original exact rational checks. |
| `explore.py` | Shared read-only helpers, data exports and interactive plots. |
| `streamlit_app.py` | Browser application. |
| `demo.ipynb` | Colab/Binder/Jupyter demonstration. |
| `data/`, `figures/` | Original data, reports and English manuscript figures. |
| `tests/`, `.github/workflows/ci.yml` | Numerical and application checks. |
| `requirements*.txt` | Pinned dependencies. |
| `runtime.txt` | Python 3.12 for Binder. |
| `CITATION.cff` | Software citation metadata. |
| `DEPLOYMENT.md` | Deployment instructions. |
| `VALIDATION.md` | Completed local checks and checks remaining after upload. |

## Reproducibility and citation

The three original Python scripts and numerical files are preserved from the
English code package. Finite recurrences use `numpy.longdouble`, whose effective
precision depends on the operating system. Saved arrays and CSV files use
double precision; simulated costs are integers.

The manuscript run used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and
Matplotlib 3.10.8. Seeds and finite reference values are recorded in
`data/validation.json`. The interactive layer uses the same numerical functions.

Use `CITATION.cff` to cite the software, and record the commit identifier when
reporting a reproduced experiment. Manuscript citation details and any archival
DOI can be added when available. The manuscript itself is not included here.
