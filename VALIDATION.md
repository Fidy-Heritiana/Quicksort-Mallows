# Validation of the prepared GitHub package

Validation date: 2026-11-12.

## Completed locally

- All three original English Python scripts and all nine original data files
  are byte-for-byte identical to the earlier English code package.
- `verify_identities.py`: passed with exact rational arithmetic.
- `verify_details.py`: passed with exact rational arithmetic.
- `python -m unittest discover -s tests -v`: all five tests passed. These cover
  the four application pages, parameter submission, saved data, reproducibility,
  input validation and the independently known n = 3 distribution.
- The two parameter-range endpoints, phi = 0.01 and phi = 0.99, were checked at
  n = 1,000 with positive finite variances. A single-size selection on the
  convergence page was also checked successfully.
- The notebook format was validated, and all nine Python cells were executed
  sequentially in-process with IPython. This verifies the Python code and widget
  construction, but is not a hosted Colab/Binder browser test.
- The installed dependency set passed `pip check`.

## To complete after upload

- Run GitHub Actions in the actual repository.
- Deploy Streamlit Community Cloud and inspect the browser rendering, plot
  controls and downloads. Automated application tests passed locally; a visual
  browser inspection could not be completed in the preparation environment.
- Open the notebook on Colab and Binder and test the widget buttons.

No public Streamlit URL or archival DOI is claimed before those services have
actually been configured. Numerical experiments and finite checks complement
the manuscript's analytic proofs.
