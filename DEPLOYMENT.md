# Deploy the interactive project

Repository: https://github.com/Fidy-Heritiana/Quicksort-Mallows

## 1. Upload the project

Extract the prepared archive. Upload the **contents** of `Quicksort-Mallows/`
to the repository root on branch `main`, replacing the initial README.
Do not upload the ZIP itself or put the whole project inside another
`Quicksort-Mallows/` directory.

The root must contain `streamlit_app.py`, `explore.py`, `reproduce.py`,
`requirements.txt`, `demo.ipynb`, and the `data/` and `figures/` directories.
Include `.github/`, `.streamlit/` and `.gitignore`. In a Linux file manager,
Ctrl+H reveals names beginning with a dot; in macOS Finder use Command+Shift+.

Suggested commit message: `Add reproducible experiments and interactive demos`.

## 2. Run the checks

Open **Actions** and select **Numerical and application checks**. A successful
run confirms that the finite checks and application tests passed in GitHub's
environment. If no workflow appears, check that `.github/workflows/ci.yml` was
uploaded to that exact path. Open the failed step's log if a run fails.

## 3. Deploy Streamlit

Open https://share.streamlit.io/ and connect the GitHub account that owns the
repository. Choose **Create app**, then the option for an existing app.

| Field | Value |
| --- | --- |
| Repository | `Fidy-Heritiana/Quicksort-Mallows` |
| Branch | `main` |
| Main file path | `streamlit_app.py` |
| Python version, under Advanced settings | `3.12` |
| Custom subdomain | Optional; choose an available name |

No API keys or application secrets are required. Click **Deploy**. Installation
and the first launch may take a few minutes. Copy the actual application URL
from the successful deployment page; a proposed subdomain is not guaranteed to
be available. Add the actual URL to the README and **About → Website**.

## 4. Check the public experience

- Manuscript data: select each parameter and inspect both plots.
- New simulation: run n = 200, φ = 0.5, 2,000 samples, then change the seed.
- Moment convergence: compare parameters and switch quantities.
- Download costs and check that the CSV records φ, n and the seed.
- Open Colab and run the notebook from top to bottom.
- Test Binder separately; its first environment build can be slow or temporarily
  unavailable depending on the public service.

## 5. Complete research metadata

Review the software authors in `CITATION.cff` and choose a project license
before announcing a release. Add manuscript citation details when available.
Add an archival DOI only after an archive has actually been created.

## Official documentation

- [Upload files](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository)
- [Deploy on Streamlit](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy)
- [Application dependencies](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/app-dependencies)
- [Binder](https://mybinder.readthedocs.io/en/latest/introduction.html)
