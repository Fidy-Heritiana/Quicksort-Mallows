"""Interactive companion to the Quicksort–Mallows research code."""
import json

import numpy as np
import pandas as pd
import streamlit as st

from explore import (
    METRICS, ROOT, SAVED_PARAMETERS, convergence_figure, experiment_metadata,
    load_manuscript, moment_table, normal_figures, run_experiment,
    sample_table, summary_table,
)

REPO_URL = "https://github.com/Fidy-Heritiana/Quicksort-Mallows"

st.set_page_config(page_title="Quicksort · Mallows", page_icon="↗", layout="wide")


@st.cache_data(max_entries=3, show_spinner=False)
def saved_experiment(phi):
    return load_manuscript(phi)


@st.cache_data(max_entries=12, show_spinner=False)
def new_experiment(phi, n, repetitions, seed):
    return run_experiment(phi, n, repetitions, seed)


@st.cache_data(max_entries=1, show_spinner=False)
def all_saved_moments():
    return pd.concat([moment_table(load_manuscript(q)) for q in SAVED_PARAMETERS], ignore_index=True)


def show_experiment(result, prefix):
    metadata = experiment_metadata(result)
    st.caption(
        f"{result['source']} · φ = {result['phi']} · n = {result['n']:,} · "
        f"{result['repetitions']:,} samples · seed = {result['seed']}"
    )
    a, b, c = st.columns(3)
    a.metric("Finite mean μₙ", f"{metadata['finite_mean']:,.3f}")
    b.metric("Finite variance σₙ²", f"{metadata['finite_variance']:,.3f}")
    c.metric("Simulated mean", f"{np.mean(result['costs']):,.3f}")
    st.write("Costs are centred and scaled using the finite theoretical mean and variance.")
    hist, qq = normal_figures(result)
    left, right = st.columns(2)
    with left:
        st.plotly_chart(hist, width="stretch", key=prefix + "_hist")
    with right:
        st.plotly_chart(qq, width="stretch", key=prefix + "_qq")
    st.dataframe(summary_table(result), hide_index=True, width="stretch")
    st.caption("The empirical variance uses the unbiased estimator (ddof = 1). "
               "At small n the cost distribution is discrete; the normal approximation is asymptotic.")
    csv = sample_table(result).to_csv(index=False).encode("utf-8")
    col1, col2 = st.columns(2)
    col1.download_button("Download simulated costs (CSV)", csv,
                         file_name=f"costs_phi_{result['phi']}_n_{result['n']}_seed_{result['seed']}.csv",
                         mime="text/csv", key=prefix + "_csv")
    col2.download_button("Download experiment parameters (JSON)",
                         json.dumps(metadata, indent=2), file_name="experiment_parameters.json",
                         mime="application/json", key=prefix + "_json")


st.sidebar.title("Quicksort · Mallows")
page = st.sidebar.radio("Explore", ["Manuscript data", "New simulation",
                                    "Moment convergence", "Model and verification"], key="page")
st.sidebar.caption("First pivot · stable partition · fixed 0 < φ < 1")
st.sidebar.markdown(f"[Source code]({REPO_URL}) · [Notebook]({REPO_URL}/blob/main/demo.ipynb)")

st.title("Quicksort under the Mallows distribution")
st.write("Explore how the input distribution changes the number of comparisons.")

if page == "Manuscript data":
    st.subheader("Explore the manuscript experiment")
    st.write("The original samples and finite-size calculations are included in the repository. "
             "Changing this selector loads a saved experiment.")
    phi = st.selectbox("Mallows parameter φ", SAVED_PARAMETERS, index=1, key="saved_phi")
    result = saved_experiment(phi)
    show_experiment(result, "saved")
    with st.expander("Original manuscript figures"):
        st.image(str(ROOT / "figures/variance_and_qq.png"),
                 caption="Original variance and Q–Q figure, with English labels.")
        st.image(str(ROOT / "figures/normal_histograms.png"),
                 caption="Original normalised histograms, with English labels.")

elif page == "New simulation":
    st.subheader("Run your own experiment")
    st.write("Choose parameters and press Run experiment. The finite moments are computed "
             "from the full recurrences, and independent costs are sampled from the recursive law.")
    with st.form("experiment_form"):
        left, right = st.columns(2)
        phi = left.slider("Mallows parameter φ", 0.01, 0.99, 0.50, 0.01, key="custom_phi")
        n = right.number_input("Input size n", min_value=3, max_value=1000,
                               value=200, step=1, key="custom_n")
        repetitions = left.slider("Number of simulations", 500, 20000, 2000, 500, key="custom_reps")
        seed = right.number_input("Random seed", min_value=0, max_value=2**32-1,
                                  value=20260922, step=1, key="custom_seed")
        submitted = st.form_submit_button("Run experiment", type="primary")
    st.caption("The range 0.01 ≤ φ ≤ 0.99 excludes the deterministic endpoint φ = 0 "
               "and the distinct uniform regime φ = 1. New simulations do not change the manuscript data.")
    if submitted:
        with st.spinner("Computing finite moments and drawing costs…"):
            result = new_experiment(f"{phi:.2f}", int(n), int(repetitions), int(seed))
        st.session_state["custom_result"] = result
    if "custom_result" in st.session_state:
        result = st.session_state["custom_result"]
        st.subheader("Last completed experiment")
        show_experiment(result, "custom")
        table = moment_table(result)
        st.download_button("Download finite moments for this experiment (CSV)",
                           table.to_csv(index=False), file_name="finite_moments_custom.csv",
                           mime="text/csv", key="custom_moments_csv")
    else:
        st.info("Press Run experiment to calculate and display your results.")

elif page == "Moment convergence":
    st.subheader("Finite moments and their asymptotic limits")
    table = all_saved_moments()
    selected = st.multiselect("Mallows parameters", SAVED_PARAMETERS,
                              default=list(SAVED_PARAMETERS), key="moment_phis")
    metric = st.selectbox("Quantity", list(METRICS),
                          format_func=lambda key: METRICS[key][0], index=1, key="moment_metric")
    n_min, n_max = st.slider("Input-size range", 3, 2000, (20, 2000), key="moment_range")
    selected_table = table[table["phi"].isin([float(q) for q in selected])
                           & table["n"].between(n_min, n_max)]
    if selected_table.empty:
        st.info("Select at least one Mallows parameter.")
    else:
        st.plotly_chart(convergence_figure(selected_table, metric), width="stretch")
        st.download_button("Download displayed moment data (CSV)",
                           selected_table.to_csv(index=False), file_name="finite_moments_selection.csv",
                           mime="text/csv")
        inspect_n = st.number_input("Inspect a size n", min_value=n_min, max_value=n_max,
                                     value=n_max, step=1, key=f"inspect_{n_min}_{n_max}")
        st.dataframe(selected_table[selected_table["n"] == inspect_n], hide_index=True,
                     width="stretch")
    st.caption("These are deterministic finite-recurrence calculations in floating-point arithmetic. "
               "The fourth standardised centred moment tends to 3 (ordinary kurtosis, not excess kurtosis).")

else:
    st.subheader("Model and interpretation")
    st.markdown("**Algorithm.** The first element is the pivot; stable partition preserves "
                "the relative order in each sublist. Cₙ counts key comparisons.")
    st.latex(r"\mathbb{P}_\phi(\pi)=\frac{\phi^{\operatorname{inv}(\pi)}}{Z_n(\phi)},\qquad 0<\phi<1.")
    st.latex(r"C_n\overset{d}=n-1+C_J+C'_{n-1-J},\quad "
             r"\mathbb{P}(J=j)=\frac{(1-\phi)\phi^j}{1-\phi^n}.")
    st.write("The two subproblem costs are independent conditional on J. The simulator draws "
             "these recursive splits directly instead of explicitly storing permutations.")
    st.markdown("**Finite expectation.**")
    st.latex(r"\mu_n=\sum_{k=1}^{n-1}(n-k)\frac{(1-\phi)(1+\phi^k)}{1-\phi^{k+1}}.")
    st.markdown("**Asymptotic results studied in the manuscript, for fixed φ.**")
    st.latex(r"\mu_n=\frac{1-\phi}{2}n^2+O(n),\qquad "
             r"\sigma_n^2=\frac{\phi(1-\phi)}{3}n^3+O(n^2).")
    st.latex(r"\frac{C_n-\mu_n}{\sigma_n}\xrightarrow{d}\mathcal{N}(0,1).")
    st.write("The plots illustrate finite-size behaviour. Simulations and finite checks "
             "complement the analytic proofs; they do not establish the limit theorem.")
    st.markdown("**Moment conventions.**")
    st.latex(r"f_{r,n}=\mathbb E[(C_n)_{\underline r}],\quad "
             r"m_{r,n}=\mathbb E[C_n^r],\quad M_r(n)=\mathbb E[(C_n-\mu_n)^r].")
    st.latex(r"\left.\partial_u^r\sum_{n\ge0}\mathbb E[u^{C_n}]z^n\right|_{u=1}"
             r"=\sum_{n\ge0}f_{r,n}z^n,\qquad \sigma_n^2=f_{2,n}+f_{1,n}-f_{1,n}^2.")
    st.write("The numerical arrays in this application store Mᵣ(n), the centred ordinary moments.")
    st.subheader("Exact finite-size checks")
    st.code("python verify_identities.py\npython verify_details.py", language="bash")
    st.write("The checks use rational arithmetic and direct enumeration of all permutations "
             "through size 7, together with recurrence checks through size 10 and order 6.")
    for name in ("verification_identities.json", "verification_details.json"):
        with st.expander(f"Bundled report: {name}"):
            st.json(json.loads((ROOT / "data" / name).read_text()))

st.divider()
st.caption("Analytic study · reproducible computations · interactive experiments")
