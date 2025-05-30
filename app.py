import streamlit as st
import GEOparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(page_title="GeneExpression-HealthHack", layout="wide")
st.title(" Gene Expression Analyzer")

# User input
geo_id = st.text_input("Enter GEO Series ID (e.g., GSE42872):")

if geo_id:
    with st.spinner(f"Fetching dataset {geo_id}..."):
        gse = GEOparse.get_GEO(geo=geo_id, destdir="./")
        samples = gse.gsms
        data_matrix = pd.DataFrame({k: v.table["VALUE"] for k, v in samples.items()})
        data_matrix.columns = list(samples.keys())

    st.success(f" Loaded {len(samples)} samples.")

    st.subheader("Select Sample Groups for Comparison")
    group1 = st.multiselect("Healthy Samples", options=list(data_matrix.columns))
    group2 = st.multiselect("Diseased Samples", options=list(data_matrix.columns))

    if group1 and group2:
        st.info(" Running T-test...")
        t_stat, p_vals = ttest_ind(data_matrix[group1], data_matrix[group2], axis=1)
        log_fc = np.log2(data_matrix[group2].mean(axis=1) + 1) - np.log2(data_matrix[group1].mean(axis=1) + 1)

        results = pd.DataFrame({
            "logFC": log_fc,
            "p-value": p_vals,
            "-log10(p-value)": -np.log10(p_vals)
        })

        st.success(" Analysis complete!")

        # Volcano plot
        st.subheader("Volcano Plot")
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.scatterplot(data=results, x="logFC", y="-log10(p-value)", hue=results["p-value"] < 0.05, palette={True: "red", False: "gray"}, ax=ax)
        ax.set_title("Volcano Plot")
        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10 p-value")
        ax.grid(True)
        st.pyplot(fig)

        # Download button
        csv = results.to_csv(index=False).encode('utf-8')
        st.download_button(
            " Download DEG Results as CSV",
            data=csv,
            file_name="deg_results.csv",
            mime="text/csv"
        )
