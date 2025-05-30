import streamlit as st
import GEOparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

# Set page config
st.set_page_config(page_title="Genescope", layout="wide")

# Title - Centered
st.markdown("<h1 style='text-align: center;'> Gene Expression Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<hr>", unsafe_allow_html=True)

# Sidebar: About & Guide
with st.sidebar:
    st.title(" Navigation")

    st.markdown("##  About")
    st.info("""
    This app helps researchers compare gene expression levels between two sets of biological samples (e.g., healthy vs diseased).  
    It uses publicly available datasets from the NCBI GEO database.
    """)

    st.markdown("##  User Guide")
    st.markdown("""
    1. Enter a GEO Series ID (e.g., `GSE42872`)  
    2. Select at least 2 samples for each group  
    3. View statistical results  
    4. Explore the volcano plot  
    5. Download CSV of results
    """)

    st.markdown("##  Project Details")
    st.markdown("""
Title : Gene Expression Analyzer  
Track : HealthTech & Bioinformatics  
Hackathon : Codeholics Hack 4 Mini 2.0  
Team  : ZEN-PAL
    """)

    demo_mode = st.checkbox("Use Demo Dataset")
    geo_id = st.text_input("Enter GEO Series ID", placeholder="e.g., GSE42872",
                           help="Example GEO ID: GSE42872 (from https://www.ncbi.nlm.nih.gov/geo/)")

# Landing Page when no input
if not geo_id and not demo_mode:
    st.markdown("<h2 style='text-align: center;'> Gene Expression Analyzer</h2>", unsafe_allow_html=True)
    st.markdown("""
This project is built for **Codeholics Hack 4 Mini 2.0** under the **HealthTech & Bioinformatics** theme.  
Analyze gene expression data using public datasets from **NCBI GEO**.

###  What can you do?
- Compare healthy vs diseased samples
- Identify differentially expressed genes (DEGs)
- Visualize volcano plots
- Export analysis results

 Enter a **GEO Series ID** in the sidebar or enable **Demo Mode** to get started.
""")
    st.stop()

# Optional Project Overview
if st.checkbox(" Show Full Project Overview"):
    st.markdown("##  Project Overview")
    st.success("""
###  Project Title: Genescope

Hackathon: Codeholics Hack 4 Mini 2.0  
Track    : HealthTech & Bioinformatics  
Team Name: ZEN-PAL

####  Objective
Build a web-based tool to:
- Fetch gene expression data from NCBI GEO
- Let users define sample groups
- Perform differential expression analysis (T-tests)
- Visualize results with volcano plots
- Download CSV of significant genes

####  Tech Stack
- Python + Streamlit
- GEOparse (data)
- Pandas / NumPy (processing)
- SciPy (stats)
- Matplotlib + Seaborn (plots)
- GitHub (repo/deploy)

####  Resources
- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- [Streamlit Docs](https://docs.streamlit.io/)
- [GitHub Repo](https://github.com/yourusername/yourrepo)
    """)
    st.stop()

# Tabbed Layout for Workflow
st.markdown("###  Analysis Workflow")
st.markdown("""
1️ Load GEO Dataset →  
2️ Select Healthy & Diseased Samples →  
3️ Run T-Test →  
4️ View Volcano Plot & Export Results  
""")

if geo_id:
    with st.spinner(f" Fetching dataset `{geo_id}`..."):
        gse = GEOparse.get_GEO(geo=geo_id, destdir="./")
        samples = gse.gsms
        data_matrix = pd.DataFrame({k: v.table["VALUE"] for k, v in samples.items()})
        data_matrix.columns = list(samples.keys())

    st.success(f" Loaded {len(samples)} samples successfully!")

    st.markdown("###  Select Sample Groups for Comparison")
    group1 = st.multiselect(" Healthy Samples", options=list(data_matrix.columns), key="healthy",
                            help="Pick samples that represent the healthy control group")
    group2 = st.multiselect(" Diseased Samples", options=list(data_matrix.columns), key="diseased",
                            help="Pick samples that represent the disease condition")

    if group1 and group2:
        if len(group1) < 2 or len(group2) < 2:
            st.warning(" Please select **at least 2 samples in each group** for meaningful statistical comparison.")
            st.stop()

        st.info(" Performing statistical analysis...")

        # Perform T-test and logFC
        t_stat, p_vals = ttest_ind(data_matrix[group1], data_matrix[group2], axis=1, nan_policy='omit')
        log_fc = np.log2(data_matrix[group2].mean(axis=1) + 1) - np.log2(data_matrix[group1].mean(axis=1) + 1)

        probe_ids = data_matrix.index.tolist()

        results = pd.DataFrame({
            "Gene": probe_ids,
            "logFC": log_fc,
            "p-value": p_vals,
            "-log10(p-value)": -np.log10(p_vals)
        }).reset_index(drop=True)

        st.success(" Differential expression analysis complete!")

        # Identify top 10 DEGs
        top_genes = results.sort_values("p-value").head(10)

        # Columns layout
        left_col, right_col = st.columns([1, 1.2], gap="large")

        with left_col:
            st.markdown("### Statistical Summary")
            st.dataframe(results.head(10), use_container_width=True)

            csv = results.to_csv(index=False).encode('utf-8')
            st.download_button(
                " Download Results as CSV",
                data=csv,
                file_name="deg_results.csv",
                mime="text/csv"
            )

        with right_col:
            st.markdown("###  Volcano Plot")
            fig, ax = plt.subplots(figsize=(8, 6))
            plot = sns.scatterplot(
                data=results,
                x="logFC",
                y="-log10(p-value)",
                hue=results["p-value"] < 0.05,
                palette={True: "red", False: "gray"},
                ax=ax
            )

            # Add gene labels for top 10
            for _, row in top_genes.iterrows():
                ax.text(row["logFC"], row["-log10(p-value)"], row["Gene"], fontsize=7, alpha=0.8)

            ax.set_title("Volcano Plot")
            ax.set_xlabel("log2 Fold Change")
            ax.set_ylabel("-log10 p-value")
            ax.axhline(y=-np.log10(0.05), color='blue', linestyle='--', linewidth=1)
            ax.grid(True)
            st.pyplot(fig)
else:
    st.info("ℹ Please enter a GEO Series ID in the sidebar to begin.")
