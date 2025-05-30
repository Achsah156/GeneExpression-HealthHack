import streamlit as st
import GEOparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
from streamlit_lottie import st_lottie
import requests

# Page config
st.set_page_config(page_title="GeneExpression-HealthHack", layout="wide")

# Lottie animation loader
def load_lottieurl(url):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()

# Load animations
lottie_bio = load_lottieurl("https://lottie.host/4930dbf7-4ae7-4ee9-bd8d-0613a1d07ba6/UhloOa8pGg.json")
lottie_plot = load_lottieurl("https://lottie.host/c84e52a1-d313-443a-9e13-d6f58d0db3cd/KLOmFiC2Jh.json")

# Sidebar navigation
st.sidebar.title("🔬 Navigation")
page = st.sidebar.radio("Go to:", ["Workspace", "User Guidelines", "About the Project"])

# ==========================================
# PAGE: ABOUT
# ==========================================
if page == "About the Project":
    st.title("🧬 GeneExpression-HealthHack")
    st.subheader("Simplifying bioinformatics for everyone.")
    st.markdown("""
    This web app helps researchers and students analyze gene expression datasets from NCBI GEO — without needing to write any code.  
    Simply enter a GEO Series ID, choose your sample groups, and get:
    - Volcano plots
    - Statistical results
    - CSV download of DEGs

    **Built by Team CodeHackers** for Codeholics Hack 4 Mini 2.0 💡
    """)
    st_lottie(lottie_bio, height=250)

# ==========================================
# PAGE: USER GUIDELINES
# ==========================================
elif page == "User Guidelines":
    st.title("📘 How to Use This App")

    with st.expander("Step-by-step Guide"):
        st.markdown("""
        1. **Enter a valid GEO Series ID** (e.g., `GSE42872`)  
        2. Wait for the dataset to load  
        3. **Select at least 2 samples** for Healthy and Diseased groups  
        4. Click to run the analysis  
        5. View the volcano plot  
        6. **Download results as CSV**  

        ⚠️ Avoid overlapping samples in both groups.  
        ⚠️ At least 2 samples in each group are required.
        """)

    st.image("volcano_placeholder.png", caption="Sample Volcano Plot", use_column_width=True)

# ==========================================
# PAGE: WORKSPACE
# ==========================================
elif page == "Workspace":
    st.title("🧪 Gene Expression Workspace")

    geo_id = st.text_input("Enter GEO Series ID (e.g., GSE42872):")

    if geo_id:
        with st.spinner(f"📥 Fetching dataset `{geo_id}`..."):
            try:
                gse = GEOparse.get_GEO(geo=geo_id, destdir="./")
                samples = gse.gsms
                data_matrix = pd.DataFrame({k: v.table["VALUE"] for k, v in samples.items()})
                data_matrix.columns = list(samples.keys())
                st.success(f"✅ Loaded {len(samples)} samples.")
            except Exception as e:
                st.error(f"❌ Could not load dataset: {e}")
                st.stop()

        st.subheader("🎯 Select Sample Groups")
        group1 = st.multiselect("Healthy Samples", options=data_matrix.columns)
        group2 = st.multiselect("Diseased Samples", options=data_matrix.columns)

        if group1 and group2:
            if set(group1).intersection(set(group2)):
                st.error("❌ A sample cannot be in both groups.")
                st.stop()

            if len(group1) < 2 or len(group2) < 2:
                st.warning("⚠️ Please select at least 2 samples in each group.")
                st.stop()

            with st.spinner("🧪 Running statistical analysis..."):
                t_stat, p_vals = ttest_ind(data_matrix[group1], data_matrix[group2], axis=1, nan_policy='omit')
                log_fc = np.log2(data_matrix[group2].mean(axis=1) + 1) - np.log2(data_matrix[group1].mean(axis=1) + 1)

                results = pd.DataFrame({
                    "logFC": log_fc,
                    "p-value": p_vals,
                    "-log10(p-value)": -np.log10(p_vals)
                })

                results_clean = results.replace([np.inf, -np.inf], np.nan).dropna()

            if len(results_clean) > 0:
                st.subheader("📈 Volcano Plot")
                fig, ax = plt.subplots(figsize=(8, 6))
                sns.scatterplot(
                    data=results_clean,
                    x="logFC",
                    y="-log10(p-value)",
                    hue=results_clean["p-value"] < 0.05,
                    palette={True: "red", False: "gray"},
                    ax=ax
                )
                ax.set_title("Volcano Plot")
                ax.set_xlabel("log2 Fold Change")
                ax.set_ylabel("-log10 p-value")
                ax.grid(True)
                st.pyplot(fig)

                st.markdown("📥 **Download Results**")
                csv = results_clean.to_csv(index=False).encode('utf-8')
                st.download_button("Download CSV", csv, "deg_results.csv", "text/csv")
            else:
                st.warning("⚠️ No valid gene results to plot. Showing sample volcano illustration.")
                st.image("volcano_placeholder.png", caption="Sample Volcano Plot", use_column_width=True)
                st_lottie(lottie_plot, height=300)

# ==========================================
# FOOTER
# ==========================================
st.markdown("---")
col1, col2 = st.columns([1, 3])
with col1:
    st_lottie(load_lottieurl("https://lottie.host/f80d3c6b-8fc5-4c40-804e-e23c377c11d7/nL1vYQWOK4.json"), height=60)
with col2:
    st.markdown("🔧 Built by **Team CodeHackers** | Codeholics Hack 4 Mini 2.0 | © 2025")
