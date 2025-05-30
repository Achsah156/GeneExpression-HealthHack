import streamlit as st
import GEOparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(page_title="GeneExpression-HealthHack", layout="wide")
st.title("🧬 Gene Expression Analyzer")

geo_id = st.text_input("Enter GEO Series ID (e.g., GSE42872):")

if geo_id:
    with st.spinner(f"🔍 Fetching dataset {geo_id}..."):
        try:
            gse = GEOparse.get_GEO(geo=geo_id, destdir="./")
            samples = gse.gsms
            data_matrix = pd.DataFrame({k: v.table["VALUE"] for k, v in samples.items()})
            data_matrix.columns = list(samples.keys())
            st.success(f"✅ Loaded {len(samples)} samples.")
        except Exception as e:
            st.error(f"❌ Failed to load dataset: {e}")
            st.stop()

    st.subheader("📌 Select Sample Groups for Comparison")
    group1 = st.multiselect("Healthy Samples (Group 1)", options=list(data_matrix.columns))
    group2 = st.multiselect("Diseased Samples (Group 2)", options=list(data_matrix.columns))

    if group1 and group2:
        if len(set(group1).intersection(set(group2))) > 0:
            st.error("❌ A sample can't be in both groups. Please select different ones.")
            st.stop()

        try:
            st.info("🧪 Performing statistical analysis...")
            t_stat, p_vals = ttest_ind(data_matrix[group1], data_matrix[group2], axis=1)
            log_fc = np.log2(data_matrix[group2].mean(axis=1) + 1) - np.log2(data_matrix[group1].mean(axis=1) + 1)

            results = pd.DataFrame({
                "logFC": log_fc,
                "p-value": p_vals,
                "-log10(p-value)": -np.log10(p_vals)
            })

            # Clean data: remove NaNs and infs
            results_clean = results.replace([np.inf, -np.inf], np.nan).dropna()

            st.write(f"🔍 Total genes: {len(results)} | Valid genes after cleaning: {len(results_clean)}")

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

                # CSV download
                csv = results_clean.to_csv(index=False).encode("utf-8")
                st.download_button("📥 Download DEG Results as CSV", csv, "deg_results.csv", "text/csv")
            else:
                st.warning("⚠️ No valid gene data to display. Try different samples.")

        except Exception as e:
            st.error(f"❌ Error during analysis: {e}")

