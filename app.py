import streamlit as st
import GEOparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(page_title="GeneExpression-HealthHack", layout="wide")
st.title("🧬 Gene Expression Analyzer – HealthTech Hackathon")

geo_id = st.text_input("🔍 Enter GEO Series ID (e.g., GSE42872):")

if geo_id:
    with st.spinner("📥 Downloading dataset..."):
        try:
            gse = GEOparse.get_GEO(geo=geo_id, destdir="./")
            samples = gse.gsms
            data_matrix = pd.DataFrame({k: v.table["VALUE"] for k, v in samples.items()})
            data_matrix.columns = list(samples.keys())
            st.success(f"✅ Loaded {len(samples)} samples.")
        except Exception as e:
            st.error(f"❌ Error loading GEO data: {e}")
            st.stop()

    st.subheader("🧬 Select Sample Groups to Compare")
    group1 = st.multiselect("Healthy Group", options=data_matrix.columns)
    group2 = st.multiselect("Diseased Group", options=data_matrix.columns)

    if group1 and group2:
        if set(group1).intersection(set(group2)):
            st.error("❌ A sample can't be in both groups.")
            st.stop()

        try:
            # Run T-test
            t_stat, p_vals = ttest_ind(data_matrix[group1], data_matrix[group2], axis=1, nan_policy='omit')
            log_fc = np.log2(data_matrix[group2].mean(axis=1) + 1) - np.log2(data_matrix[group1].mean(axis=1) + 1)

            results = pd.DataFrame({
                "logFC": log_fc,
                "p-value": p_vals,
                "-log10(p-value)": -np.log10(p_vals)
            })

            # Clean invalid rows
            results_clean = results.replace([np.inf, -np.inf], np.nan).dropna()

            st.subheader("📋 Result Summary (Top 5 Genes)")
            st.dataframe(results_clean.head())

            st.write(f"📊 Total Genes: {len(results)} | Plottable: {len(results_clean)}")

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

                # Download results
                csv = results_clean.to_csv(index=False).encode("utf-8")
                st.download_button("📥 Download Results", csv, "deg_results.csv", "text/csv")
            else:
                st.warning("⚠️ No valid gene results found to plot.")
        except Exception as e:
            st.error(f"❌ Error in statistical analysis: {e}")



