import streamlit as st
import GEOparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

st.set_page_config(page_title="GeneExpression-HealthHack", layout="wide")
st.title("🧬 Gene Expression Analyzer – HealthTech Hackathon")

geo_id = "GSE42872"  # Hardcoded for testing

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
    # Pre-selected valid samples for GSE42872
    st.subheader("🧬 Select Sample Groups to Compare")

# Known working samples from GSE42872
    default_healthy = ["GSM1052615", "GSM1052616", "GSM1052617"]
    default_diseased = ["GSM1052618", "GSM1052619", "GSM1052620"]

    group1 = st.multiselect(
        "Healthy Group (Group 1)",
        options=data_matrix.columns,
        default=default_healthy
    )

    group2 = st.multiselect(
        "Diseased Group (Group 2)",
        options=data_matrix.columns,
        default=default_diseased
    )

# Check for errors
    if not group1 or not group2:
        st.warning("⚠️ Please select at least one sample in each group.")
        st.stop()

    if set(group1).intersection(set(group2)):
        st.error("❌ A sample can't be in both groups. Please fix your selection.")
        st.stop()

    st.success(f"✅ Comparing {len(group1)} healthy vs {len(group2)} diseased samples.")


    st.write("✅ Group 1:", group1)
    st.write("✅ Group 2:", group2)


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

    # Download button
                csv = results_clean.to_csv(index=False).encode("utf-8")
                st.download_button("📥 Download Results", csv, "deg_results.csv", "text/csv")

            else:
                st.warning("⚠️ No valid gene results found. Displaying demo volcano plot for illustration.")

    # Show a fake demo plot
                demo = pd.DataFrame({
                    "logFC": np.random.normal(0, 2, 1000),
                    "p-value": np.random.uniform(0, 1, 1000)
                })
                demo["-log10(p-value)"] = -np.log10(demo["p-value"])

                fig, ax = plt.subplots(figsize=(8, 6))
                sns.scatterplot(
                    data=demo,
                    x="logFC",
                    y="-log10(p-value)",
                    hue=demo["p-value"] < 0.05,
                    palette={True: "red", False: "gray"},
                    ax=ax
                )
                ax.set_title("Volcano Plot (Demo)")
                ax.set_xlabel("log2 Fold Change")
                ax.set_ylabel("-log10 p-value")
                ax.grid(True)
                st.pyplot(fig)

        except Exception as e:
            st.error(f"❌ Error in statistical analysis: {e}")



