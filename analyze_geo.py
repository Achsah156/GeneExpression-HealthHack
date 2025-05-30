import GEOparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

print("📥 Downloading dataset GSE42872...")
gse = GEOparse.get_GEO("GSE42872", destdir="./")

# Extract expression matrix
samples = gse.gsms
print(f"✅ Loaded {len(samples)} samples.")

# Create data matrix from sample VALUEs
data_matrix = pd.DataFrame({k: v.table["VALUE"] for k, v in samples.items()})
data_matrix.columns = list(samples.keys())

# Simulate healthy vs diseased (first 3 vs next 3)
group1 = list(data_matrix.columns)[:3]
group2 = list(data_matrix.columns)[3:6]

print(f"🧪 Running T-test on groups:\nHealthy: {group1}\nDiseased: {group2}")

# Calculate differential expression
t_stat, p_vals = ttest_ind(data_matrix[group1], data_matrix[group2], axis=1)
log_fc = np.log2(data_matrix[group2].mean(axis=1) + 1) - np.log2(data_matrix[group1].mean(axis=1) + 1)

# Create result table
results = pd.DataFrame({
    "logFC": log_fc,
    "p-value": p_vals,
    "-log10(p-value)": -np.log10(p_vals)
})

# Save results
results.to_csv("deg_results_new.csv")
print("✅ DEG results saved to deg_results.csv")

# Volcano plot
plt.figure(figsize=(8, 6))
sns.scatterplot(data=results, x="logFC", y="-log10(p-value)", hue=results["p-value"] < 0.05, palette={True: "red", False: "gray"})
plt.title("Volcano Plot – GSE42872")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 p-value")
plt.grid(True)
plt.tight_layout()
plt.savefig("volcano_plot.png")
plt.show()
