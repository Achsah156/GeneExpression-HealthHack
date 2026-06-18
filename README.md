<div align="center">

<img src="https://img.shields.io/badge/HealthTech-%2300C49A?style=for-the-badge&logo=dna&logoColor=white" />
<img src="https://img.shields.io/badge/Bioinformatics-%230A66C2?style=for-the-badge&logo=azuredevops&logoColor=white" />
<img src="https://img.shields.io/badge/Built%20with%20Python-3776AB?style=for-the-badge&logo=python&logoColor=white" />
<img src="https://img.shields.io/badge/Streamlit-FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white" />

# 🧬 Genescope — Gene Expression Analyzer

**Decode the Language of Genes. Instantly.**

_Codeholics Hack 4 Mini 2.0 · HealthTech & Bioinformatics Track_

[![Live App](https://img.shields.io/badge/🚀%20Live%20App-Try%20Now-brightgreen?style=for-the-badge)](https://geneexpression-healthhack-rpqrcgnmw5ustscwaj6tba.streamlit.app)
[![GitHub Stars](https://img.shields.io/github/stars/yourusername/genescope?style=for-the-badge&logo=github)](https://github.com/yourusername/genescope)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)](LICENSE)

</div>

---

## 🔬 What is Genescope?

**Genescope** is a researcher-friendly web tool for **differential gene expression (DEG) analysis** — no coding required. It connects directly to the [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) database, letting scientists go from raw dataset ID to statistically validated, publication-ready results in minutes.

Built for the **HealthTech & Bioinformatics** track of **Codeholics Hack 4 Mini 2.0**, Genescope bridges the gap between raw genomic data and meaningful biological insight.

> 💡 **No GEO dataset yet?** Enable **Demo Mode** to explore all features with a built-in test dataset — no setup needed.

---

## ✨ Features

| Feature | Description |
|---|---|
| 📂 **GEO Dataset Loader** | Load any NCBI GEO dataset instantly via its GSE ID |
| ⚖️ **Group Comparison** | Compare healthy vs. diseased sample groups with ease |
| 📊 **Statistical Analysis** | Auto-compute log₂ fold change and p-values |
| 🌋 **Interactive Volcano Plot** | Visualize DEGs with a clean, interactive plot |
| 💾 **CSV Export** | Download your differentially expressed genes for downstream analysis |
| 🧠 **Biological Insights Panel** | Get contextual biological interpretation of results |
| 📈 **DEG Summary Charts** | At-a-glance charts summarizing expression patterns |
| 📱 **Responsive Design** | Fully functional on mobile, tablet, and desktop |
| 🎓 **In-App Guide** | Built-in tutorial so anyone can get started quickly |
| 🎭 **Demo Mode** | Explore the full experience without uploading any data |

---

## 🖥️ Live Demo

👉 **[Launch Genescope →](https://geneexpression-healthhack-rpqrcgnmw5ustscwaj6tba.streamlit.app)**

> Try it with any valid **GSE ID** from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/), or enable **Demo Mode** from the sidebar for an instant walkthrough.

---

## 🛠️ Tech Stack

```
┌─────────────────────────────────────────────────────────────────┐
│                        GENESCOPE STACK                          │
├──────────────┬──────────────────────────────────────────────────┤
│  Frontend    │  Streamlit · HTML/CSS · Seaborn · Matplotlib     │
│  Backend     │  Python · GEOparse · Pandas · NumPy · SciPy      │
│  Deployment  │  Streamlit Cloud · GitHub                        │
└──────────────┴──────────────────────────────────────────────────┘
```

---

## ⚙️ Getting Started Locally

### Prerequisites

- Python 3.8+
- pip

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/genescope.git
cd genescope

# 2. Install dependencies
pip install -r requirements.txt

# 3. Launch the app
streamlit run app.py
```

The app will open automatically at `http://localhost:8501`.

### Quick Start

1. Enter a **GSE ID** (e.g., `GSE12345`) in the sidebar
2. Select your **healthy** and **diseased** sample groups
3. Click **Analyze** to generate your volcano plot and DEG table
4. **Download** results as CSV or explore the Biological Insights panel

> ✅ Or just check **Demo Mode** in the sidebar to skip straight to the results!

---

## 📁 Folder Structure

```
genescope/
├── app.py               # Main Streamlit application
├── requirements.txt     # Python dependencies
├── README.md            # Project documentation
└── data/                # (Optional) Cached GEO datasets
```

---

## 📸 Screenshots

> _Volcano plot · DEG summary charts · Biological Insights panel_

<!-- Add screenshots here once available -->
<!-- ![Volcano Plot](assets/volcano_plot.png) -->
<!-- ![DEG Summary](assets/deg_summary.png) -->

---

## 🧑‍💻 Team ZEN-PAL

| Name | Role |
|---|---|
| **Achsah Grace Karri** | Team Member |
| **Srishanth** | Team Member |
| **Laxmi** | Team Member |
| **RamPrasad** | Team Member |

Built with ❤️ for **Codeholics Hack 4 Mini 2.0** · Theme: _HealthTech & Bioinformatics_

---

## 📚 Resources & References

- [NCBI GEO Database](https://www.ncbi.nlm.nih.gov/geo/) — Gene Expression Omnibus
- [Streamlit Documentation](https://docs.streamlit.io/) — Web app framework
- [GEOparse Library](https://geoparse.readthedocs.io/) — Python GEO data parser
- [Seaborn](https://seaborn.pydata.org/) — Statistical data visualization

---

## 🤝 Contributing

Contributions, issues, and feature requests are welcome!

1. Fork the repository
2. Create your feature branch: `git checkout -b feature/amazing-feature`
3. Commit your changes: `git commit -m 'Add amazing feature'`
4. Push to the branch: `git push origin feature/amazing-feature`
5. Open a Pull Request

---

## 📄 License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for more information.

---

<div align="center">

### ⭐ If Genescope helped you, give it a star!

_Star · Fork · Share with researchers and developers_

**Made with 🧬 by Team ZEN-PAL**

</div>
