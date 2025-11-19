# DGE-Analysis-of-GSE152641-Using-DESeq2-with-PCA-Heatmap-and-Volcano-Plot-Visualization
DGE Analysis of GSE152641 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization
# 🦠 DESeq2 DGE Pipeline: SARS-CoV-2 Host Response (GSE152641)

This R script automates a complete bioinformatics workflow for **Differential Gene Expression (DGE)** analysis of the **GSE152641** dataset, comparing **SARS-CoV-2 (COVID-19) infected patients** against **Healthy Controls**.

The pipeline uses the **`DESeq2`** package for count-based DGE, automatically handles data alignment and GEO file retrieval, and generates essential quality control (QC) and results visualizations.

## 🚀 Key Features

* **Automated Data Handling:** Automatically fetches metadata and downloads the supplementary count matrix from **GEO (GSE152641)**.
* **DESeq2 Analysis:** Implements the standard `DESeq2` workflow for robust DGE testing.
* **Critical Data Alignment:** Includes necessary steps to align sample names between the downloaded count matrix and the GEO metadata, a common challenge in public data analysis.
* **QC Transformation:** Applies **Variance Stabilizing Transformation (VST)** for accurate distance-based plotting.
* **Integrated Visualization:** Generates three essential plots for interpreting the results: **PCA**, **Heatmap**, and **Volcano Plot**.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE152641 | Transcriptomic Host Response in COVID-19 vs. Other Viral Infections. |
| **DGE Tool** | `DESeq2` | Statistical method optimized for count data (RNA-Seq). |
| **Comparison** | COVID-19 vs. Normal (Healthy Control) | Identifies genes regulated by the host response to SARS-CoV-2 infection. |
| **Significance** | $\text{pCutoff} = 0.05$, $\text{FCcutoff} = 1.5$ | Used for highlighting significant findings in the Volcano Plot. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically checks for and installs the necessary Bioconductor and CRAN packages:
* `GEOquery` (For data download)
* `DESeq2` (For DGE analysis)
* `pheatmap` (For Heatmap visualization)
* `EnhancedVolcano` (For Volcano Plot visualization)
* `ggplot2` (For PCA visualization)

### ⚙️ Execution

1.  **Download** the `DGE Analysis of GSE152641 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization.R` file.
2.  **Optional:** Modify the `output_dir` variable (Step 2b) to your preferred saving location.
    ```R
    output_dir <- "D:/DOWNLOADS" # Change this path
    ```
3.  **Execute** the script in your R environment:
    ```R
    source("DGE Analysis of GSE152641 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization.R")
    ```
    *Note: The script will automatically download the necessary supplementary files from GEO into a subdirectory named `GSE152641`.*

---

## 📁 Output Files (3 Plots)

The script automatically saves the following files to the specified `output_dir` (default: `D:/DOWNLOADS`).

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `GSE152641_PCA_plot.png` | QC / Results | **Principal Component Analysis (PCA)** plot demonstrating global clustering and separation of COVID vs. Normal samples. |
| `GSE152641_Heatmap_Top50.png` | QC | **Heatmap of the Top 50 Most Variable Genes** (based on VST-transformed data) to visualize sample grouping quality. |
| `GSE152641_Volcano_plot.png` | Results | **Volcano Plot** generated using `EnhancedVolcano`, showing the $\log_2 \text{Fold Change}$ vs. $P_{\text{value}}$, highlighting significant and highly changed genes. |
