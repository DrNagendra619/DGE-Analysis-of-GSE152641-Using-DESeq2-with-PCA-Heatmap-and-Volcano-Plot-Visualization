# -------------------------------------------------------------------
# Project: DGE Analysis of GSE152641 Using DESeq2 with PCA, 
#          Heatmap, and Volcano Plot Visualization
#
# Dataset: GSE152641
# Title:   Transcriptomic Similarities and Differences in Host Response 
#          between SARS-CoV-2 and Other Viral Infection
# -------------------------------------------------------------------

# Define the GEO Accession ID (This variable will be used throughout)
geo_id <- "GSE152641"

#  Step 1: Package Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "DESeq2", "pheatmap", "EnhancedVolcano", "ggplot2"), update = FALSE)

#  Step 2: Load Libraries
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)


# Step 2b: Define Output Path and Create Directory
# Define your desired output directory.
# R uses forward slashes '/' for paths, even on Windows.
output_dir <- "D:/DOWNLOADS"

# Create the directory if it doesn't already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# Step 3: Fetch Metadata & Download Count Matrix
# This dataset (GSE152641) has separate metadata and count files.

# 3a. Fetch Metadata (pdata) from GEO using the geo_id
gse <- getGEO(geo_id, GSEMatrix = TRUE)
gse <- gse[[1]]
pdata <- pData(gse) # pdata = Phenotype (metadata)

cat(paste("Metadata (pdata) dimensions for", geo_id, ":\n"))
print(dim(pdata)) # Should be 86 samples


# 3b. Download the Supplementary Count Matrix file using the geo_id
# This downloads files into a new directory named after the geo_id (e.g., "GSE152641")
getGEOSuppFiles(geo_id)

# 3c. Read the Count Matrix (exprSet)
# We will build the file path using the geo_id variable for consistency.
count_dir <- geo_id 
count_filename <- paste0(geo_id, "_Inflammatix_COVID19_counts_entrez.csv.gz")

# Combine them into a full path
count_file <- file.path(count_dir, count_filename)

cat("\nReading count file:", count_file, "\n")

# read.csv can read .gz files directly. 
# row.names = 1 means the first column (gene IDs) will be used as rownames.
exprSet <- read.csv(count_file, row.names = 1)

# DESeq2 requires integer counts
exprSet <- round(exprSet)

cat(paste("\nExpression matrix dimensions for", geo_id, ":\n"))
print(dim(exprSet))


# Step 4: Define Sample Condition (CRITICAL)
# We inspect 'pdata' for this dataset (GSE152641)

# 4a. Find the correct column in pdata (it's 'source_name_ch1' for GSE152641)
# You can check with: table(pdata$source_name_ch1)

# 4b. Create the 'condition' column
# We use 'grepl' to find any row that contains the text "healthy control"
pdata$condition <- ifelse(grepl("healthy control", pdata$source_name_ch1, ignore.case = TRUE), 
                          "Normal", 
                          "COVID")

# Convert the condition to a factor (required by DESeq2)
pdata$condition <- as.factor(pdata$condition)

# CRITICAL CHECK: This table MUST show two groups (24 Normal, 62 COVID for GSE152641).
cat("\nSample Group Distribution (Must show 2 groups):\n")
print(table(pdata$condition))


# Step 5: Prepare DESeq2 Dataset (Aligning Data)
# We must align pdata and exprSet for GSE152641.

# 5a. Find the "bridge" column.
# For GSE152641, 'pdata$title' (e.g., "IMX_sample00001") matches
# 'colnames(exprSet)' (e.g., "IMX_sample00001").
# We use this to align them.

if ("title" %in% colnames(pdata)) {
  rownames(pdata) <- pdata$title
} else {
  stop(paste("Could not find 'title' column in pdata for", geo_id, "to match sample names."))
}

# 5b. Find common samples and subset
common_samples <- intersect(rownames(pdata), colnames(exprSet))
cat(paste("\nFound", length(common_samples), "common samples for", geo_id, ".\n"))

# Subset both the matrix and metadata to match AND be in the same order
exprSet_subset <- exprSet[, common_samples]
pdata_subset <- pdata[common_samples, ]

# 5c. CRITICAL CHECK: Ensure order matches
stopifnot(all(colnames(exprSet_subset) == rownames(pdata_subset)))

# 5d. Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = exprSet_subset,
                              colData = pdata_subset,
                              design = ~ condition)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]
cat("\nDESeqDataSet object created successfully for", geo_id, ".\n")


# Step 6: Run DESeq2 Pipeline
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
summary(res)


# Step 7: Variance Stabilizing Transformation (VST)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
cat("\nVST data dimensions:\n")
print(dim(assay(vsd)))


# Step 8: PCA Plot (Dimensional Reduction) & Save
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create the ggplot object
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  ggtitle(paste("PCA:", geo_id, "- COVID vs Normal"))

# Define the filename and save the plot
pca_filename <- file.path(output_dir, paste0(geo_id, "_PCA_plot.png"))
ggsave(filename = pca_filename, plot = pca_plot)
cat("\nPCA plot saved to:", pca_filename, "\n")

# Print the plot to the RStudio plot pane
print(pca_plot)


# Step 9: Heatmap of Top 50 Variable Genes & Save
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

# Define the filename for the heatmap
heatmap_filename <- file.path(output_dir, paste0(geo_id, "_Heatmap_Top50.png"))

# pheatmap can save directly with the 'filename' argument
# It will also still print to the plot pane
pheatmap(assay(vsd)[topVarGenes, ],
         cluster_rows = TRUE,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "condition", drop = FALSE]),
         main = paste(geo_id, ": Top 50 Variable Genes"),
         filename = heatmap_filename)

cat("Heatmap saved to:", heatmap_filename, "\n")


# Step 10: Volcano Plot & Save
# Create the EnhancedVolcano plot (it's a ggplot object)
volcano_plot <- EnhancedVolcano(res,
                                lab = rownames(res),
                                x = 'log2FoldChange',
                                y = 'pvalue',
                                pCutoff = 0.05,
                                FCcutoff = 1.5,
                                title = paste(geo_id, ': Volcano Plot - COVID vs Normal'),
                                subtitle = 'Differential Gene Expression',
                                legendLabels = c('NS', 'Log2FC', 'p-value', 'p-value & Log2FC'))

# Define the filename and save the plot
volcano_filename <- file.path(output_dir, paste0(geo_id, "_Volcano_plot.png"))
ggsave(filename = volcano_filename, plot = volcano_plot)
cat("Volcano plot saved to:", volcano_filename, "\n")

# Print the plot to the RStudio plot pane
print(volcano_plot)
