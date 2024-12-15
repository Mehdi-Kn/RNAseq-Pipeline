# Differential Expression Analysis Script

# Load DESeq2 package
if (!requireNamespace("DESeq2", quietly = TRUE)) {
    BiocManager::install("DESeq2")
}
library(DESeq2)

# Load normalized counts
normalized_counts <- read.table("results/counts/normalized_counts.txt", header = TRUE, row.names = 1)

# Define sample information (ensure it matches normalization step)
sample_info <- data.frame(
    row.names = colnames(normalized_counts),
    condition = c("control", "treated", "control", "treated")  # Example conditions
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = normalized_counts,
                              colData = sample_info,
                              design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Save results to CSV
write.csv(as.data.frame(res_ordered), file = "results/differential_expression/results.csv")

# Handy tip: Filter out genes with very low counts to reduce noise!
