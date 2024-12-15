# Normalization Script

# Load DESeq2 package
if (!requireNamespace("DESeq2", quietly = TRUE)) {
    BiocManager::install("DESeq2")
}
library(DESeq2)

# Function to create count matrix
create_count_matrix <- function(counts_dir) {
    count_files <- list.files(counts_dir, pattern = "_counts\\.txt$", full.names = TRUE)
    counts_list <- lapply(count_files, function(file) {
        counts <- read.table(file, header = TRUE, row.names = 1)
        return(counts[,1])
    })
    count_matrix <- do.call(cbind, counts_list)
    colnames(count_matrix) <- tools::file_path_sans_ext(basename(count_files))
    return(count_matrix)
}

# Create count matrix
counts_dir <- "results/counts/"
count_matrix <- create_count_matrix(counts_dir)

# Define sample information (modify according to your experimental design)
sample_info <- data.frame(
    row.names = colnames(count_matrix),
    condition = c("control", "treated", "control", "treated")  # Example conditions
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)

# Perform normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Save normalized counts
write.table(normalized_counts, "results/counts/normalized_counts.txt", sep = "\t", quote = FALSE, col.names = NA)

# Quick tip: Always visualize normalized counts to ensure proper scaling!
