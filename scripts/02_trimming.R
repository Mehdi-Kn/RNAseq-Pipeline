# Read Trimming and Filtering Script

# Load ShortRead package
if (!requireNamespace("ShortRead", quietly = TRUE)) {
    BiocManager::install("ShortRead")
}
library(ShortRead)

# Function to trim reads based on quality
trim_reads <- function(input_fastq, output_fastq, quality_threshold = 20) {
    fq <- readFastq(input_fastq)
    fq_trimmed <- trimTailw(fq, threshold = quality_threshold)
    writeFastq(fq_trimmed, output_fastq)
}

# Define directories
input_dir <- "data/raw/"
output_dir <- "data/processed/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# List all FASTQ files
fastq_files <- list.files(input_dir, pattern = "\\.fastq$", full.names = TRUE)

# Apply trimming to each file
for (file in fastq_files) {
    output_file <- file.path(output_dir, basename(file))
    trim_reads(file, output_file)
}
# Pro tip: Rename files if needed to keep track of different samples!
