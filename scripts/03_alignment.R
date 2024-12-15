# Alignment/Mapping Script

# Load Rsubread package
if (!requireNamespace("Rsubread", quietly = TRUE)) {
    BiocManager::install("Rsubread")
}
library(Rsubread)

# Define reference genome and build index if not already done
index_basename <- "reference/genome_index"
genome_fasta <- "reference/genome.fa"

if (!file.exists(paste0(index_basename, ".1.bt2"))) {
    buildindex(basename = index_basename, reference = genome_fasta)
    # Note: Indexing only needs to be done once unless your reference changes!
}

# Function to align reads
align_reads <- function(input_fastq, output_bam, index) {
    align(index = index,
          readfile1 = input_fastq,
          output_file = output_bam,
          nthreads = 4)
}

# Define directories
input_dir <- "data/processed/"
output_dir <- "results/alignment/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Align all processed FASTQ files
processed_files <- list.files(input_dir, pattern = "\\.fastq$", full.names = TRUE)

for (file in processed_files) {
    bam_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file)), ".bam"))
    align_reads(file, bam_file, index_basename)
}
 # Handy tip: Monitor CPU usage to optimize thread allocation!
