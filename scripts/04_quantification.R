# Quantification Script

# Load Rsubread package
if (!requireNamespace("Rsubread", quietly = TRUE)) {
    BiocManager::install("Rsubread")
}
library(Rsubread)

# Define annotation file
gtf_file <- "reference/annotation.gtf"

# Function to count reads per gene
count_features <- function(bam_file, gtf, output_file) {
    fc <- featureCounts(files = bam_file,
                       annot.ext = gtf,
                       isGTFAnnotationFile = TRUE,
                       GTF.featureType = "exon",
                       GTF.attrType = "gene_id",
                       nthreads = 4)
    write.table(fc$counts, file = output_file, sep = "\t", quote = FALSE, col.names = NA)
}

# Define directories
bam_dir <- "results/alignment/"
counts_dir <- "results/counts/"
dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)

# Quantify all BAM files
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

for (bam in bam_files) {
    output_file <- file.path(counts_dir, paste0(tools::file_path_sans_ext(basename(bam)), "_counts.txt"))
    count_features(bam, gtf_file, output_file)
    # Tip: Keep your count files organized by sample for easy access!
}
