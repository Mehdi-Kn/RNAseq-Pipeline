# RNA-Seq Pipeline Scripts

Below are the R scripts for each step of the RNA-Seq pipeline. Save each script in the `scripts/` directory of your repository with the respective filenames.
---
### `scripts/01_qc.R`

```r
# Quality Control Script
# Load necessary library
if (!requireNamespace("fastqcr", quietly = TRUE)) {
    install.packages("fastqcr")
}
library(fastqcr)

# Define input and output directories
input_dir <- "data/raw/"
output_dir <- "results/qc_reports/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Run QC
qc_report <- fastqc(dir = input_dir, output = output_dir)

# Save QC report
saveRDS(qc_report, file = file.path(output_dir, "qc_report.rds"))

# Tip: Open the QC reports in your browser to quickly spot any issues!
