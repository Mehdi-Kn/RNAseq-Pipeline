# Functional Analysis and Visualization Script

# Load necessary packages
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {  # Replace with your organism
    BiocManager::install("org.Hs.eg.db")
}
library(clusterProfiler)
library(org.Hs.eg.db)

# Load differential expression results
res <- read.csv("results/differential_expression/results.csv", row.names = 1)

# Filter significant genes (adjust thresholds as needed)
sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]

# Convert gene symbols to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# Perform GO enrichment analysis
go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",  # Biological Process
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

# Save GO enrichment results
write.csv(as.data.frame(go_enrich), file = "results/functional_analysis/go_enrichment.csv")

# Plot GO enrichment
pdf("results/functional_analysis/go_dotplot.pdf")
dotplot(go_enrich, showCategory = 10) + ggtitle("GO Enrichment - Biological Process")
dev.off()

# Handy tip: Customize your plots with different themes to match your presentation style!
