library(DESeq2)

# ------------------------------------------------------------
# Differential expression analysis with yeast spike-in
# normalization for human TT-seq data
# ------------------------------------------------------------

# 1. Load count matrices
human_raw <- read.delim("human_counts_gene.txt", row.names = 1)
human_counts <- human_raw[, 6:ncol(human_raw)]

yeast_raw <- read.delim("yeast_counts_gene.txt", row.names = 1)
yeast_counts <- yeast_raw[, 6:ncol(yeast_raw)]

# Make sure sample names and sample order are identical
stopifnot(identical(colnames(human_counts), colnames(yeast_counts)))

# 2. Define sample metadata
coldata <- data.frame(
  condition = factor(c("control", "control", "treatment", "treatment")),
  row.names = colnames(human_counts)
)

# 3. Construct DESeq2 objects
dds_human <- DESeqDataSetFromMatrix(
  countData = human_counts,
  colData = coldata,
  design = ~ condition
)

dds_spike <- DESeqDataSetFromMatrix(
  countData = yeast_counts,
  colData = coldata,
  design = ~ condition
)

# Remove genes with zero counts across all samples
dds_human <- dds_human[rowSums(counts(dds_human)) > 0, ]
dds_spike <- dds_spike[rowSums(counts(dds_spike)) > 0, ]

# 4. Estimate size factors from yeast spike-in
dds_spike <- estimateSizeFactors(dds_spike)
sizeFactors(dds_human) <- sizeFactors(dds_spike)

# 5. Run differential expression analysis
dds_human <- DESeq(dds_human)
res <- results(dds_human, contrast = c("condition", "treatment", "control"))
res_df <- as.data.frame(res)

# Optional: save DE results
write.table(
  res_df,
  file = "DESeq2_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# 6. Define significance thresholds
padj_cutoff <- 0.05
lfc_cutoff <- 0.585   # ~1.5-fold change

sig_up <- which(
  !is.na(res_df$padj) &
    res_df$padj < padj_cutoff &
    res_df$log2FoldChange > lfc_cutoff
)

sig_down <- which(
  !is.na(res_df$padj) &
    res_df$padj < padj_cutoff &
    res_df$log2FoldChange < -lfc_cutoff
)

up_count <- length(sig_up)
down_count <- length(sig_down)
total_tested <- sum(!is.na(res_df$padj))

up_pct <- round(up_count / total_tested * 100, 1)
down_pct <- round(down_count / total_tested * 100, 1)

# 7. Draw MA plot directly in R
# Draw MA plot in R
res_df$baseMean_plot <- res_df$baseMean + 1e-1

with(res_df, {
  plot(
    baseMean_plot, log2FoldChange,
    log = "x",
    pch = 20,
    cex = 0.5,
    col = "gray60",
    xlab = "Mean of normalized counts",
    ylab = expression(log[2]~"fold change"),
    main = "TT-seq (swap/26CTD)",
    ylim = c(-4, 4)
  )
  
  points(baseMean_plot[sig_down], log2FoldChange[sig_down],
         pch = 20, cex = 0.5, col = "blue")
  
  points(baseMean_plot[sig_up], log2FoldChange[sig_up],
         pch = 20, cex = 0.5, col = "red")
  
  abline(h = 0, col = "black")
  abline(h = c(-lfc_cutoff, lfc_cutoff), col = "black", lty = 2)
  
  text(5, -3.5,
       labels = paste0("Downregulated = ", down_count, " (", down_pct, "%)"),
       pos = 4)
  
  text(5, 3.5,
       labels = paste0("Upregulated = ", up_count, " (", up_pct, "%)"),
       pos = 4)
})
