library(tidyverse)

# ------------------------------------------------------------
# Readthrough index analysis for two conditions: yCTD and Swap
#
# This script:
# 1. Loads TES+150 and TES-150 signal tables from plus/minus strands
# 2. Annotates genomic intervals with gene IDs
# 3. Merges strand-specific signal tables
# 4. Calculates gene-level readthrough index
# 5. Compares distributions between conditions
# 6. Generates ECDF and boxplot visualizations
# ------------------------------------------------------------

# -----------------------------
# 1. Load and annotate TES+150 signal
# -----------------------------
prom_plus <- read.delim("up150.plus.tab", header = TRUE, check.names = FALSE)
colnames(prom_plus) <- c("chr", "start", "end", "sample1", "sample2")

bed_prom_plus <- read.delim("up150.plus.clean.bed6", header = FALSE) %>%
  dplyr::transmute(chr = V1, start = V2, end = V3, gene_id = V4)

bed_prom_plus_unique <- bed_prom_plus %>%
  dplyr::distinct(chr, start, end, .keep_all = TRUE)

prom_plus_annot <- prom_plus %>%
  dplyr::inner_join(bed_prom_plus_unique, by = c("chr", "start", "end")) %>%
  dplyr::select(gene_id, dplyr::everything())

prom_minus <- read.delim("up150.minus.tab", header = TRUE, check.names = FALSE)
colnames(prom_minus) <- c("chr", "start", "end", "sample1", "sample2")

bed_prom_minus <- read.delim("up150.minus.clean.bed6", header = FALSE) %>%
  dplyr::transmute(chr = V1, start = V2, end = V3, gene_id = V4)

bed_prom_minus_unique <- bed_prom_minus %>%
  dplyr::distinct(chr, start, end, .keep_all = TRUE)

prom_minus_annot <- prom_minus %>%
  dplyr::inner_join(bed_prom_minus_unique, by = c("chr", "start", "end")) %>%
  dplyr::select(gene_id, dplyr::everything())

# -----------------------------
# 2. Load and annotate TES-150 signal
# -----------------------------
body_plus <- read.delim("dn150.plus.tab", header = TRUE, check.names = FALSE)
colnames(body_plus) <- c("chr", "start", "end", "sample1", "sample2")

bed_body_plus <- read.delim("dn150.plus.clean.bed6", header = FALSE) %>%
  dplyr::transmute(chr = V1, start = V2, end = V3, gene_id = V4)

bed_body_plus_unique <- bed_body_plus %>%
  dplyr::distinct(chr, start, end, .keep_all = TRUE)

body_plus_annot <- body_plus %>%
  dplyr::inner_join(bed_body_plus_unique, by = c("chr", "start", "end")) %>%
  dplyr::select(gene_id, dplyr::everything())

body_minus <- read.delim("dn150.minus.tab", header = TRUE, check.names = FALSE)
colnames(body_minus) <- c("chr", "start", "end", "sample1", "sample2")

bed_body_minus <- read.delim("dn150.minus.clean.bed6", header = FALSE) %>%
  dplyr::transmute(chr = V1, start = V2, end = V3, gene_id = V4)

bed_body_minus_unique <- bed_body_minus %>%
  dplyr::distinct(chr, start, end, .keep_all = TRUE)

body_minus_annot <- body_minus %>%
  dplyr::inner_join(bed_body_minus_unique, by = c("chr", "start", "end")) %>%
  dplyr::select(gene_id, dplyr::everything())

# -----------------------------
# 3. Merge strand-specific annotated tables
# -----------------------------
prom_all <- dplyr::bind_rows(
  prom_plus_annot %>% dplyr::mutate(source = "plus150_forward"),
  prom_minus_annot %>% dplyr::mutate(source = "plus150_reverse")
) %>%
  dplyr::arrange(chr, start)

body_all <- dplyr::bind_rows(
  body_plus_annot %>% dplyr::mutate(source = "minus150_forward"),
  body_minus_annot %>% dplyr::mutate(source = "minus150_reverse")
) %>%
  dplyr::arrange(chr, start)

prom_out <- prom_all %>%
  dplyr::select(gene_id, chr, start, end, sample1, sample2, source)

body_out <- body_all %>%
  dplyr::select(gene_id, chr, start, end, sample1, sample2, source)

# Save annotated merged tables
write.table(
  prom_out,
  file = "average_TES_plus150_counts_annotated_merged.tab",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

write.table(
  body_out,
  file = "average_TES_minus150_counts_annotated_merged.tab",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# Save forward-strand annotated tables
write.table(
  prom_plus_annot,
  file = "average_TES_plus150_counts_annotated_fstand.tab",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  body_plus_annot,
  file = "average_TES_minus150_counts_annotated_fstand.tab",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 4. Convert merged tables to long format
# -----------------------------
prom_long <- prom_out %>%
  tidyr::pivot_longer(
    cols = c(sample1, sample2),
    names_to = "sample",
    values_to = "prom_signal"
  )

body_long <- body_out %>%
  tidyr::pivot_longer(
    cols = c(sample1, sample2),
    names_to = "sample",
    values_to = "body_signal"
  )

# -----------------------------
# 5. Summarize signal per gene and sample
# -----------------------------
# If multiple intervals map to the same gene, average the signal.

prom_sum <- prom_long %>%
  dplyr::group_by(gene_id, sample) %>%
  dplyr::summarise(
    prom_signal = mean(prom_signal, na.rm = TRUE),
    .groups = "drop"
  )

body_sum <- body_long %>%
  dplyr::group_by(gene_id, sample) %>%
  dplyr::summarise(
    body_signal = mean(body_signal, na.rm = TRUE),
    .groups = "drop"
  )

signal_df <- dplyr::inner_join(prom_sum, body_sum, by = c("gene_id", "sample"))

# -----------------------------
# 6. Calculate readthrough index
# -----------------------------
# Readthrough index is defined here as:
# TES-150 signal / TES+150 signal

rt_df <- signal_df %>%
  dplyr::filter(prom_signal > 0, body_signal > 0) %>%
  dplyr::mutate(
    RI = body_signal / prom_signal,
    log2RI = log2(RI)
  )

summary(rt_df$log2RI)

# Rename samples for plotting
rt_df <- rt_df %>%
  dplyr::mutate(
    sample = dplyr::recode(
      sample,
      "sample1" = "yCTD",
      "sample2" = "Swap"
    )
  )

# -----------------------------
# 7. Retain genes quantified in both conditions
# -----------------------------
wide_df <- rt_df %>%
  dplyr::select(gene_id, sample, log2RI) %>%
  tidyr::pivot_wider(
    names_from = sample,
    values_from = log2RI
  ) %>%
  tidyr::drop_na(yCTD, Swap)

clean_df <- wide_df %>%
  tidyr::pivot_longer(
    cols = c(yCTD, Swap),
    names_to = "sample",
    values_to = "log2RI"
  )

# -----------------------------
# 8. Plot ECDF
# -----------------------------
p_ecdf <- ggplot(clean_df, aes(x = log2RI, color = sample)) +
  stat_ecdf(linewidth = 1) +
  coord_cartesian(xlim = c(-4, 10)) +
  scale_x_continuous(breaks = seq(-4, 10, by = 4)) +
  scale_color_manual(
    name = "Group",
    values = c(
      "yCTD" = "#459af8",
      "Swap" = "orange"
    )
  ) +
  labs(
    x = "Log2 Readthrough Index",
    y = "Fraction of Genes",
    title = "RNAPII Readthrough Index"
  ) +
  theme_classic(base_size = 14)

print(p_ecdf)

# -----------------------------
# 9. Statistical comparison
# -----------------------------
medians <- clean_df %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(median_log2RI = median(log2RI), .groups = "drop")

print(medians)

# Use paired Wilcoxon test if the two conditions are paired by gene
p_value <- wilcox.test(
  wide_df$Swap,
  wide_df$yCTD,
  paired = TRUE
)$p.value

cat("Paired Wilcoxon p-value:", signif(p_value, 3), "\n")

# -----------------------------
# 10. Draw boxplot
# -----------------------------
plot_df <- clean_df %>%
  dplyr::mutate(sample = factor(sample, levels = c("yCTD", "Swap")))

cols <- c("#1f77b4", "#ff7f0e")
ylim_range <- c(-4.5, 4)

par(mar = c(6, 6, 4, 1))

boxplot(
  log2RI ~ sample,
  data = plot_df,
  col = cols,
  border = "black",
  outline = FALSE,
  lwd = 2,
  las = 1,
  ylab = "Log2 Readthrough Index",
  main = "Readthrough Index",
  ylim = ylim_range,
  boxwex = 0.55,
  medlwd = 2
)
