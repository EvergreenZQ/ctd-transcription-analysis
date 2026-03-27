library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# ------------------------------------------------------------
# Analysis of RNAPII traveling ratio using promoter/body CPM
# signals across two conditions: CTD and swap
# ------------------------------------------------------------

# -----------------------------
# 1. Load promoter CPM signal
# -----------------------------
prom <- read.delim("promoter.cpm.tab", header = TRUE, check.names = FALSE)
colnames(prom) <- c("chr", "start", "end", "sample1", "sample2")

bed_prom <- read.delim("promoter.sorted.bed", header = FALSE) %>%
  transmute(chr = V1, start = V2, end = V3, gene_id = V4) %>%
  distinct(chr, start, end, .keep_all = TRUE)

prom_annot <- prom %>%
  inner_join(bed_prom, by = c("chr", "start", "end")) %>%
  dplyr::select(gene_id, dplyr::everything())

# -----------------------------
# 2. Load gene body CPM signal
# -----------------------------
body <- read.delim("body.cpm.tab", header = TRUE, check.names = FALSE)
colnames(body) <- c("chr", "start", "end", "sample1", "sample2")

bed_body <- read.delim("body.sorted.bed", header = FALSE) %>%
  transmute(chr = V1, start = V2, end = V3, gene_id = V4) %>%
  distinct(chr, start, end, .keep_all = TRUE)

body_annot <- body %>%
  inner_join(bed_body, by = c("chr", "start", "end")) %>%
  dplyr::select(gene_id, dplyr::everything())

# Save annotated signal tables
write.table(
  prom_annot,
  file = "prom_counts_annotated.tab",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  body_annot,
  file = "body_counts_annotated.tab",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# 3. Convert to long format
# -----------------------------
# Each row represents one gene in one sample.

prom_long <- prom_annot %>%
  pivot_longer(
    cols = c(sample1, sample2),
    names_to = "sample",
    values_to = "prom_signal"
  )

body_long <- body_annot %>%
  pivot_longer(
    cols = c(sample1, sample2),
    names_to = "sample",
    values_to = "body_signal"
  )

# -----------------------------
# 4. Summarize CPM per gene
# -----------------------------
# If multiple promoter or gene body regions map to the same gene,
# average the CPM values for each sample.

prom_sum <- prom_long %>%
  group_by(gene_id, sample) %>%
  summarise(prom_signal = mean(prom_signal, na.rm = TRUE), .groups = "drop")

body_sum <- body_long %>%
  group_by(gene_id, sample) %>%
  summarise(body_signal = mean(body_signal, na.rm = TRUE), .groups = "drop")

signal_df <- inner_join(prom_sum, body_sum, by = c("gene_id", "sample"))

# -----------------------------
# 5. Calculate traveling ratio
# -----------------------------
# Traveling ratio is defined here as:
# promoter signal / gene body signal
#
# Only genes with positive promoter and body CPM are retained.

tr_df <- signal_df %>%
  filter(prom_signal > 0, body_signal > 0) %>%
  mutate(
    PI = prom_signal / body_signal,
    log2PI = log2(PI),
    sample = recode(sample, "sample1" = "CTD", "sample2" = "swap")
  )

summary(tr_df$log2PI)

# Kolmogorov-Smirnov test comparing overall distributions
ks_ecdf <- ks.test(
  tr_df$log2PI[tr_df$sample == "CTD"],
  tr_df$log2PI[tr_df$sample == "swap"]
)
print(ks_ecdf)

# -----------------------------
# 6. Keep genes quantified in both conditions
# -----------------------------
wide_df <- tr_df %>%
  dplyr::select(gene_id, sample, log2PI) %>%
  tidyr::pivot_wider(names_from = sample, values_from = log2PI) %>%
  tidyr::drop_na(CTD, swap)

clean_df <- wide_df %>%
  pivot_longer(
    cols = c(CTD, swap),
    names_to = "sample",
    values_to = "log2PI"
  ) %>%
  mutate(sample = factor(sample, levels = c("CTD", "swap")))

# -----------------------------
# 7. Plot ECDF
# -----------------------------
p_ecdf <- ggplot(clean_df, aes(x = log2PI, color = sample)) +
  stat_ecdf(linewidth = 1) +
  coord_cartesian(xlim = c(-8, 10)) +
  scale_x_continuous(breaks = seq(-4, 10, by = 4)) +
  scale_color_manual(
    name = "Group",
    values = c("CTD" = "#459af8", "swap" = "orange")
  ) +
  labs(
    x = "Log2 Traveling Ratio",
    y = "Fraction of Genes",
    title = "RNAPII Traveling Ratio"
  ) +
  theme_classic(base_size = 14)

print(p_ecdf)

ggsave("RNAPII_traveling_ratio_ECDF.pdf", p_ecdf, width = 6, height = 5)
ggsave("RNAPII_traveling_ratio_ECDF.png", p_ecdf, width = 6, height = 5, dpi = 300)

# -----------------------------
# 8. Compare group medians
# -----------------------------
medians <- clean_df %>%
  group_by(sample) %>%
  summarise(median_log2PI = median(log2PI), .groups = "drop")

print(medians)

wilcox_res <- wilcox.test(log2PI ~ sample, data = clean_df)
print(wilcox_res)

# -----------------------------
# 9. Plot boxplot
# -----------------------------
p_box <- ggplot(clean_df, aes(x = sample, y = log2PI, fill = sample)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.7) +
  scale_fill_manual(values = c("CTD" = "#1f77b4", "swap" = "#ff7f0e")) +
  coord_cartesian(ylim = c(-4.5, 4)) +
  labs(
    x = NULL,
    y = "Log2 Traveling Ratio",
    title = "Traveling Ratio"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

print(p_box)

ggsave("RNAPII_traveling_ratio_boxplot.pdf", p_box, width = 4.5, height = 5)
ggsave("RNAPII_traveling_ratio_boxplot.png", p_box, width = 4.5, height = 5, dpi = 300)

# -----------------------------
# 10. Classify genes by change magnitude
# -----------------------------
# delta = swap - CTD
# Genes are classified into strong / medium / weak groups based on
# the absolute delta quantiles.

delta_df <- wide_df %>%
  mutate(
    delta = swap - CTD,
    abs_delta = abs(delta)
  )

summary(delta_df$delta)

q80 <- quantile(delta_df$abs_delta, 0.80, na.rm = TRUE)
q50 <- quantile(delta_df$abs_delta, 0.50, na.rm = TRUE)

delta_df <- delta_df %>%
  mutate(
    strength = case_when(
      abs_delta >= q80 ~ "strong",
      abs_delta >= q50 ~ "medium",
      TRUE ~ "weak"
    ),
    direction = case_when(
      delta > 0 ~ "up_in_swap",
      delta < 0 ~ "down_in_swap",
      TRUE ~ "no_change"
    )
  )

print(table(delta_df$strength))
print(table(delta_df$strength, delta_df$direction))

# Export gene lists
write.table(
  delta_df %>% filter(strength == "strong") %>% pull(gene_id),
  file = "TR_strong_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  delta_df %>% filter(strength == "strong", direction == "up_in_swap") %>% pull(gene_id),
  file = "TR_strong_up_in_swap_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  delta_df %>% filter(strength == "strong", direction == "down_in_swap") %>% pull(gene_id),
  file = "TR_strong_down_in_swap_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# -----------------------------
# 11. GO enrichment analysis
# -----------------------------
# Remove Ensembl version suffixes before ID conversion.

delta_df <- delta_df %>%
  mutate(gene_id_clean = sub("\\.\\d+$", "", gene_id))

bg_ensg <- delta_df$gene_id_clean
strong_ensg <- delta_df %>% filter(strength == "strong") %>% pull(gene_id_clean)
up_ensg <- delta_df %>% filter(strength == "strong", direction == "up_in_swap") %>% pull(gene_id_clean)
down_ensg <- delta_df %>% filter(strength == "strong", direction == "down_in_swap") %>% pull(gene_id_clean)

bg_map <- bitr(bg_ensg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
strong_map <- bitr(strong_ensg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
up_map <- bitr(up_ensg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_map <- bitr(down_ensg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

bg_entrez <- unique(bg_map$ENTREZID)
strong_entrez <- unique(strong_map$ENTREZID)
up_entrez <- unique(up_map$ENTREZID)
down_entrez <- unique(down_map$ENTREZID)

cat(
  "Background mapped:", length(bg_entrez), "\n",
  "Strong mapped:", length(strong_entrez), "\n",
  "Up mapped:", length(up_entrez), "\n",
  "Down mapped:", length(down_entrez), "\n"
)

ego_down_bp <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1,
  readable = TRUE
)

print(head(as.data.frame(ego_down_bp)))

p_go <- dotplot(ego_down_bp, showCategory = 10) +
  ggtitle("GO BP enrichment: strong down in swap") +
  theme(plot.title = element_text(hjust = 0.5))

print(p_go)

ggsave("GO_BP_strong_down_in_swap.pdf", p_go, width = 8, height = 6)
ggsave("GO_BP_strong_down_in_swap.png", p_go, width = 8, height = 6, dpi = 300)