# Consensus clustering and external validation

PROJECT_DIR <- "."

gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GSEABase)
  library(ConsensusClusterPlus)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(vegan)
  library(reshape2)
  library(grid)
  library(gridExtra)
})

workDir  <- file.path(PROJECT_DIR, "results/bulk/discovery")
exprFile <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_expr_HighLow.csv")
geneFile <- file.path(PROJECT_DIR, "results/bulk/discovery/03_HVG_blue_module_gene_list.txt")

dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
setwd(workDir)
dir.create("01_Result", showWarnings = FALSE)
dir.create("02_Figures", showWarnings = FALSE)
dir.create("03_Tables", showWarnings = FALSE)

col_c1 <- "#b43665"
col_c2 <- "#6fa6cf"
cluster_colors_2 <- c("1" = col_c1, "2" = col_c2)
cluster_colors_named <- c("C1" = col_c1, "C2" = col_c2)

hm_colors <- colorRampPalette(c("white", "#f4d6e2", col_c1))(100)

maxK       <- 9
reps_num   <- 100
pItem_num  <- 0.8
pFeature_num <- 1
seed_num   <- 2025

expr_raw <- read.csv(exprFile,
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

expr_mat <- as.matrix(expr_raw)
mode(expr_mat) <- "numeric"

gene_list <- read.table(geneFile,
                        header = FALSE,
                        stringsAsFactors = FALSE,
                        sep = "\t",
                        quote = "",
                        fill = TRUE)

rloop_genes <- unique(na.omit(as.character(gene_list[, 1])))
rloop_genes <- rloop_genes[rloop_genes != ""]

matched_genes <- intersect(rloop_genes, rownames(expr_mat))
missing_genes <- setdiff(rloop_genes, rownames(expr_mat))

write.table(data.frame(Missing_Genes = missing_genes),
            file = "03_Tables/Missing_Rloop_genes.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

expr_rloop <- expr_mat[matched_genes, , drop = FALSE]

gene_var <- apply(expr_rloop, 1, var, na.rm = TRUE)
zero_var_genes <- names(gene_var[gene_var == 0])

if (length(zero_var_genes) > 0) {
  expr_rloop <- expr_rloop[!(rownames(expr_rloop) %in% zero_var_genes), , drop = FALSE]
}

write.csv(expr_rloop, file = "03_Tables/Rloop_expr_used_for_clustering.csv", quote = FALSE)

results <- ConsensusClusterPlus(
  d = expr_rloop,
  maxK = maxK,
  reps = reps_num,
  pItem = pItem_num,
  pFeature = pFeature_num,
  title = file.path(workDir, "01_Result"),
  clusterAlg = "pam",
  distance = "euclidean",
  seed = seed_num,
  plot = "pdf"
)

Kvec <- 2:maxK
PAC <- sapply(Kvec, function(k) {
  M <- results[[k]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  Fn(0.9) - Fn(0.1)
})

pac_df <- data.frame(
  K = Kvec,
  PAC = PAC
)

write.csv(pac_df, file = "03_Tables/PAC_values.csv", row.names = FALSE)

optK_pac <- Kvec[which.min(PAC)]

optK_val <- 2

cluster_num <- results[[optK_val]]$consensusClass
cluster_df <- data.frame(
  sample = names(cluster_num),
  Cluster_num = cluster_num,
  Cluster = paste0("C", cluster_num),
  stringsAsFactors = FALSE
)

write.csv(cluster_df,
          file = sprintf("03_Tables/ConsensusCluster_K%d_assignment.csv", optK_val),
          row.names = FALSE)

print(table(cluster_df$Cluster))

p_pac <- ggplot(pac_df, aes(x = K, y = PAC)) +
  geom_line(color = "#444444", linewidth = 1.1) +
  geom_point(size = 3.2, color = col_c1) +
  geom_vline(xintercept = optK_val, linetype = "solid", color = col_c1, linewidth = 0.9) +
  geom_vline(xintercept = optK_pac, linetype = "dashed", color = col_c2, linewidth = 0.9) +
  annotate("text", x = optK_val, y = max(pac_df$PAC) * 0.96,
           label = paste0("Selected K = ", optK_val),
           color = col_c1, hjust = -0.1, vjust = 1, size = 4.2, fontface = "bold") +
  annotate("text", x = optK_pac, y = max(pac_df$PAC) * 0.84,
           label = paste0("Min PAC K = ", optK_pac),
           color = col_c2, hjust = -0.1, vjust = 1, size = 4.0, fontface = "bold") +
  scale_x_continuous(breaks = Kvec) +
  labs(title = "PAC Curve for Consensus Clustering",
       x = "Number of clusters (K)",
       y = "Proportion of ambiguous clustering") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black")
  )

ggsave("02_Figures/Figure_S1_PAC_curve.pdf", p_pac, width = 7.5, height = 5.8)

cdf_data <- data.frame()
for (k in 2:maxK) {
  m <- results[[k]]$consensusMatrix
  vals <- m[lower.tri(m, diag = FALSE)]
  tmp <- data.frame(
    K = factor(rep(k, length(vals))),
    Consensus = vals
  )
  cdf_data <- rbind(cdf_data, tmp)
}

p_cdf <- ggplot(cdf_data, aes(x = Consensus, color = K)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(title = "Consensus Cumulative Distribution Function",
       x = "Consensus index",
       y = "Cumulative distribution",
       color = "K") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black")
  ) +
  scale_color_brewer(palette = "Set1")

ggsave("02_Figures/Figure_S2_CDF_curve.pdf", p_cdf, width = 8, height = 6)

consensus_mat <- results[[optK_val]]$consensusMatrix
colnames(consensus_mat) <- names(cluster_num)
rownames(consensus_mat) <- names(cluster_num)

ann_col <- data.frame(Cluster = factor(paste0("C", cluster_num), levels = c("C1", "C2")))
rownames(ann_col) <- names(cluster_num)

ann_colors <- list(Cluster = cluster_colors_named)

pdf("02_Figures/Figure_3A_Consensus_heatmap_K2.pdf", width = 9.2, height = 8.5)
pheatmap(
  consensus_mat,
  color = colorRampPalette(c("white", "#d9e7f3", col_c2, col_c1))(100),
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  main = "Consensus Matrix Heatmap (K = 2)",
  fontsize = 12
)
dev.off()

cluster_tab <- do.call(cbind, lapply(2:maxK, function(k) {
  results[[k]]$consensusClass
}))
colnames(cluster_tab) <- paste0("K", 2:maxK)

final_cluster <- cluster_tab[, paste0("K", optK_val)]
cluster_tab <- cluster_tab[order(final_cluster), , drop = FALSE]

melt_data <- reshape2::melt(cluster_tab)
colnames(melt_data) <- c("Sample", "K", "Cluster")

nclust_max <- max(cluster_tab)
pal_dynamic <- colorRampPalette(c(col_c2, "#d9d9d9", col_c1))(max(6, nclust_max))

p_track <- ggplot(melt_data, aes(x = K, y = Sample, fill = factor(Cluster))) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = pal_dynamic, name = "Cluster") +
  labs(title = paste0("Cluster Membership Tracking Across K (Selected K = ", optK_val, ")"),
       x = "Number of clusters",
       y = "Samples") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black")
  )

ggsave("02_Figures/Figure_S3_Tracking_plot.pdf", p_track, width = 8.8, height = 7.6)

pca_res <- prcomp(t(expr_rloop), center = TRUE, scale. = FALSE)
pca_data <- as.data.frame(pca_res$x[, 1:2])
pca_data$sample <- rownames(pca_data)

pca_data <- pca_data %>%
  left_join(cluster_df[, c("sample", "Cluster")], by = "sample")

percent_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
cum_var <- sum(percent_var)

adonis_result <- adonis2(dist(pca_data[, c("PC1", "PC2")]) ~ Cluster, data = pca_data, permutations = 999)
p_perm <- adonis_result$`Pr(>F)`[1]

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3.2, alpha = 0.88) +
  stat_ellipse(level = 0.95, linewidth = 0.9, aes(fill = Cluster), geom = "polygon", alpha = 0.12, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors_named) +
  scale_fill_manual(values = cluster_colors_named) +
  labs(title = "Principal Component Analysis of R-loop Expression Pattern",
       subtitle = paste0("PERMANOVA p = ", format.pval(p_perm, digits = 3),
                         " | Total variance explained = ", cum_var, "%"),
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", color = "black")
  )

ggsave("02_Figures/Figure_3B_PCA_K2.pdf", p_pca, width = 7.5, height = 6.2)

gene_var2 <- apply(expr_rloop, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_var2, decreasing = TRUE))[1:min(100, length(gene_var2))]

expr_top <- expr_rloop[top_genes, cluster_df$sample, drop = FALSE]
expr_top_scaled <- t(scale(t(expr_top)))
expr_top_scaled[is.na(expr_top_scaled)] <- 0

ann_col2 <- data.frame(Cluster = factor(cluster_df$Cluster, levels = c("C1", "C2")))
rownames(ann_col2) <- cluster_df$sample

pdf("02_Figures/Figure_3C_TopVariableGenes_heatmap.pdf", width = 8.8, height = 9.5)
pheatmap(
  expr_top_scaled,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  annotation_col = ann_col2,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  border_color = NA,
  fontsize = 11,
  main = "Top Variable Genes Across Clusters"
)
dev.off()

report_lines <- c(
  "Gene-set consensus clustering summary",
  "===================================",
  paste0("Expression file: ", exprFile),
  paste0("Gene list file: ", geneFile),
  paste0("Original expression matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples"),
  paste0("Total genes in file: ", length(rloop_genes)),
  paste0("Matched genes: ", length(matched_genes)),
  paste0("Missing genes: ", length(missing_genes)),
  paste0("Zero-variance genes removed: ", length(zero_var_genes)),
  paste0("Final clustering matrix: ", nrow(expr_rloop), " genes x ", ncol(expr_rloop), " samples"),
  paste0("Consensus clustering parameters: maxK=", maxK,
         ", reps=", reps_num,
         ", pItem=", pItem_num,
         ", pFeature=", pFeature_num,
         ", algorithm=pam, distance=euclidean"),
  paste0("Selected K: ", optK_val),
  paste0("PAC-optimal K: ", optK_pac),
  "",
  "Cluster distribution:",
  capture.output(print(table(cluster_df$Cluster)))
)

writeLines(report_lines, con = "03_Tables/Clustering_summary_report.txt")

gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(limma)
})

workDir <- file.path(PROJECT_DIR, "results/bulk/discovery")
setwd(workDir)

assign_file <- file.path(PROJECT_DIR, "results/bulk/discovery/03_Tables/ConsensusCluster_K2_assignment.csv")
expr_file   <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_expr_HighLow.csv")
clin_file   <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_HRp_HERn_clinical.csv")

out_dir <- file.path(workDir, "04_C1C2_DEG_Analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assign_df <- read.csv(assign_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "Cluster") %in% colnames(assign_df))) {
  stop("Assignment file must contain columns: sample, Cluster")
}

assign_df$sample  <- trimws(as.character(assign_df$sample))
assign_df$Cluster <- trimws(as.character(assign_df$Cluster))

assign_df$sample_base <- sub("_[^_]+$", "", assign_df$sample)

assign_df$new_sample <- paste0(assign_df$sample_base, "_", assign_df$Cluster)

assign_out <- assign_df[, c("new_sample", "Cluster_num", "Cluster", "sample", "sample_base")]
colnames(assign_out) <- c("sample", "Cluster_num", "Cluster", "old_sample", "sample_base")

write.csv(
  assign_out,
  file = file.path(out_dir, "ConsensusCluster_K2_GSE96058.csv"),
  row.names = FALSE
)

expr_df <- read.csv(expr_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

gene_col <- colnames(expr_df)[1]
colnames(expr_df)[1] <- "GeneSymbol"

expr_mat <- expr_df

old_expr_samples <- colnames(expr_mat)[-1]

map_df <- assign_df[, c("sample", "new_sample")]
sample_map <- setNames(map_df$new_sample, map_df$sample)

matched_new_names <- sample_map[old_expr_samples]

keep_cols <- !is.na(matched_new_names)

expr_mat2 <- expr_mat[, c(TRUE, keep_cols), drop = FALSE]
new_names <- c("GeneSymbol", unname(matched_new_names[keep_cols]))
colnames(expr_mat2) <- new_names

write.csv(
  expr_mat2,
  file = file.path(out_dir, "GSE96058_expr_C1C2.csv"),
  row.names = FALSE
)

clinical <- read.csv(clin_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "Rloop_ssGSEA_Score") %in% colnames(clinical))) {
  stop("Clinical file must contain columns: sample, Rloop_ssGSEA_Score")
}

clinical$sample <- trimws(as.character(clinical$sample))

clinical2 <- merge(
  clinical,
  assign_df[, c("sample_base", "Cluster")],
  by.x = "sample",
  by.y = "sample_base",
  all.x = FALSE,
  all.y = FALSE
)

clinical2$sample_C1C2 <- paste0(clinical2$sample, "_", clinical2$Cluster)

clinical_out <- clinical2
colnames(clinical_out)[colnames(clinical_out) == "sample_C1C2"] <- "sample_new"

write.csv(
  clinical_out,
  file = file.path(out_dir, "GSE96058_clinical_C1C2.csv"),
  row.names = FALSE
)

plot_df <- clinical2 %>%
  mutate(
    Cluster = factor(Cluster, levels = c("C1", "C2")),
    Rloop_ssGSEA_Score = as.numeric(Rloop_ssGSEA_Score)
  ) %>%
  filter(!is.na(Cluster), !is.na(Rloop_ssGSEA_Score))

if (length(unique(plot_df$Cluster)) < 2) {
  stop("Need both C1 and C2 for violin plot.")
}

p_violin <- ggplot(plot_df, aes(x = Cluster, y = Rloop_ssGSEA_Score, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", color = "black",
               outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.4, color = "grey25") +
  scale_fill_manual(values = c("C1" = "#6fa6cf", "C2" = "#b43665")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 5
  ) +
  labs(
    x = "Consensus Cluster",
    y = "R-loop ssGSEA Score",
    title = "R-loop score difference between C1 and C2"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = file.path(out_dir, "GSE96058_Rloop_C1C2_violin.pdf"),
  plot = p_violin,
  width = 5.5,
  height = 5.5
)

expr_for_deg <- expr_mat2
rownames(expr_for_deg) <- expr_for_deg$GeneSymbol
expr_for_deg <- expr_for_deg[, -1, drop = FALSE]

expr_for_deg <- as.matrix(expr_for_deg)
mode(expr_for_deg) <- "numeric"

sample_names <- colnames(expr_for_deg)
sample_cluster <- sub("^.*_(C1|C2)$", "\\1", sample_names)
group <- factor(sample_cluster, levels = c("C1", "C2"))

if (length(unique(group)) < 2) {
  stop("Need both C1 and C2 for DEG analysis.")
}

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_for_deg, design)
contrast.matrix <- makeContrasts(C2 - C1, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg_res <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg_res$GeneSymbol <- rownames(deg_res)
deg_res <- deg_res[, c("GeneSymbol", setdiff(colnames(deg_res), "GeneSymbol"))]

write.csv(
  deg_res,
  file = file.path(out_dir, "GSE96058_C2_vs_C1_limma_results.csv"),
  row.names = FALSE
)

p_cutoff <- 0.05
logfc_cutoff <- 1

deg_res$group <- dplyr::case_when(
  deg_res$logFC > logfc_cutoff & deg_res$adj.P.Val < p_cutoff ~ "up",
  deg_res$logFC < -logfc_cutoff & deg_res$adj.P.Val < p_cutoff ~ "down",
  TRUE ~ "none"
)

deg_up <- deg_res %>%
  filter(group == "up") %>%
  arrange(desc(logFC))

deg_down <- deg_res %>%
  filter(group == "down") %>%
  arrange(logFC)

deg_sig <- deg_res %>%
  filter(group != "none")

write.table(
  deg_up$GeneSymbol,
  file = file.path(out_dir, "GSE96058_C2_vs_C1_up_genes.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  deg_down$GeneSymbol,
  file = file.path(out_dir, "GSE96058_C2_vs_C1_down_genes.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  deg_sig$GeneSymbol,
  file = file.path(out_dir, "GSE96058_C2_vs_C1_sig_genes.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.csv(
  deg_sig,
  file = file.path(out_dir, "GSE96058_C2_vs_C1_sig_DEG_with_logFC.csv"),
  row.names = FALSE
)

vol_df <- deg_res
vol_df$pvalue        <- vol_df$P.Value
vol_df$log2FoldChange <- vol_df$logFC

vol_df$group <- dplyr::case_when(
  vol_df$log2FoldChange >  logfc_cutoff & vol_df$adj.P.Val < p_cutoff ~ "up",
  vol_df$log2FoldChange < -logfc_cutoff & vol_df$adj.P.Val < p_cutoff ~ "down",
  TRUE ~ "none"
)

vol_df$neglog10p <- -log10(vol_df$adj.P.Val)
vol_df$group <- factor(vol_df$group, levels = c("up", "down", "none"))

mycol <- c("#b43665", "#6fa6cf", "#d8d8d8")

fc_abs_max <- ceiling(max(abs(vol_df$log2FoldChange), na.rm = TRUE))
fc_abs_max <- max(fc_abs_max, 2)
fc_breaks  <- pretty(c(-fc_abs_max, fc_abs_max), n = 5)

neglogp_max <- ceiling(max(vol_df$neglog10p, na.rm = TRUE))
neglogp_max <- max(neglogp_max, 5)

neglogp_breaks <- pretty(c(0, neglogp_max), n = 5)

p_vol <- ggplot(
  data = vol_df,
  aes(x = log2FoldChange, y = neglog10p, color = group)
) +
  geom_point(size = 1.8, alpha = 0.65) +
  scale_colour_manual(
    name   = "",
    values = scales::alpha(mycol, 0.75),
    labels = c(
      paste0("Up (", sum(vol_df$group == "up"), ")"),
      paste0("Down (", sum(vol_df$group == "down"), ")"),
      "NS"
    )
  ) +

  scale_x_continuous(
    limits = c(-fc_abs_max, fc_abs_max),
    breaks = fc_breaks
  ) +

  scale_y_continuous(
    limits = c(0, neglogp_max * 1.05),
    breaks = neglogp_breaks,
    expand = expansion(mult = c(0, 0.05))
  ) +

  geom_vline(
    xintercept = c(-logfc_cutoff, logfc_cutoff),
    linewidth  = 0.6,
    color      = "black",
    linetype   = "dashed"
  ) +

  geom_hline(
    yintercept = -log10(p_cutoff),
    linewidth  = 0.6,
    color      = "black",
    linetype   = "dashed"
  ) +
  labs(
    title = "Volcano plot: C2 vs C1",
    x     = expression(log[2]~"Fold Change"),
    y     = expression(-log[10]~"(FDR)")
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title        = element_text(size = 13, face = "bold", color = "black"),
    axis.text         = element_text(size = 12, color = "black"),
    axis.ticks        = element_line(color = "black"),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position   = c(0.88, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.4),
    legend.text       = element_text(size = 11)
  ) +
  coord_flip()

ggsave(
  filename = file.path(out_dir, "GSE96058_C2_vs_C1_flipped_volcano.pdf"),
  plot     = p_vol,
  width    = 6,
  height   = 4
)

summary_df <- data.frame(
  Matched_Expression_Samples = ncol(expr_for_deg),
  C1_Samples = sum(group == "C1"),
  C2_Samples = sum(group == "C2"),
  Total_Genes = nrow(expr_for_deg),
  Significant_Up = nrow(deg_up),
  Significant_Down = nrow(deg_down),
  P_cutoff = p_cutoff,
  logFC_cutoff = logfc_cutoff,
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  file = file.path(out_dir, "GSE96058_C2_vs_C1_summary.csv"),
  row.names = FALSE
)

gc()

suppressPackageStartupMessages({
  library(VennDiagram)
  library(RColorBrewer)
  library(grid)
})

workDir <- file.path(PROJECT_DIR, "results/bulk/discovery")

up_file <- file.path(PROJECT_DIR, "results/bulk/discovery/04_C1C2_DEG_Analysis/GSE96058_C2_vs_C1_up_genes.txt")
hvg_file <- file.path(PROJECT_DIR, "results/bulk/discovery/03_HVG_blue_module_gene_list.txt")

out_dir <- file.path(workDir, "05_gene")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

c2_up_genes <- read.table(up_file, header = FALSE, stringsAsFactors = FALSE)$V1
hvg_genes   <- read.table(hvg_file, header = FALSE, stringsAsFactors = FALSE)$V1

c2_up_genes <- unique(trimws(as.character(c2_up_genes)))
hvg_genes   <- unique(trimws(as.character(hvg_genes)))

c2_up_genes <- c2_up_genes[!is.na(c2_up_genes) & c2_up_genes != ""]
hvg_genes   <- hvg_genes[!is.na(hvg_genes) & hvg_genes != ""]

intersection_genes <- intersect(c2_up_genes, hvg_genes)

write.table(
  intersection_genes,
  file = file.path(out_dir, "C2_up_HVG_intersection_genes.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

venn.plot <- venn.diagram(
  x = list(
    `C2 up genes` = c2_up_genes,
    `HVG genes`   = hvg_genes
  ),
  filename = NULL,
  category.names = c("C2 up genes", "HVG genes"),
  euler.d = FALSE,
  scaled = FALSE,
  fill = c("#b43665", "#6fa6cf"),
  alpha = 0.5,
  lwd = 3,
  lty = "solid",
  col = c("#b43665", "#6fa6cf"),
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  cat.fontfamily = "sans",
  cat.col = c("#b43665", "#6fa6cf"),
  print.mode = "raw",
  margin = 0.08
)

pdf(
  file = file.path(out_dir, "C2_up_HVG_intersection_venn.pdf"),
  width = 8,
  height = 8
)
grid.draw(venn.plot)
dev.off()

gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GSEABase)
  library(ConsensusClusterPlus)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(vegan)
  library(reshape2)
  library(grid)
  library(gridExtra)
})

workDir  <- file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538")
exprFile <- file.path(PROJECT_DIR, "results/bulk/GSE81538/GSE81538_expr_HighLow.csv")
geneFile <- file.path(PROJECT_DIR, "results/bulk/discovery/03_HVG_blue_module_gene_list.txt")

dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
setwd(workDir)
dir.create("01_Result", showWarnings = FALSE)
dir.create("02_Figures", showWarnings = FALSE)
dir.create("03_Tables", showWarnings = FALSE)

col_c1 <- "#b43665"
col_c2 <- "#6fa6cf"
cluster_colors_2 <- c("1" = col_c1, "2" = col_c2)
cluster_colors_named <- c("C1" = col_c1, "C2" = col_c2)

hm_colors <- colorRampPalette(c("white", "#f4d6e2", col_c1))(100)

maxK       <- 9
reps_num   <- 100
pItem_num  <- 0.8
pFeature_num <- 1
seed_num   <- 2025

expr_raw <- read.csv(exprFile,
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

expr_mat <- as.matrix(expr_raw)
mode(expr_mat) <- "numeric"

gene_list <- read.table(geneFile,
                        header = FALSE,
                        stringsAsFactors = FALSE,
                        sep = "\t",
                        quote = "",
                        fill = TRUE)

rloop_genes <- unique(na.omit(as.character(gene_list[, 1])))
rloop_genes <- rloop_genes[rloop_genes != ""]

matched_genes <- intersect(rloop_genes, rownames(expr_mat))
missing_genes <- setdiff(rloop_genes, rownames(expr_mat))

write.table(data.frame(Missing_Genes = missing_genes),
            file = "03_Tables/Missing_Rloop_genes.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

expr_rloop <- expr_mat[matched_genes, , drop = FALSE]

gene_var <- apply(expr_rloop, 1, var, na.rm = TRUE)
zero_var_genes <- names(gene_var[gene_var == 0])

if (length(zero_var_genes) > 0) {
  expr_rloop <- expr_rloop[!(rownames(expr_rloop) %in% zero_var_genes), , drop = FALSE]
}

write.csv(expr_rloop, file = "03_Tables/Rloop_expr_used_for_clustering.csv", quote = FALSE)

results <- ConsensusClusterPlus(
  d = expr_rloop,
  maxK = maxK,
  reps = reps_num,
  pItem = pItem_num,
  pFeature = pFeature_num,
  title = file.path(workDir, "01_Result"),
  clusterAlg = "pam",
  distance = "euclidean",
  seed = seed_num,
  plot = "pdf"
)

Kvec <- 2:maxK
PAC <- sapply(Kvec, function(k) {
  M <- results[[k]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  Fn(0.9) - Fn(0.1)
})

pac_df <- data.frame(
  K = Kvec,
  PAC = PAC
)

write.csv(pac_df, file = "03_Tables/PAC_values.csv", row.names = FALSE)

optK_pac <- Kvec[which.min(PAC)]

optK_val <- 2

cluster_num <- results[[optK_val]]$consensusClass
cluster_df <- data.frame(
  sample = names(cluster_num),
  Cluster_num = cluster_num,
  Cluster = paste0("C", cluster_num),
  stringsAsFactors = FALSE
)

write.csv(cluster_df,
          file = sprintf("03_Tables/ConsensusCluster_K%d_assignment.csv", optK_val),
          row.names = FALSE)

print(table(cluster_df$Cluster))

p_pac <- ggplot(pac_df, aes(x = K, y = PAC)) +
  geom_line(color = "#444444", linewidth = 1.1) +
  geom_point(size = 3.2, color = col_c1) +
  geom_vline(xintercept = optK_val, linetype = "solid", color = col_c1, linewidth = 0.9) +
  geom_vline(xintercept = optK_pac, linetype = "dashed", color = col_c2, linewidth = 0.9) +
  annotate("text", x = optK_val, y = max(pac_df$PAC) * 0.96,
           label = paste0("Selected K = ", optK_val),
           color = col_c1, hjust = -0.1, vjust = 1, size = 4.2, fontface = "bold") +
  annotate("text", x = optK_pac, y = max(pac_df$PAC) * 0.84,
           label = paste0("Min PAC K = ", optK_pac),
           color = col_c2, hjust = -0.1, vjust = 1, size = 4.0, fontface = "bold") +
  scale_x_continuous(breaks = Kvec) +
  labs(title = "PAC Curve for Consensus Clustering",
       x = "Number of clusters (K)",
       y = "Proportion of ambiguous clustering") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black")
  )

ggsave("02_Figures/Figure_S1_PAC_curve.pdf", p_pac, width = 7.5, height = 5.8)

cdf_data <- data.frame()
for (k in 2:maxK) {
  m <- results[[k]]$consensusMatrix
  vals <- m[lower.tri(m, diag = FALSE)]
  tmp <- data.frame(
    K = factor(rep(k, length(vals))),
    Consensus = vals
  )
  cdf_data <- rbind(cdf_data, tmp)
}

p_cdf <- ggplot(cdf_data, aes(x = Consensus, color = K)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(title = "Consensus Cumulative Distribution Function",
       x = "Consensus index",
       y = "Cumulative distribution",
       color = "K") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black")
  ) +
  scale_color_brewer(palette = "Set1")

ggsave("02_Figures/Figure_S2_CDF_curve.pdf", p_cdf, width = 8, height = 6)

consensus_mat <- results[[optK_val]]$consensusMatrix
colnames(consensus_mat) <- names(cluster_num)
rownames(consensus_mat) <- names(cluster_num)

ann_col <- data.frame(Cluster = factor(paste0("C", cluster_num), levels = c("C1", "C2")))
rownames(ann_col) <- names(cluster_num)

ann_colors <- list(Cluster = cluster_colors_named)

pdf("02_Figures/Figure_3A_Consensus_heatmap_K2.pdf", width = 9.2, height = 8.5)
pheatmap(
  consensus_mat,
  color = colorRampPalette(c("white", "#d9e7f3", col_c2, col_c1))(100),
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  main = "Consensus Matrix Heatmap (K = 2)",
  fontsize = 12
)
dev.off()

cluster_tab <- do.call(cbind, lapply(2:maxK, function(k) {
  results[[k]]$consensusClass
}))
colnames(cluster_tab) <- paste0("K", 2:maxK)

final_cluster <- cluster_tab[, paste0("K", optK_val)]
cluster_tab <- cluster_tab[order(final_cluster), , drop = FALSE]

melt_data <- reshape2::melt(cluster_tab)
colnames(melt_data) <- c("Sample", "K", "Cluster")

nclust_max <- max(cluster_tab)
pal_dynamic <- colorRampPalette(c(col_c2, "#d9d9d9", col_c1))(max(6, nclust_max))

p_track <- ggplot(melt_data, aes(x = K, y = Sample, fill = factor(Cluster))) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = pal_dynamic, name = "Cluster") +
  labs(title = paste0("Cluster Membership Tracking Across K (Selected K = ", optK_val, ")"),
       x = "Number of clusters",
       y = "Samples") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black")
  )

ggsave("02_Figures/Figure_S3_Tracking_plot.pdf", p_track, width = 8.8, height = 7.6)

pca_res <- prcomp(t(expr_rloop), center = TRUE, scale. = FALSE)
pca_data <- as.data.frame(pca_res$x[, 1:2])
pca_data$sample <- rownames(pca_data)

pca_data <- pca_data %>%
  left_join(cluster_df[, c("sample", "Cluster")], by = "sample")

percent_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
cum_var <- sum(percent_var)

adonis_result <- adonis2(dist(pca_data[, c("PC1", "PC2")]) ~ Cluster, data = pca_data, permutations = 999)
p_perm <- adonis_result$`Pr(>F)`[1]

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3.2, alpha = 0.88) +
  stat_ellipse(level = 0.95, linewidth = 0.9, aes(fill = Cluster), geom = "polygon", alpha = 0.12, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors_named) +
  scale_fill_manual(values = cluster_colors_named) +
  labs(title = "Principal Component Analysis of R-loop Expression Pattern",
       subtitle = paste0("PERMANOVA p = ", format.pval(p_perm, digits = 3),
                         " | Total variance explained = ", cum_var, "%"),
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", color = "black")
  )

ggsave("02_Figures/Figure_3B_PCA_K2.pdf", p_pca, width = 7.5, height = 6.2)

gene_var2 <- apply(expr_rloop, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_var2, decreasing = TRUE))[1:min(100, length(gene_var2))]

expr_top <- expr_rloop[top_genes, cluster_df$sample, drop = FALSE]
expr_top_scaled <- t(scale(t(expr_top)))
expr_top_scaled[is.na(expr_top_scaled)] <- 0

ann_col2 <- data.frame(Cluster = factor(cluster_df$Cluster, levels = c("C1", "C2")))
rownames(ann_col2) <- cluster_df$sample

pdf("02_Figures/Figure_3C_TopVariableGenes_heatmap.pdf", width = 8.8, height = 9.5)
pheatmap(
  expr_top_scaled,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  annotation_col = ann_col2,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  border_color = NA,
  fontsize = 11,
  main = "Top Variable Genes Across Clusters"
)
dev.off()

gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

workDir <- file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538")
setwd(workDir)
assign_file <- file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538/03_Tables/ConsensusCluster_K2_assignment.csv")
expr_file   <- file.path(PROJECT_DIR, "results/bulk/GSE81538/GSE81538_expr_HighLow.csv")
clin_file   <- file.path(PROJECT_DIR, "results/bulk/GSE81538/GSE81538_HRp_HERn_clinical.csv")

out_dir <- file.path(workDir, "C1C2_Rloop_Violin")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assign_df <- read.csv(assign_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "Cluster") %in% colnames(assign_df))) {
  stop("Assignment file must contain columns: sample, Cluster")
}

assign_df$sample  <- trimws(as.character(assign_df$sample))
assign_df$Cluster <- trimws(as.character(assign_df$Cluster))

assign_df$sample_base <- sub("_[^_]+$", "", assign_df$sample)

assign_df$new_sample <- paste0(assign_df$sample_base, "_", assign_df$Cluster)

assign_out <- assign_df[, c("new_sample", "Cluster_num", "Cluster", "sample", "sample_base")]
colnames(assign_out) <- c("sample", "Cluster_num", "Cluster", "old_sample", "sample_base")

write.csv(
  assign_out,
  file = file.path(out_dir, "ConsensusCluster_K2_GSE81538.csv"),
  row.names = FALSE
)

expr_df <- read.csv(expr_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

colnames(expr_df)[1] <- "GeneSymbol"

old_expr_samples <- colnames(expr_df)[-1]

map_df <- assign_df[, c("sample", "new_sample")]
sample_map <- setNames(map_df$new_sample, map_df$sample)

matched_new_names <- sample_map[old_expr_samples]
keep_cols <- !is.na(matched_new_names)

expr_mat2 <- expr_df[, c(TRUE, keep_cols), drop = FALSE]
colnames(expr_mat2) <- c("GeneSymbol", unname(matched_new_names[keep_cols]))

write.csv(
  expr_mat2,
  file = file.path(out_dir, "GSE81538_expr_C1C2.csv"),
  row.names = FALSE
)

clinical <- read.csv(clin_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "Rloop_ssGSEA_Score") %in% colnames(clinical))) {
  stop("Clinical file must contain columns: sample, Rloop_ssGSEA_Score")
}

clinical$sample <- trimws(as.character(clinical$sample))

clinical2 <- merge(
  clinical,
  assign_df[, c("sample_base", "Cluster")],
  by.x = "sample",
  by.y = "sample_base",
  all.x = FALSE,
  all.y = FALSE
)

plot_df <- clinical2 %>%
  mutate(
    Cluster = factor(Cluster, levels = c("C1", "C2")),
    Rloop_ssGSEA_Score = as.numeric(Rloop_ssGSEA_Score)
  ) %>%
  filter(!is.na(Cluster), !is.na(Rloop_ssGSEA_Score))

p_violin <- ggplot(plot_df, aes(x = Cluster, y = Rloop_ssGSEA_Score, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", color = "black",
               outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.4, color = "grey25") +
  scale_fill_manual(values = c("C1" = "#187d79", "C2" = "#e58027")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 5
  ) +
  labs(
    x = "Consensus Cluster",
    y = "R-loop ssGSEA Score",
    title = "R-loop score difference between C1 and C2"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = file.path(out_dir, "GSE81538_Rloop_C1C2_violin.pdf"),
  plot = p_violin,
  width = 5.5,
  height = 5.5
)

gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(GSEABase)
  library(ConsensusClusterPlus)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(vegan)
  library(reshape2)
  library(grid)
  library(gridExtra)
})

workDir  <- file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC")
exprFile <- file.path(PROJECT_DIR, "results/bulk/METBRIC/METBRIC_expr_HighLow.csv")
geneFile <- file.path(PROJECT_DIR, "results/bulk/discovery/03_HVG_blue_module_gene_list.txt")

dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
setwd(workDir)
dir.create("01_Result", showWarnings = FALSE)
dir.create("02_Figures", showWarnings = FALSE)
dir.create("03_Tables", showWarnings = FALSE)

col_c1 <- "#b43665"
col_c2 <- "#6fa6cf"
cluster_colors_2 <- c("1" = col_c1, "2" = col_c2)
cluster_colors_named <- c("C1" = col_c1, "C2" = col_c2)

hm_colors <- colorRampPalette(c("white", "#f4d6e2", col_c1))(100)

maxK       <- 9
reps_num   <- 100
pItem_num  <- 0.8
pFeature_num <- 1
seed_num   <- 2025

expr_raw <- read.csv(exprFile,
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

expr_mat <- as.matrix(expr_raw)
mode(expr_mat) <- "numeric"

gene_list <- read.table(geneFile,
                        header = FALSE,
                        stringsAsFactors = FALSE,
                        sep = "\t",
                        quote = "",
                        fill = TRUE)

rloop_genes <- unique(na.omit(as.character(gene_list[, 1])))
rloop_genes <- rloop_genes[rloop_genes != ""]

matched_genes <- intersect(rloop_genes, rownames(expr_mat))
missing_genes <- setdiff(rloop_genes, rownames(expr_mat))

write.table(data.frame(Missing_Genes = missing_genes),
            file = "03_Tables/Missing_Rloop_genes.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

expr_rloop <- expr_mat[matched_genes, , drop = FALSE]

gene_var <- apply(expr_rloop, 1, var, na.rm = TRUE)
zero_var_genes <- names(gene_var[gene_var == 0])

if (length(zero_var_genes) > 0) {
  expr_rloop <- expr_rloop[!(rownames(expr_rloop) %in% zero_var_genes), , drop = FALSE]
}

write.csv(expr_rloop, file = "03_Tables/Rloop_expr_used_for_clustering.csv", quote = FALSE)

results <- ConsensusClusterPlus(
  d = expr_rloop,
  maxK = maxK,
  reps = reps_num,
  pItem = pItem_num,
  pFeature = pFeature_num,
  title = file.path(workDir, "01_Result"),
  clusterAlg = "pam",
  distance = "euclidean",
  seed = seed_num,
  plot = "pdf"
)

Kvec <- 2:maxK
PAC <- sapply(Kvec, function(k) {
  M <- results[[k]]$consensusMatrix
  Fn <- ecdf(M[lower.tri(M)])
  Fn(0.9) - Fn(0.1)
})

pac_df <- data.frame(
  K = Kvec,
  PAC = PAC
)

write.csv(pac_df, file = "03_Tables/PAC_values.csv", row.names = FALSE)

optK_pac <- Kvec[which.min(PAC)]

optK_val <- 2

cluster_num <- results[[optK_val]]$consensusClass
cluster_df <- data.frame(
  sample = names(cluster_num),
  Cluster_num = cluster_num,
  Cluster = paste0("C", cluster_num),
  stringsAsFactors = FALSE
)

write.csv(cluster_df,
          file = sprintf("03_Tables/ConsensusCluster_K%d_assignment.csv", optK_val),
          row.names = FALSE)

print(table(cluster_df$Cluster))

p_pac <- ggplot(pac_df, aes(x = K, y = PAC)) +
  geom_line(color = "#444444", linewidth = 1.1) +
  geom_point(size = 3.2, color = col_c1) +
  geom_vline(xintercept = optK_val, linetype = "solid", color = col_c1, linewidth = 0.9) +
  geom_vline(xintercept = optK_pac, linetype = "dashed", color = col_c2, linewidth = 0.9) +
  annotate("text", x = optK_val, y = max(pac_df$PAC) * 0.96,
           label = paste0("Selected K = ", optK_val),
           color = col_c1, hjust = -0.1, vjust = 1, size = 4.2, fontface = "bold") +
  annotate("text", x = optK_pac, y = max(pac_df$PAC) * 0.84,
           label = paste0("Min PAC K = ", optK_pac),
           color = col_c2, hjust = -0.1, vjust = 1, size = 4.0, fontface = "bold") +
  scale_x_continuous(breaks = Kvec) +
  labs(title = "PAC Curve for Consensus Clustering",
       x = "Number of clusters (K)",
       y = "Proportion of ambiguous clustering") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black")
  )

ggsave("02_Figures/Figure_S1_PAC_curve.pdf", p_pac, width = 7.5, height = 5.8)

cdf_data <- data.frame()
for (k in 2:maxK) {
  m <- results[[k]]$consensusMatrix
  vals <- m[lower.tri(m, diag = FALSE)]
  tmp <- data.frame(
    K = factor(rep(k, length(vals))),
    Consensus = vals
  )
  cdf_data <- rbind(cdf_data, tmp)
}

p_cdf <- ggplot(cdf_data, aes(x = Consensus, color = K)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  labs(title = "Consensus Cumulative Distribution Function",
       x = "Consensus index",
       y = "Cumulative distribution",
       color = "K") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black")
  ) +
  scale_color_brewer(palette = "Set1")

ggsave("02_Figures/Figure_S2_CDF_curve.pdf", p_cdf, width = 8, height = 6)

consensus_mat <- results[[optK_val]]$consensusMatrix
colnames(consensus_mat) <- names(cluster_num)
rownames(consensus_mat) <- names(cluster_num)

ann_col <- data.frame(Cluster = factor(paste0("C", cluster_num), levels = c("C1", "C2")))
rownames(ann_col) <- names(cluster_num)

ann_colors <- list(Cluster = cluster_colors_named)

pdf("02_Figures/Figure_3A_Consensus_heatmap_K2.pdf", width = 9.2, height = 8.5)
pheatmap(
  consensus_mat,
  color = colorRampPalette(c("white", "#d9e7f3", col_c2, col_c1))(100),
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  main = "Consensus Matrix Heatmap (K = 2)",
  fontsize = 12
)
dev.off()

cluster_tab <- do.call(cbind, lapply(2:maxK, function(k) {
  results[[k]]$consensusClass
}))
colnames(cluster_tab) <- paste0("K", 2:maxK)

final_cluster <- cluster_tab[, paste0("K", optK_val)]
cluster_tab <- cluster_tab[order(final_cluster), , drop = FALSE]

melt_data <- reshape2::melt(cluster_tab)
colnames(melt_data) <- c("Sample", "K", "Cluster")

nclust_max <- max(cluster_tab)
pal_dynamic <- colorRampPalette(c(col_c2, "#d9d9d9", col_c1))(max(6, nclust_max))

p_track <- ggplot(melt_data, aes(x = K, y = Sample, fill = factor(Cluster))) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = pal_dynamic, name = "Cluster") +
  labs(title = paste0("Cluster Membership Tracking Across K (Selected K = ", optK_val, ")"),
       x = "Number of clusters",
       y = "Samples") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black")
  )

ggsave("02_Figures/Figure_S3_Tracking_plot.pdf", p_track, width = 8.8, height = 7.6)

pca_res <- prcomp(t(expr_rloop), center = TRUE, scale. = FALSE)
pca_data <- as.data.frame(pca_res$x[, 1:2])
pca_data$sample <- rownames(pca_data)

pca_data <- pca_data %>%
  left_join(cluster_df[, c("sample", "Cluster")], by = "sample")

percent_var <- round(100 * summary(pca_res)$importance[2, 1:2], 1)
cum_var <- sum(percent_var)

adonis_result <- adonis2(dist(pca_data[, c("PC1", "PC2")]) ~ Cluster, data = pca_data, permutations = 999)
p_perm <- adonis_result$`Pr(>F)`[1]

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3.2, alpha = 0.88) +
  stat_ellipse(level = 0.95, linewidth = 0.9, aes(fill = Cluster), geom = "polygon", alpha = 0.12, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors_named) +
  scale_fill_manual(values = cluster_colors_named) +
  labs(title = "Principal Component Analysis of R-loop Expression Pattern",
       subtitle = paste0("PERMANOVA p = ", format.pval(p_perm, digits = 3),
                         " | Total variance explained = ", cum_var, "%"),
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", color = "black")
  )

ggsave("02_Figures/Figure_3B_PCA_K2.pdf", p_pca, width = 7.5, height = 6.2)

gene_var2 <- apply(expr_rloop, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_var2, decreasing = TRUE))[1:min(100, length(gene_var2))]

expr_top <- expr_rloop[top_genes, cluster_df$sample, drop = FALSE]
expr_top_scaled <- t(scale(t(expr_top)))
expr_top_scaled[is.na(expr_top_scaled)] <- 0

ann_col2 <- data.frame(Cluster = factor(cluster_df$Cluster, levels = c("C1", "C2")))
rownames(ann_col2) <- cluster_df$sample

pdf("02_Figures/Figure_3C_TopVariableGenes_heatmap.pdf", width = 8.8, height = 9.5)
pheatmap(
  expr_top_scaled,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  annotation_col = ann_col2,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  border_color = NA,
  fontsize = 11,
  main = "Top Variable Genes Across Clusters"
)
dev.off()

gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

workDir <- file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC")
setwd(workDir)
assign_file <- file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC/03_Tables/ConsensusCluster_K2_assignment.csv")
expr_file   <- file.path(PROJECT_DIR, "results/bulk/METBRIC/METBRIC_expr_HighLow.csv")
clin_file   <- file.path(PROJECT_DIR, "results/bulk/METBRIC/05_METBRIC_HRp_HERn_clinical.csv")

out_dir <- file.path(workDir, "C1C2_Rloop_Violin")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

assign_df <- read.csv(assign_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "Cluster") %in% colnames(assign_df))) {
  stop("Assignment file must contain columns: sample, Cluster")
}

assign_df$sample  <- trimws(as.character(assign_df$sample))
assign_df$Cluster <- trimws(as.character(assign_df$Cluster))

assign_df$sample_base <- sub("_[^_]+$", "", assign_df$sample)

assign_df$new_sample <- paste0(assign_df$sample_base, "_", assign_df$Cluster)

assign_out <- assign_df[, c("new_sample", "Cluster_num", "Cluster", "sample", "sample_base")]
colnames(assign_out) <- c("sample", "Cluster_num", "Cluster", "old_sample", "sample_base")

write.csv(
  assign_out,
  file = file.path(out_dir, "ConsensusCluster_K2_METABRIC.csv"),
  row.names = FALSE
)

expr_df <- read.csv(expr_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

colnames(expr_df)[1] <- "GeneSymbol"

old_expr_samples <- colnames(expr_df)[-1]

map_df <- assign_df[, c("sample", "new_sample")]
sample_map <- setNames(map_df$new_sample, map_df$sample)

matched_new_names <- sample_map[old_expr_samples]
keep_cols <- !is.na(matched_new_names)

expr_mat2 <- expr_df[, c(TRUE, keep_cols), drop = FALSE]
colnames(expr_mat2) <- c("GeneSymbol", unname(matched_new_names[keep_cols]))

write.csv(
  expr_mat2,
  file = file.path(out_dir, "METABRIC_expr_C1C2.csv"),
  row.names = FALSE
)

clinical <- read.csv(clin_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample", "Rloop_ssGSEA_Score") %in% colnames(clinical))) {
  stop("Clinical file must contain columns: sample, Rloop_ssGSEA_Score")
}

clinical$sample <- trimws(as.character(clinical$sample))

clinical2 <- merge(
  clinical,
  assign_df[, c("sample_base", "Cluster")],
  by.x = "sample",
  by.y = "sample_base",
  all.x = FALSE,
  all.y = FALSE
)

plot_df <- clinical2 %>%
  mutate(
    Cluster = factor(Cluster, levels = c("C1", "C2")),
    Rloop_ssGSEA_Score = as.numeric(Rloop_ssGSEA_Score)
  ) %>%
  filter(!is.na(Cluster), !is.na(Rloop_ssGSEA_Score))

p_violin <- ggplot(plot_df, aes(x = Cluster, y = Rloop_ssGSEA_Score, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "white", color = "black",
               outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.4, color = "grey25") +
  scale_fill_manual(values = c("C1" = "#187d79", "C2" = "#e58027")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 5
  ) +
  labs(
    x = "Consensus Cluster",
    y = "R-loop ssGSEA Score",
    title = "R-loop score difference between C1 and C2"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  filename = file.path(out_dir, "METABRIC_Rloop_C1C2_violin.pdf"),
  plot = p_violin,
  width = 5.5,
  height = 5.5
)

suppressPackageStartupMessages({
  library(GSEABase)
  library(ggvenn)
  library(ggplot2)
})

RLOOP_GMT <- file.path(PROJECT_DIR, "results/bulk/GSE96058/01_Rloop_regulators.gmt")
RL65_TXT  <- file.path(PROJECT_DIR, "results/bulk/discovery/05_gene/Final_validated_Rloop_signature_genes.txt")
OUTDIR    <- file.path(PROJECT_DIR, "results/bulk/RLSig65")

rloop_genes <- unique(unlist(lapply(
  as.list(getGmt(RLOOP_GMT, geneIdType = SymbolIdentifier())), geneIds
)))

rl65_genes <- unique(trimws(readLines(RL65_TXT, warn = FALSE)))
rl65_genes <- rl65_genes[rl65_genes != ""]

overlap_genes <- intersect(rloop_genes, rl65_genes)

writeLines(sort(overlap_genes),
           file.path(OUTDIR, "Rloop_RLsig65_overlap_genes.txt"))

p <- ggvenn(
  list("Original R-loop" = rloop_genes, "RL-Sig65" = rl65_genes),
  fill_color    = c("#4DBBD5", "#E64B35"),
  fill_alpha    = 0.55,
  stroke_size   = 0.7,
  stroke_color  = "white",
  set_name_size = 5,
  text_size     = 5,
  show_percentage = FALSE
) +
  labs(title    = "Overlap between Original R-loop regulators and RL-Sig65",
       subtitle = sprintf("R-loop: %d  |  RL-Sig65: %d  |  Overlap: %d",
                          length(rloop_genes), length(rl65_genes), length(overlap_genes))) +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(OUTDIR, "Figure_Rloop_vs_RLsig65_Venn.pdf"),
       p, width = 6, height = 5.5)

gc()

suppressPackageStartupMessages({
  library(dplyr)
})

file1 <- file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538/C1C2_Rloop_Violin/GSE81538_expr_C1C2.csv")
out1  <- file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538/C1C2_Rloop_Violin/GSE81538_expr_C1C2_HL.csv")

file2 <- file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC/C1C2_Rloop_Violin/METABRIC_expr_C1C2.csv")
out2  <- file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC/C1C2_Rloop_Violin/METABRIC_expr_C1C2_HL.csv")

file3 <- file.path(PROJECT_DIR, "results/bulk/discovery/04_C1C2_DEG_Analysis/GSE96058_expr_C1C2.csv")
out3  <- file.path(PROJECT_DIR, "results/bulk/discovery/04_C1C2_DEG_Analysis/GSE96058_expr_C1C2_LH.csv")

rename_suffix <- function(df, old1, new1, old2, new2) {
  colnames(df) <- gsub(paste0("_", old1, "$"), paste0("_", new1), colnames(df))
  colnames(df) <- gsub(paste0("_", old2, "$"), paste0("_", new2), colnames(df))
  return(df)
}

df1 <- read.csv(file1, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
df1 <- rename_suffix(df1, old1 = "C1", new1 = "Rloop-high",
                     old2 = "C2", new2 = "Rloop-low")
write.csv(df1, out1, row.names = FALSE, quote = FALSE)

df2 <- read.csv(file2, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
df2 <- rename_suffix(df2, old1 = "C1", new1 = "Rloop-high",
                     old2 = "C2", new2 = "Rloop-low")
write.csv(df2, out2, row.names = FALSE, quote = FALSE)

df3 <- read.csv(file3, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
df3 <- rename_suffix(df3, old1 = "C1", new1 = "Rloop-low",
                     old2 = "C2", new2 = "Rloop-high")
write.csv(df3, out3, row.names = FALSE, quote = FALSE)
