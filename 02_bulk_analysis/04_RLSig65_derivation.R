# RL-Sig65 derivation and cross-dataset validation

PROJECT_DIR <- "."

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(GSVA)
  library(GSEABase)
  library(patchwork)
  library(ggExtra)
  library(gridExtra)
  library(grid)
  library(scales)
})

expr_files <- c(
  file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538/C1C2_Rloop_Violin/GSE81538_expr_C1C2.csv"),
  file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC/C1C2_Rloop_Violin/METABRIC_expr_C1C2.csv")
)

clin_files <- c(
  file.path(PROJECT_DIR, "results/bulk/GSE81538/GSE81538_HRp_HERn_clinical.csv"),
  file.path(PROJECT_DIR, "results/bulk/METBRIC/05_METBRIC_HRp_HERn_clinical.csv")
)

assign_files <- c(
  file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538/03_Tables/ConsensusCluster_K2_assignment.csv"),
  file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC/03_Tables/ConsensusCluster_K2_assignment.csv")
)

gene_file <- file.path(PROJECT_DIR, "results/bulk/discovery/05_gene/C2_up_HVG_intersection_genes.txt")

outdir <- file.path(PROJECT_DIR, "results/bulk/discovery/05_gene")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

candidate_genes <- fread(gene_file, header = FALSE, data.table = FALSE)[, 1]
candidate_genes <- unique(candidate_genes)

calc_validation <- function(expr_file, clin_file, assign_file, candidate_genes) {

  expr <- fread(expr_file, data.table = FALSE, check.names = FALSE)
  rownames(expr) <- expr[, 1]
  expr <- expr[, -1, drop = FALSE]

  clinical <- fread(clin_file, data.table = FALSE)
  clinical$sample <- trimws(as.character(clinical$sample))

  if (!"Rloop_ssGSEA_Score" %in% colnames(clinical)) {
    stop(paste("missingRloop_ssGSEA_Score:", clin_file))
  }

  assign_df <- fread(assign_file, data.table = FALSE)
  assign_df$sample_base <- sub("_[^_]+$", "", assign_df$sample)

  clin_cluster <- merge(
    clinical,
    assign_df[, c("sample_base", "Cluster")],
    by.x = "sample",
    by.y = "sample_base",
    all = FALSE
  )

  c1_score <- mean(clin_cluster$Rloop_ssGSEA_Score[clin_cluster$Cluster == "C1"], na.rm = TRUE)
  c2_score <- mean(clin_cluster$Rloop_ssGSEA_Score[clin_cluster$Cluster == "C2"], na.rm = TRUE)

  if (c1_score > c2_score) {
    high_cluster <- "C1"
    low_cluster  <- "C2"
  } else {
    high_cluster <- "C2"
    low_cluster  <- "C1"
  }

  common_genes <- intersect(candidate_genes, rownames(expr))
  expr_sub <- expr[common_genes, , drop = FALSE]

  sample_names <- colnames(expr_sub)
  high_samples <- grep(paste0("_", high_cluster, "$"), sample_names, value = TRUE)
  low_samples  <- grep(paste0("_", low_cluster, "$"), sample_names, value = TRUE)

  if (length(high_samples) == 0 || length(low_samples) == 0) {
    stop(paste("groupsample:", expr_file))
  }

  high_mean <- rowMeans(expr_sub[, high_samples, drop = FALSE], na.rm = TRUE)
  low_mean  <- rowMeans(expr_sub[, low_samples, drop = FALSE], na.rm = TRUE)
  MeanDiff <- high_mean - low_mean

  sample_base_vec <- sub("_(C1|C2)$", "", sample_names)
  rloop_score_vec <- clinical$Rloop_ssGSEA_Score[match(sample_base_vec, clinical$sample)]

  valid_idx <- !is.na(rloop_score_vec)
  expr_for_cor <- expr_sub[, valid_idx, drop = FALSE]
  rloop_for_cor <- rloop_score_vec[valid_idx]

  cor_res <- sapply(rownames(expr_for_cor), function(g) {
    cor(as.numeric(expr_for_cor[g, ]), rloop_for_cor, method = "spearman", use = "complete.obs")
  })

  cor_pval <- sapply(rownames(expr_for_cor), function(g) {
    tryCatch(
      cor.test(as.numeric(expr_for_cor[g, ]), rloop_for_cor, method = "spearman")$p.value,
      error = function(e) NA
    )
  })

  res <- data.frame(
    Gene = common_genes,
    MeanDiff = as.numeric(MeanDiff[common_genes]),
    Correlation = as.numeric(cor_res[common_genes]),
    Cor_Pvalue = as.numeric(cor_pval[common_genes]),
    Pass_MeanDiff = as.numeric(MeanDiff[common_genes]) > 0,
    Pass_Cor = as.numeric(cor_res[common_genes]) > 0.3,
    stringsAsFactors = FALSE
  )

  return(res)
}

res_list <- list()
dataset_names <- c("GSE81538", "METABRIC")

for (i in seq_along(expr_files)) {
  res <- calc_validation(expr_files[i], clin_files[i], assign_files[i], candidate_genes)
  colnames(res)[2:6] <- paste0(colnames(res)[2:6], "_", dataset_names[i])
  res_list[[dataset_names[i]]] <- res
}

merged_res <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), res_list)

meandiff_pass_cols <- grep("^Pass_MeanDiff_", colnames(merged_res), value = TRUE)
cor_pass_cols <- grep("^Pass_Cor_", colnames(merged_res), value = TRUE)

merged_res$Both_MeanDiff_pass <- rowSums(merged_res[, meandiff_pass_cols, drop = FALSE], na.rm = TRUE) == length(meandiff_pass_cols)
merged_res$Both_Cor_pass <- rowSums(merged_res[, cor_pass_cols, drop = FALSE], na.rm = TRUE) == length(cor_pass_cols)
merged_res$Final_Validated <- merged_res$Both_MeanDiff_pass & merged_res$Both_Cor_pass

write.csv(
  merged_res,
  file = file.path(outdir, "Candidate_genes_validation_results.csv"),
  row.names = FALSE
)

meandiff_cols <- grep("^MeanDiff_", colnames(merged_res), value = TRUE)
cor_cols <- grep("^Correlation_", colnames(merged_res), value = TRUE)

meandiff_mat <- as.matrix(merged_res[, meandiff_cols, drop = FALSE])
cor_mat <- as.matrix(merged_res[, cor_cols, drop = FALSE])

rownames(meandiff_mat) <- merged_res$Gene
rownames(cor_mat) <- merged_res$Gene

colnames(meandiff_mat) <- gsub("MeanDiff_", "", colnames(meandiff_mat))
colnames(cor_mat) <- gsub("Correlation_", "", colnames(cor_mat))

ord <- order(-merged_res$Final_Validated, -merged_res$Both_MeanDiff_pass, -rowMeans(meandiff_mat, na.rm = TRUE))
meandiff_mat <- meandiff_mat[ord, , drop = FALSE]
cor_mat <- cor_mat[ord, , drop = FALSE]

anno_row <- data.frame(
  Validated = ifelse(merged_res$Final_Validated[ord], "Pass", "Fail")
)
rownames(anno_row) <- rownames(meandiff_mat)

col_fun_md <- colorRamp2(
  c(min(meandiff_mat, na.rm = TRUE), 0, max(meandiff_mat, na.rm = TRUE)),
  c("#2166AC", "white", "#B2182B")
)

ht_md <- Heatmap(
  meandiff_mat,
  name = "MeanDiff",
  col = col_fun_md,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  column_title = "Mean difference\n(R-loop-high vs R-loop-low)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  left_annotation = rowAnnotation(
    Status = anno_row$Validated,
    col = list(Status = c("Pass" = "#1B9E77", "Fail" = "#D95F02")),
    show_legend = TRUE,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  ),
  heatmap_legend_param = list(
    title = "MeanDiff",
    direction = "horizontal",
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_width = unit(3.5, "cm")
  ),
  border = TRUE
)

col_fun_cor <- colorRamp2(
  c(min(cor_mat, na.rm = TRUE), 0, max(cor_mat, na.rm = TRUE)),
  c("#5E3C99", "white", "#E66101")
)

ht_cor <- Heatmap(
  cor_mat,
  name = "Spearman R",
  col = col_fun_cor,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_title = "Correlation with\nR-loop score",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Spearman R",
    direction = "horizontal",
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_width = unit(3.5, "cm")
  ),
  border = TRUE
)

pdf(file.path(outdir, "Validation_heatmap_MeanDiff_Correlation.pdf"), width = 4.8, height = 8.5)
draw(
  ht_md + ht_cor,
  heatmap_legend_side = "bottom",
  padding = unit(c(2, 2, 2, 10), "mm")
)
dev.off()

validation_summary <- data.frame(
  Category = factor(
    c("Total\ncandidates",
      "Pass MeanDiff\nin both cohorts",
      "Pass correlation\nin both cohorts",
      "Pass both\ncriteria"),
    levels = c("Total\ncandidates",
               "Pass MeanDiff\nin both cohorts",
               "Pass correlation\nin both cohorts",
               "Pass both\ncriteria")
  ),
  Count = c(
    nrow(merged_res),
    sum(merged_res$Both_MeanDiff_pass, na.rm = TRUE),
    sum(merged_res$Both_Cor_pass, na.rm = TRUE),
    sum(merged_res$Final_Validated, na.rm = TRUE)
  )
)

validation_summary$Percent <- round(validation_summary$Count / validation_summary$Count[1] * 100, 1)
validation_summary$Label <- paste0(validation_summary$Count, "\n(", validation_summary$Percent, "%)")

p_bar <- ggplot(validation_summary, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.6) +
  geom_text(aes(label = Label), vjust = -0.35, size = 4.5, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = c(
    "Total\ncandidates" = "#BDBDBD",
    "Pass MeanDiff\nin both cohorts" = "#66C2A5",
    "Pass correlation\nin both cohorts" = "#FC8D62",
    "Pass both\ncriteria" = "#8DA0CB"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Cross-cohort validation of candidate genes",
    x = "",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11.5, color = "black", face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  file.path(outdir, "Validation_summary_barplot.pdf"),
  p_bar,
  width = 8,
  height = 6
)

joint_df <- data.frame(
  MeanDiff_Status = c("Fail", "Fail", "Pass", "Pass"),
  Correlation_Status = c("Fail", "Pass", "Fail", "Pass"),
  Count = c(
    sum(!merged_res$Both_MeanDiff_pass & !merged_res$Both_Cor_pass, na.rm = TRUE),
    sum(!merged_res$Both_MeanDiff_pass &  merged_res$Both_Cor_pass, na.rm = TRUE),
    sum( merged_res$Both_MeanDiff_pass & !merged_res$Both_Cor_pass, na.rm = TRUE),
    sum( merged_res$Both_MeanDiff_pass &  merged_res$Both_Cor_pass, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

joint_df$MeanDiff_Status <- factor(joint_df$MeanDiff_Status, levels = c("Fail", "Pass"))
joint_df$Correlation_Status <- factor(joint_df$Correlation_Status, levels = c("Fail", "Pass"))
joint_df$Percent <- round(joint_df$Count / sum(joint_df$Count) * 100, 1)

joint_df$Is_Final <- (joint_df$MeanDiff_Status == "Pass" & joint_df$Correlation_Status == "Pass")

joint_df$Fill_Color <- ifelse(joint_df$Is_Final, "#FFD700",
                              ifelse(joint_df$Count == 0, "#F0F0F0", "#CCCCCC"))

p_square <- ggplot(joint_df, aes(x = Correlation_Status, y = MeanDiff_Status)) +
  geom_tile(aes(fill = I(Fill_Color)), color = "black", linewidth = 1.2, width = 0.88, height = 0.88) +
  geom_text(
    aes(label = paste0(Count, "\n(", Percent, "%)")),
    size = 5.5,
    fontface = "bold",
    color = ifelse(joint_df$Is_Final, "grey20", "grey40")
  ) +
  scale_x_discrete(
    labels = c("Fail" = "Fail", "Pass" = "Pass"),
    expand = c(0.15, 0.15)
  ) +
  scale_y_discrete(
    labels = c("Fail" = "Fail", "Pass" = "Pass"),
    expand = c(0.15, 0.15)
  ) +
  labs(
    title = "Joint validation criteria",
    x = "Positive correlation in both cohorts",
    y = "Positive mean difference\nin both cohorts"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 8)),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(
  file.path(outdir, "Validation_joint_square_tile.pdf"),
  p_square,
  width = 6,
  height = 6
)

validated_genes <- merged_res$Gene[merged_res$Final_Validated]

write.table(
  validated_genes,
  file = file.path(outdir, "Final_validated_Rloop_signature_genes.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

dataset_info <- list(
  Discovery = file.path(PROJECT_DIR, "results/bulk/discovery/04_C1C2_DEG_Analysis/GSE96058_expr_C1C2_LH.csv"),
  GSE81538  = file.path(PROJECT_DIR, "results/bulk/discovery/GSE81538/C1C2_Rloop_Violin/GSE81538_expr_C1C2_HL.csv"),
  METABRIC  = file.path(PROJECT_DIR, "results/bulk/discovery/METABRIC/C1C2_Rloop_Violin/METABRIC_expr_C1C2_HL.csv"),
  TCGA      = file.path(PROJECT_DIR, "results/bulk/TCGA/HRpos_HER2neg_tumor_only.txt"),
  GSE25066  = file.path(PROJECT_DIR, "results/bulk/GSE25066/GSE25066_expr_HighLow.csv")
)

dataset_order <- c("Discovery", "GSE81538", "METABRIC", "TCGA", "GSE25066")

gmt_file <- file.path(PROJECT_DIR, "results/bulk/GSE96058/01_Rloop_regulators.gmt")

validated_gene_file <- file.path(outdir, "Final_validated_Rloop_signature_genes.txt")

PLOT_COLORS <- c(
  "Discovery" = "#b43665",
  "GSE81538"  = "#4779bd",
  "METABRIC"  = "#187d79",
  "TCGA"      = "#d98c2b",
  "GSE25066"  = "#7a5ba8"
)

SSGSEA_ALPHA     <- 0.25
SSGSEA_NORMALIZE <- FALSE

POINT_SIZE    <- 1.2
POINT_ALPHA   <- 0.50
SMOOTH_ALPHA  <- 0.15
SMOOTH_LW     <- 0.9
SIDE_SCALE    <- 0.22
HIST_BINWIDTH <- 0.06
BASE_SIZE     <- 12

OUT_CSV_SUMMARY <- "Rescored_signature_vs_Rloop_correlation.csv"
OUT_PDF_SCATTER <- "Rescored_signature_vs_Rloop_correlation_scatter_5datasets.pdf"
PDF_WIDTH       <- 24
PDF_HEIGHT      <- 5.2

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggside)
  library(dplyr)
  library(gridExtra)
  library(data.table)
  library(GSEABase)
  library(GSVA)
  library(scales)
})

if (!exists("validated_genes")) {
  validated_genes <- fread(validated_gene_file, header = FALSE, data.table = FALSE)[, 1]
}
validated_genes <- unique(validated_genes)

geneSets <- getGmt(gmt_file, geneIdType = SymbolIdentifier())

gene_sets_list <- list()
for (i in seq_along(geneSets)) {
  gs <- geneSets[[i]]
  gene_sets_list[[setName(gs)]] <- geneIds(gs)
}

rloop_all_genes <- unique(unlist(gene_sets_list))

clean_sample_name4 <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^\"|\"$", "", x)
  x <- sub("_C1$",         "", x)
  x <- sub("_C2$",         "", x)
  x <- sub("_Rloop-high$", "", x)
  x <- sub("_Rloop-low$",  "", x)
  x <- sub("_tumor$",      "", x)
  return(x)
}

read_expr_for_scoring <- function(expr_file) {
  expr <- fread(expr_file, data.table = FALSE, check.names = FALSE)
  rownames(expr) <- make.unique(as.character(expr[, 1]))
  expr <- expr[, -1, drop = FALSE]
  expr_mat <- as.matrix(expr)
  mode(expr_mat) <- "numeric"
  colnames(expr_mat) <- clean_sample_name4(colnames(expr_mat))
  return(expr_mat)
}

minmax_norm <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  mn <- min(x, na.rm = TRUE)
  mx <- max(x, na.rm = TRUE)
  if (mx == mn) return(rep(0, length(x)))
  (x - mn) / (mx - mn)
}

format_pval_sci <- function(p, prefix = "P") {
  if (is.na(p) || is.null(p))        return(paste0(prefix, " = NA"))
  if (p == 0 || p < 1e-300)          return(paste0(prefix, " < 1e-300"))
  if (p < 1e-50)                     return(paste0(prefix, " < 1e-50"))
  if (p < 1e-10)                     return(paste0(prefix, " < 1e-10"))
  if (p < 0.0001)                    return(paste0(prefix, " < 0.0001"))
  if (p < 0.001)                     return(sprintf(paste0(prefix, " = %.2e"), p))
  if (p < 0.05)                      return(sprintf(paste0(prefix, " = %.3f"), p))
  sprintf(paste0(prefix, " = %.2f"), p)
}

analyze_dataset_rescore <- function(expr_file, dataset_name, validated_genes, rloop_all_genes) {

  expr_mat <- read_expr_for_scoring(expr_file)

  rloop_matched <- intersect(rloop_all_genes, rownames(expr_mat))
  sig_matched   <- intersect(validated_genes, rownames(expr_mat))

  if (length(rloop_matched) < 5)
    stop(paste0(dataset_name, " R-loopgene < 5, checkgene"))
  if (length(sig_matched) < 5)
    stop(paste0(dataset_name, " validated signaturegene < 5, checkgene"))

  score_gene_sets <- list(
    Rloop_Signature = rloop_matched,
    Final_Signature = sig_matched
  )

  ssgsea_param <- ssgseaParam(
    exprData  = expr_mat,
    geneSets  = score_gene_sets,
    alpha     = SSGSEA_ALPHA,
    normalize = SSGSEA_NORMALIZE
  )
  gsvaResult <- gsva(ssgsea_param, verbose = FALSE)

  df <- data.frame(
    Sample              = colnames(gsvaResult),
    Rloop_score_raw     = as.numeric(gsvaResult["Rloop_Signature", ]),
    Signature_score_raw = as.numeric(gsvaResult["Final_Signature",  ]),
    stringsAsFactors    = FALSE
  )

  df$Rloop_score     <- minmax_norm(df$Rloop_score_raw)
  df$Signature_score <- minmax_norm(df$Signature_score_raw)
  df <- df[complete.cases(df[, c("Rloop_score", "Signature_score")]), , drop = FALSE]

  if (nrow(df) < 10)
    stop(paste0(dataset_name, " analysissample < 10"))

  df$Rloop_z     <- as.numeric(scale(df$Rloop_score))
  df$Signature_z <- as.numeric(scale(df$Signature_score))

  cor_test <- cor.test(df$Signature_score, df$Rloop_score, method = "spearman")
  rho <- as.numeric(cor_test$estimate)
  pv  <- cor_test$p.value

  return(list(
    dataset          = dataset_name,
    expr_dim         = dim(expr_mat),
    rloop_genes_used = rloop_matched,
    sig_genes_used   = sig_matched,
    plot_df          = df,
    n                = nrow(df),
    cor              = rho,
    p                = pv
  ))
}

rescore_results <- list()
for (nm in dataset_order) {
  rescore_results[[nm]] <- analyze_dataset_rescore(
    expr_file    = dataset_info[[nm]],
    dataset_name = nm,
    validated_genes = validated_genes,
    rloop_all_genes = rloop_all_genes
  )
}

all_p   <- sapply(rescore_results, function(x) x$p)
all_fdr <- p.adjust(all_p, method = "BH")
for (nm in names(rescore_results)) {
  rescore_results[[nm]]$fdr <- all_fdr[nm]
}

rescore_summary <- data.frame(
  Dataset              = names(rescore_results),
  N                    = sapply(rescore_results, function(x) x$n),
  Rloop_Genes_Used     = sapply(rescore_results, function(x) length(x$rloop_genes_used)),
  Signature_Genes_Used = sapply(rescore_results, function(x) length(x$sig_genes_used)),
  Spearman_R           = round(sapply(rescore_results, function(x) x$cor), 4),
  Pvalue_raw           = sapply(rescore_results, function(x) x$p),
  FDR_raw              = sapply(rescore_results, function(x) x$fdr),
  Pvalue_label         = sapply(rescore_results, function(x) format_pval_sci(x$p,   "P")),
  FDR_label            = sapply(rescore_results, function(x) format_pval_sci(x$fdr, "FDR")),
  stringsAsFactors     = FALSE
)

write.csv(
  rescore_summary,
  file      = file.path(outdir, OUT_CSV_SUMMARY),
  row.names = FALSE,
  quote     = FALSE
)

for (nm in names(rescore_results)) {
  write.csv(
    rescore_results[[nm]]$plot_df,
    file      = file.path(outdir, paste0("Rescored_scores_", nm, ".csv")),
    row.names = FALSE,
    quote     = FALSE
  )
}

make_plot_ggside <- function(res_obj, color) {

  df  <- res_obj$plot_df
  nm  <- res_obj$dataset
  rho <- res_obj$cor
  pv  <- res_obj$p
  fdr <- res_obj$fdr
  n   <- res_obj$n

  x_range <- range(df$Signature_score, na.rm = TRUE)
  y_range <- range(df$Rloop_score,     na.rm = TRUE)
  x_pad   <- diff(x_range) * 0.10
  y_pad   <- diff(y_range) * 0.10
  x_lim   <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  y_lim   <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  x_bk    <- pretty(df$Signature_score, n = 4)
  y_bk    <- pretty(df$Rloop_score,     n = 4)

  stat_label <- paste0(
    "italic(rho) == ", sprintf("%.2f", rho), "*\",\"~",
    sub("P ", "italic(P) ", format_pval_sci(pv,  "P")),   "*\",\"~",
    sub("FDR ", "italic(FDR) ", format_pval_sci(fdr, "FDR")), "*\",\"~",
    "italic(n) == ", n
  )

  p <- ggplot(df, aes(x = Signature_score, y = Rloop_score)) +

    geom_point(
      color  = color,
      alpha  = POINT_ALPHA,
      size   = POINT_SIZE,
      shape  = 16
    ) +

    geom_smooth(
      method    = "lm",
      se        = TRUE,
      color     = color,
      fill      = color,
      alpha     = SMOOTH_ALPHA,
      linewidth = SMOOTH_LW,
      formula   = y ~ x
    ) +

    geom_xsidehistogram(
      aes(y = after_stat(density)),
      binwidth = HIST_BINWIDTH,
      fill     = color,
      alpha    = 0.35,
      color    = "white",
      linewidth = 0.2
    ) +
    geom_xsidedensity(
      aes(y = after_stat(density)),
      color     = color,
      linewidth = 0.9,
      alpha     = 0.9
    ) +
    scale_xsidey_continuous(labels = NULL, breaks = NULL) +

    geom_ysidehistogram(
      aes(x = after_stat(density)),
      binwidth = HIST_BINWIDTH,
      fill     = color,
      alpha    = 0.35,
      color    = "white",
      linewidth = 0.2
    ) +
    geom_ysidedensity(
      aes(x = after_stat(density)),
      color     = color,
      linewidth = 0.9,
      alpha     = 0.9
    ) +
    scale_ysidex_continuous(labels = NULL, breaks = NULL) +

    scale_x_continuous(limits = x_lim, breaks = x_bk) +
    scale_y_continuous(limits = y_lim, breaks = y_bk) +

    annotate(
      "text",
      x          = x_lim[2] - diff(x_lim) * 0.02,
      y          = y_lim[1] + diff(y_lim) * 0.02,
      label      = stat_label,
      hjust      = 1,
      vjust      = 0,
      size       = 3.2,
      color      = "grey20",
      parse      = TRUE,
      lineheight = 1.5
    ) +

    labs(
      title = nm,
      x     = "Signature ssGSEA score",
      y     = "R-loop ssGSEA score"
    ) +

    theme_classic(base_size = BASE_SIZE) +
    theme(

      plot.title        = element_text(hjust = 0.5, face = "bold",
                                       size = BASE_SIZE + 1,
                                       margin = margin(b = 4)),
      axis.title        = element_text(face = "bold", size = BASE_SIZE),
      axis.text         = element_text(color = "black", size = BASE_SIZE - 2),
      panel.border      = element_rect(color = "black", fill = NA,
                                       linewidth = 0.8),
      panel.grid.major  = element_line(color = "grey94", linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      axis.line         = element_blank(),
      plot.margin       = margin(8, 8, 8, 8),

      ggside.panel.scale   = SIDE_SCALE,
      ggside.axis.line     = element_blank(),
      ggside.axis.ticks    = element_blank(),
      ggside.panel.border  = element_rect(color = "grey80", linewidth = 0.4)
    )

  return(p)
}

plot_list <- list()
for (nm in dataset_order) {
  plot_list[[nm]] <- make_plot_ggside(rescore_results[[nm]], PLOT_COLORS[[nm]])
}

combined_plot <- gridExtra::arrangeGrob(
  grobs = plot_list,
  ncol  = 5,
  nrow  = 1
)

ggsave(
  filename = file.path(outdir, OUT_PDF_SCATTER),
  plot     = combined_plot,
  width    = PDF_WIDTH,
  height   = PDF_HEIGHT,
  dpi      = 300
)

gc()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

workdir <- file.path(PROJECT_DIR, "results/bulk/discovery")

expr_files <- c(
  GSE81538 = file.path(workdir, "GSE81538/C1C2_Rloop_Violin/GSE81538_expr_C1C2_HL.csv"),
  METABRIC = file.path(workdir, "METABRIC/C1C2_Rloop_Violin/METABRIC_expr_C1C2_HL.csv"),
  GSE96058 = file.path(workdir, "04_C1C2_DEG_Analysis/GSE96058_expr_C1C2_LH.csv")
)

gmt_file <- file.path(workdir, "h.all.v2026.1.Hs.symbols.gmt")

outdir <- file.path(workdir, "06_Hallmark_mult_ssGSEA_GSEA")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (nm in names(expr_files)) {
  if (!file.exists(expr_files[[nm]])) {
    stop("expression: ", expr_files[[nm]])
  }
}
if (!file.exists(gmt_file)) {
  stop(" GMT : ", gmt_file)
}

geneSets <- getGmt(gmt_file, geneIdType = SymbolIdentifier())
gene_sets_list <- lapply(geneSets, geneIds)
names(gene_sets_list) <- gsub("^HALLMARK_", "", sapply(geneSets, setName))

term2gene <- read.gmt(gmt_file)
colnames(term2gene) <- c("gs_name", "gene_symbol")
hallmark_names_full <- unique(term2gene$gs_name)

read_expr_matrix <- function(expr_file, dataset_name) {
  expr_df <- fread(expr_file, data.table = FALSE, check.names = FALSE)
  rownames(expr_df) <- expr_df[, 1]
  expr_df <- expr_df[, -1, drop = FALSE]

  expr_mat <- as.matrix(expr_df)
  mode(expr_mat) <- "numeric"

  if (any(duplicated(rownames(expr_mat)))) {
    expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  }

  colnames(expr_mat) <- paste0(dataset_name, "__", colnames(expr_mat))

  return(expr_mat)
}

run_ssgsea_one_cohort <- function(expr_mat, dataset_name, gene_sets_list, outdir) {

  ssgsea_param <- ssgseaParam(
    exprData  = expr_mat,
    geneSets  = gene_sets_list,
    alpha     = 0.25,
    normalize = FALSE
  )

  gsva_res <- gsva(ssgsea_param, verbose = FALSE)
  gsva_z   <- t(scale(t(gsva_res)))
  gsva_z[is.na(gsva_z)] <- 0

  write.csv(as.data.frame(gsva_res),
            file = file.path(outdir, paste0(dataset_name, "_ssGSEA_raw.csv")),
            quote = FALSE)

  write.csv(as.data.frame(gsva_z),
            file = file.path(outdir, paste0(dataset_name, "_ssGSEA_zscore.csv")),
            quote = FALSE)

  return(gsva_z)
}

run_gsea_one_cohort <- function(expr_mat, dataset_name, term2gene, outdir) {

  s <- colnames(expr_mat)
  group <- dplyr::case_when(
    grepl("_Rloop-high$", s) ~ "High",
    grepl("_Rloop-low$",  s) ~ "Low",
    TRUE ~ NA_character_
  )

  if (any(is.na(group))) {
    bad <- s[is.na(group)]
    stop(paste0(dataset_name, " samplegroup:\n",
                paste(head(bad, 10), collapse = "\n")))
  }

  group <- factor(group, levels = c("Low", "High"))

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit  <- lmFit(expr_mat, design)
  cont <- makeContrasts(High_vs_Low = High - Low, levels = design)
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)

  deg <- topTable(fit2,
                  coef = "High_vs_Low",
                  number = Inf,
                  adjust.method = "BH",
                  sort.by = "P")
  deg$Gene <- rownames(deg)

  write.csv(deg,
            file = file.path(outdir, paste0(dataset_name, "_limma_High_vs_Low.csv")),
            row.names = FALSE, quote = FALSE)

  gl <- setNames(deg$t, deg$Gene)
  gl <- gl[!is.na(names(gl)) & !duplicated(names(gl))]
  gl <- sort(gl, decreasing = TRUE)

  gsea_res <- GSEA(
    geneList      = gl,
    TERM2GENE     = term2gene,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    verbose       = FALSE,
    seed          = TRUE
  )

  df <- gsea_res@result
  df$Pathway <- gsub("^HALLMARK_", "", df$Description)

  write.csv(df,
            file = file.path(outdir, paste0(dataset_name, "_Hallmark_GSEA_results.csv")),
            row.names = FALSE, quote = FALSE)

  return(df)
}

expr_list <- list()
for (nm in names(expr_files)) {
  expr_list[[nm]] <- read_expr_matrix(expr_files[[nm]], nm)
}

ssgsea_list <- list()
for (nm in names(expr_list)) {
  ssgsea_list[[nm]] <- run_ssgsea_one_cohort(
    expr_mat = expr_list[[nm]],
    dataset_name = nm,
    gene_sets_list = gene_sets_list,
    outdir = outdir
  )
}

common_pathways <- Reduce(intersect, lapply(ssgsea_list, rownames))
if (length(common_pathways) == 0) {
  stop(",  ssGSEA.")
}

combined_mat <- do.call(cbind, lapply(ssgsea_list, function(x) {
  x[common_pathways, , drop = FALSE]
}))
combined_mat <- as.matrix(combined_mat)
combined_mat[is.na(combined_mat)] <- 0

write.csv(as.data.frame(combined_mat),
          file = file.path(outdir, "Hallmark_ssGSEA_combined_zscore.csv"),
          quote = FALSE)

anno_df <- data.frame(
  Sample = colnames(combined_mat),
  stringsAsFactors = FALSE
)

anno_df$Cohort <- NA_character_
anno_df$Cohort[grepl("^GSE81538__", anno_df$Sample)] <- "GSE81538"
anno_df$Cohort[grepl("^METABRIC__", anno_df$Sample)] <- "METABRIC"
anno_df$Cohort[grepl("^GSE96058__", anno_df$Sample)] <- "GSE96058"

anno_df$Group <- NA_character_
anno_df$Group[grepl("_Rloop-high$", anno_df$Sample)] <- "High"
anno_df$Group[grepl("_Rloop-low$",  anno_df$Sample)] <- "Low"

keep <- !is.na(anno_df$Cohort) & !is.na(anno_df$Group)
if (!all(keep)) {
  print(anno_df$Sample[!keep])
  anno_df <- anno_df[keep, , drop = FALSE]
  combined_mat <- combined_mat[, anno_df$Sample, drop = FALSE]
}

ord <- order(
  factor(anno_df$Cohort, levels = c("GSE81538", "METABRIC", "GSE96058")),
  factor(anno_df$Group, levels = c("Low", "High")),
  anno_df$Sample
)

anno_df <- anno_df[ord, , drop = FALSE]
combined_mat <- combined_mat[, anno_df$Sample, drop = FALSE]
anno_df <- anno_df[match(colnames(combined_mat), anno_df$Sample), , drop = FALSE]

if (nrow(anno_df) != ncol(combined_mat)) {
  stop("ssGSEA annotation.")
}

write.csv(anno_df,
          file = file.path(outdir, "Hallmark_ssGSEA_sample_annotation.csv"),
          row.names = FALSE, quote = FALSE)

idx_high_96058 <- which(anno_df$Cohort == "GSE96058" & anno_df$Group == "High")
if (length(idx_high_96058) > 0) {
  pathway_mean <- rowMeans(combined_mat[, idx_high_96058, drop = FALSE], na.rm = TRUE)
  combined_mat <- combined_mat[order(pathway_mean, decreasing = TRUE), , drop = FALSE]
}

ssgsea_pathway_order <- rownames(combined_mat)

group_col <- c("Low" = "#4DBBD5", "High" = "#E64B35")
cohort_col <- c("GSE81538" = "#00A087", "METABRIC" = "#3C5488", "GSE96058" = "#F39B7F")

top_anno <- HeatmapAnnotation(
  Cohort = anno_df$Cohort,
  Group  = anno_df$Group,
  col = list(
    Cohort = cohort_col,
    Group  = group_col
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size = unit(4, "mm"),
  gap = unit(1, "mm")
)

column_split_equal <- factor(
  anno_df$Cohort,
  levels = c("GSE81538", "METABRIC", "GSE96058")
)

col_fun_ssgsea <- colorRamp2(c(-2, 0, 2), c("#6fa6cf", "white", "#b43665"))

ht_ssgsea <- Heatmap(
  combined_mat,
  name = "ssGSEA z-score",
  col = col_fun_ssgsea,

  top_annotation = top_anno,
  column_split = column_split_equal,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,

  show_column_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.5),

  column_gap = unit(3, "mm"),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),

  border = TRUE,

  heatmap_legend_param = list(
    title = "ssGSEA z-score",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8.5),
    at = c(-2, -1, 0, 1, 2)
  )
)

pdf(file.path(outdir, "Hallmark_ssGSEA_multicohort_heatmap_equal_width.pdf"),
    width = 12, height = 13)

draw(
  ht_ssgsea,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(5, 5, 5, 5), "mm")
)

dev.off()

gsea_list <- list()
for (nm in names(expr_list)) {
  gsea_list[[nm]] <- run_gsea_one_cohort(
    expr_mat = expr_list[[nm]],
    dataset_name = nm,
    term2gene = term2gene,
    outdir = outdir
  )
}

all_pw <- gsub("^HALLMARK_", "", hallmark_names_full)

nes_mat <- matrix(NA_real_, nrow = length(all_pw), ncol = length(expr_files),
                  dimnames = list(all_pw, names(expr_files)))
fdr_mat <- matrix(NA_real_, nrow = length(all_pw), ncol = length(expr_files),
                  dimnames = list(all_pw, names(expr_files)))

for (nm in names(gsea_list)) {
  df <- gsea_list[[nm]]
  if (is.null(df) || nrow(df) == 0) next
  idx <- match(df$Pathway, rownames(nes_mat))
  ok  <- !is.na(idx)
  nes_mat[idx[ok], nm] <- df$NES[ok]
  fdr_mat[idx[ok], nm] <- df$p.adjust[ok]
}

write.csv(
  data.frame(Pathway = rownames(nes_mat), nes_mat, check.names = FALSE),
  file = file.path(outdir, "Hallmark_multicohort_NES_matrix.csv"),
  row.names = FALSE, quote = FALSE
)

write.csv(
  data.frame(Pathway = rownames(fdr_mat), fdr_mat, check.names = FALSE),
  file = file.path(outdir, "Hallmark_multicohort_FDR_matrix.csv"),
  row.names = FALSE, quote = FALSE
)

if (!exists("ssgsea_pathway_order")) {
  stop(" ssgsea_pathway_order,  GSEA  ssGSEA .")
}

nes_mat <- nes_mat[ssgsea_pathway_order, , drop = FALSE]
fdr_mat <- fdr_mat[ssgsea_pathway_order, , drop = FALSE]

sig_label <- matrix("", nrow = nrow(fdr_mat), ncol = ncol(fdr_mat),
                    dimnames = dimnames(fdr_mat))
sig_label[!is.na(fdr_mat) & fdr_mat < 0.05] <- "*"

lim <- max(abs(nes_mat), na.rm = TRUE)
if (!is.finite(lim) || lim == 0) lim <- 1

col_fun_gsea <- colorRamp2(
  c(-lim, 0, lim),
  c("#2166AC", "#FDDBC7", "#D6604D")
)

ht_gsea <- Heatmap(
  nes_mat,
  name = "NES",
  col  = col_fun_gsea,

  cluster_rows = FALSE,
  cluster_columns = FALSE,

  row_names_side = "left",
  row_names_gp   = gpar(fontsize = 8.5),

  column_names_gp       = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot      = 45,
  column_names_centered = FALSE,

  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_label[i, j] == "*") {
      grid.text("*", x, y,
                gp = gpar(fontsize = 12, col = "white", fontface = "bold"))
    }
  },

  border = TRUE,

  column_title = "Hallmark pathway enrichment (R-loop High vs Low)\nNES > 0: enriched in High group | NES < 0: enriched in Low group\nRow order synchronized with ssGSEA",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),

  heatmap_legend_param = list(
    title = "NES\n(High vs Low)",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  )
)

pdf(file.path(outdir, "Hallmark_GSEA_multicohort_heatmap.pdf"),
    width = 7, height = 14)

draw(
  ht_gsea,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(5, 5, 5, 5), "mm")
)

dev.off()

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

gsea_bubble_data <- as.data.frame(nes_mat) %>%
  tibble::rownames_to_column("Pathway") %>%
  tidyr::pivot_longer(
    cols = -Pathway,
    names_to = "Cohort",
    values_to = "NES"
  ) %>%
  left_join(
    as.data.frame(fdr_mat) %>%
      tibble::rownames_to_column("Pathway") %>%
      tidyr::pivot_longer(
        cols = -Pathway,
        names_to = "Cohort",
        values_to = "adj_p_val"
      ),
    by = c("Pathway", "Cohort")
  ) %>%
  mutate(

    p_category = case_when(
      is.na(adj_p_val) ~ NA_character_,
      adj_p_val >= 0.05 ~ NA_character_,
      adj_p_val < 0.001 ~ "< 0.001",
      adj_p_val < 0.01  ~ "< 0.01",
      adj_p_val < 0.05  ~ "< 0.05",
      TRUE ~ NA_character_
    ),
    p_category = factor(
      p_category,
      levels = c("< 0.001", "< 0.01", "< 0.05")
    )
  )

gsea_bubble_data <- gsea_bubble_data %>%
  complete(
    Pathway = ssgsea_pathway_order,
    Cohort = c("GSE81538", "METABRIC", "GSE96058")
  )

gsea_bubble_data$Pathway <- factor(
  gsea_bubble_data$Pathway,
  levels = rev(ssgsea_pathway_order)
)

gsea_bubble_data$Cohort <- factor(
  gsea_bubble_data$Cohort,
  levels = c("GSE81538", "METABRIC", "GSE96058")
)

p <- ggplot(gsea_bubble_data, aes(x = Cohort, y = Pathway)) +

  geom_hline(
    yintercept = seq_along(levels(gsea_bubble_data$Pathway)),
    linetype = 1,
    linewidth = 0.2,
    color = "grey80"
  ) +

  geom_point(
    aes(
      size = p_category,
      fill = NES
    ),
    alpha = 1,
    shape = 21,
    na.rm = TRUE
  ) +

  scale_y_discrete(
    expand = c(0, 1),
    drop = FALSE
  ) +

  scale_x_discrete(
    expand = c(0, 0.5),
    drop = FALSE
  ) +

  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES",
    na.value = NA
  ) +

  scale_size_manual(
    name = "adj. p-val.",
    breaks = c("< 0.001", "< 0.01", "< 0.05"),
    values = c(
      "< 0.001" = 6,
      "< 0.01"  = 4,
      "< 0.05"  = 2
    ),
    na.value = 0
  ) +

  labs(
    title = "Hallmark GSEA (R-loop High vs Low)",
    subtitle = "Only FDR < 0.05 are shown as bubbles",
    x = "",
    y = ""
  ) +

  theme_classic() +
  theme(
    axis.text.x = element_text(
      size = 11,
      angle = 45,
      hjust = 1,
      color = "black",
      face = "bold"
    ),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.ticks.length = unit(0.5, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(
      linetype = 2,
      linewidth = 0.8,
      color = "grey90"
    ),
    legend.frame = element_rect(color = "black"),
    legend.ticks = element_line(color = "black"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40")
  )

ggsave(
  filename = file.path(outdir, "Hallmark_GSEA_multicohort_bubble_heatmap.pdf"),
  plot = p,
  width = 6,
  height = 14
)

dev.off()

write.csv(
  data.frame(Pathway = ssgsea_pathway_order),
  file = file.path(outdir, "Hallmark_pathway_order_used_for_all_heatmaps.csv"),
  row.names = FALSE,
  quote = FALSE
)

for (nm in names(gsea_list)) {
  df <- gsea_list[[nm]]
  if (!is.null(df) && nrow(df) > 0) {
  }
}
