# WGCNA and R-loop-associated gene selection

PROJECT_DIR <- "."

gc()

workDir <- file.path(PROJECT_DIR, "results/bulk/discovery")
dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
setwd(workDir)

expr_file <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_exp.csv")
clin_file <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_HRp_HERn_clinical.csv")

project_name <- "GSE96058_Rloop_WGCNA"
out_dir <- "WGCNA_Rloop_Only"
dir.create(out_dir, showWarnings = FALSE)

subdirs <- c(
  "01_Quality_Control",
  "02_Network_Analysis",
  "03_Module_Detection",
  "04_Module_Trait_Analysis",
  "05_Gene_Analysis",
  "06_Data_Output"
)
for (d in subdirs) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}

rloop_score_col <- "Rloop_ssGSEA_Score"

min_sample_n <- 30
min_module_size <- 30
merge_cut_height <- 0.25
soft_threshold_range <- 1:20
network_type <- "signed"
cor_type <- "pearson"
var_quantile_threshold <- 0.45

suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
  library(ggplot2)
  library(gridExtra)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

theme_sci <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

expr_raw <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)
exprData <- as.matrix(expr_raw)
mode(exprData) <- "numeric"

if (max(exprData, na.rm = TRUE) > 100) {
  exprData <- log2(exprData + 1)
} else {
}

exprData <- normalizeBetweenArrays(exprData)

gene_var <- apply(exprData, 1, var, na.rm = TRUE)
cut_var <- quantile(gene_var, var_quantile_threshold, na.rm = TRUE)
exprData <- exprData[gene_var > cut_var, , drop = FALSE]

datExpr <- t(exprData)

clin <- read.csv(clin_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"sample" %in% colnames(clin)) {
  stop("Clinical file lacks column: sample")
}
if (!rloop_score_col %in% colnames(clin)) {
  stop(paste("Clinical file lacks column:", rloop_score_col))
}

clin <- clin[, c("sample", rloop_score_col)]
colnames(clin) <- c("sample", "Rloop_score")

clin$Rloop_score <- suppressWarnings(as.numeric(clin$Rloop_score))

common_samples <- intersect(rownames(datExpr), clin$sample)

if (length(common_samples) < min_sample_n) {
  stop("Too few common samples for WGCNA.")
}

datExpr <- datExpr[common_samples, , drop = FALSE]
clin <- clin[match(common_samples, clin$sample), , drop = FALSE]

keep <- !is.na(clin$Rloop_score)

datExpr <- datExpr[keep, , drop = FALSE]
clin <- clin[keep, , drop = FALSE]

if (nrow(datExpr) < min_sample_n) {
  stop("Too few samples after removing NA Rloop_score.")
}

datTraits <- data.frame(Rloop_score = clin$Rloop_score)
rownames(datTraits) <- clin$sample

write.csv(datTraits,
          file.path(out_dir, "06_Data_Output", "00_Rloop_trait.csv"),
          quote = FALSE)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  datTraits <- datTraits[gsg$goodSamples, , drop = FALSE]
}

sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file.path(out_dir, "01_Quality_Control", "01_SampleClustering.pdf"), width = 10, height = 6)
plot(sampleTree, main = "Sample clustering", xlab = "", sub = "", cex = 0.7)
dev.off()

tree_samples <- sampleTree$labels

if (is.null(tree_samples)) {
  stop("sampleTree$labels is NULL, cannot align trait colors.")
}

if (!all(tree_samples %in% rownames(datTraits))) {
  missing_samples <- tree_samples[!tree_samples %in% rownames(datTraits)]
  stop(paste("These tree samples are missing in datTraits:", paste(missing_samples, collapse = ", ")))
}

trait_vec <- datTraits[tree_samples, "Rloop_score"]

if (length(trait_vec) == 0) {
  stop("trait_vec is empty after alignment.")
}
if (all(is.na(trait_vec))) {
  stop("trait_vec is all NA after alignment.")
}
if (!is.numeric(trait_vec)) {
  stop("trait_vec is not numeric.")
}
if (!all(is.finite(trait_vec))) {
  stop("trait_vec contains non-finite values.")
}

traitColors <- numbers2colors(trait_vec, signed = TRUE)

traitColorMatrix <- matrix(traitColors, ncol = 1)
colnames(traitColorMatrix) <- "Rloop_score"
rownames(traitColorMatrix) <- tree_samples

pdf(file.path(out_dir, "01_Quality_Control", "02_SampleDendrogram_Trait.pdf"), width = 12, height = 7)
plotDendroAndColors(
  dendro = sampleTree,
  colors = traitColorMatrix,
  groupLabels = "Rloop_score",
  main = "Sample dendrogram and R-loop score"
)
dev.off()

sft <- pickSoftThreshold(
  datExpr,
  powerVector = soft_threshold_range,
  networkType = network_type,
  verbose = 5
)

fit_df <- data.frame(
  Power = sft$fitIndices[,1],
  SFT.R.sq = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  mean.k = sft$fitIndices[,5]
)

softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6

p1 <- ggplot(fit_df, aes(Power, SFT.R.sq)) +
  geom_point(color = "#b43665", size = 3) +
  geom_text(aes(label = Power), vjust = -0.8, size = 3) +
  geom_hline(yintercept = 0.9, linetype = 2, color = "#6fa6cf") +
  labs(title = "Scale independence", x = "Soft threshold (power)", y = "Signed R^2") +
  theme_sci()

p2 <- ggplot(fit_df, aes(Power, mean.k)) +
  geom_point(color = "#6fa6cf", size = 3) +
  geom_text(aes(label = Power), vjust = -0.8, size = 3) +
  labs(title = "Mean connectivity", x = "Soft threshold (power)", y = "Mean connectivity") +
  theme_sci()

pdf(file.path(out_dir, "02_Network_Analysis", "01_SoftThreshold.pdf"), width = 12, height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = network_type,
  networkType = network_type,
  minModuleSize = min_module_size,
  reassignThreshold = 0,
  mergeCutHeight = merge_cut_height,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- net$colors
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

pdf(file.path(out_dir, "03_Module_Detection", "01_GeneDendrogram_ModuleColors.pdf"), width = 12, height = 8)
plotDendroAndColors(
  geneTree,
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

module_table <- sort(table(moduleColors), decreasing = TRUE)
write.csv(data.frame(Module = names(module_table), Size = as.numeric(module_table)),
          file.path(out_dir, "06_Data_Output", "01_ModuleSize.csv"),
          row.names = FALSE)

print(module_table)

pdf(file.path(out_dir, "03_Module_Detection", "01_GeneDendrogram_ModuleColors.pdf"),
    width = 12, height = 8)

for (b in seq_along(net$dendrograms)) {
  plotDendroAndColors(
    net$dendrograms[[b]],
    moduleColors[net$blockGenes[[b]]],
    "Dynamic Tree Cut",
    dendroLabels = FALSE,
    hang         = 0.03,
    addGuide     = TRUE,
    guideHang    = 0.05,
    main         = "Gene dendrogram and module colors"
  )
}

dev.off()

dissTOM_all <- 1 - cor(datExpr)
geneTree_all <- hclust(as.dist(dissTOM_all), method = "average")

pdf(file.path(out_dir, "03_Module_Detection", "01_GeneDendrogram_ModuleColors.pdf"),
    width = 12, height = 8)
plotDendroAndColors(
  geneTree_all,
  moduleColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  main         = "Gene dendrogram and module colors"
)
dev.off()

MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs0 <- orderMEs(MEs0)

module_sizes <- table(moduleColors)

moduleTraitCor <- cor(MEs0, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

write.csv(moduleTraitCor,
          file.path(out_dir, "06_Data_Output", "01_Module_Trait_Correlations.csv"))
write.csv(moduleTraitPvalue,
          file.path(out_dir, "06_Data_Output", "02_Module_Trait_Pvalues.csv"))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file.path(out_dir, "04_Module_Trait_Analysis", "01_Module_Trait_ALL_Classic.pdf"),
    width = 8, height = max(6, nrow(moduleTraitCor) * 0.35 + 3))

par(mar = c(8, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.lab.x = 0.9,
               cex.lab.y = 0.8,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships (R-loop score only)")
dev.off()

core_cor <- moduleTraitCor[, "Rloop_score", drop = FALSE]
core_pval <- moduleTraitPvalue[, "Rloop_score", drop = FALSE]

core_text <- matrix(
  paste0(
    sprintf("%.2f", core_cor),
    "\n(",
    ifelse(core_pval < 1e-100, "<1e-100",
           ifelse(core_pval < 1e-10, sprintf("%.0e", core_pval),
                  ifelse(core_pval < 0.001, sprintf("%.4f", core_pval),
                         sprintf("%.3f", core_pval)))),
    ")"
  ),
  nrow = nrow(core_cor)
)

module_names <- rownames(core_cor)
module_colors_vec <- gsub("^ME", "", module_names)

color_annotation <- data.frame(
  Module = factor(module_colors_vec, levels = unique(module_colors_vec))
)
rownames(color_annotation) <- module_names

module_color_map <- setNames(unique(module_colors_vec), unique(module_colors_vec))

pdf(file.path(out_dir, "04_Module_Trait_Analysis", "02_Module_Rloop_Enhanced.pdf"),
    width = 8, height = max(8, nrow(core_cor) * 0.45 + 3))

pheatmap(core_cor,
         color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                                    "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582",
                                    "#D6604D", "#B2182B", "#67001F"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 10,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "Module-Rloop Relationship",
         border_color = "white",
         cellwidth = 120,
         cellheight = 28,
         display_numbers = core_text,
         number_color = "black",
         fontsize_number = 8,
         annotation_row = color_annotation,
         annotation_colors = list(Module = module_color_map),
         annotation_names_row = FALSE,
         legend = TRUE,
         legend_breaks = c(-1, -0.5, 0, 0.5, 1),
         legend_labels = c("-1", "-0.5", "0", "0.5", "1"))
dev.off()

sig_cor_threshold <- 0.30
sig_p_threshold <- 0.05

rloop_sig_modules <- rownames(moduleTraitCor)[
  abs(moduleTraitCor[, "Rloop_score"]) > sig_cor_threshold &
    moduleTraitPvalue[, "Rloop_score"] < sig_p_threshold
]

if (length(rloop_sig_modules) == 0) {

  rloop_sig_modules <- rownames(moduleTraitCor)[
    order(abs(moduleTraitCor[, "Rloop_score"]), decreasing = TRUE)
  ][1:min(2, nrow(moduleTraitCor))]
}

print(rloop_sig_modules)

candidate_info <- data.frame(
  Module = rloop_sig_modules,
  Rloop_cor = moduleTraitCor[rloop_sig_modules, "Rloop_score"],
  Rloop_p = moduleTraitPvalue[rloop_sig_modules, "Rloop_score"],
  GeneCount = as.numeric(module_sizes[gsub("^ME", "", rloop_sig_modules)]),
  stringsAsFactors = FALSE
)

candidate_info <- candidate_info[order(abs(candidate_info$Rloop_cor), decreasing = TRUE), ]

write.csv(candidate_info,
          file.path(out_dir, "06_Data_Output", "03_Candidate_Modules_Summary.csv"),
          row.names = FALSE)

print(candidate_info)

bestModule <- candidate_info$Module[1]
bestModuleColor <- sub("^ME", "", bestModule)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs0, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples = nrow(datExpr)))

geneTraitSignificance <- as.data.frame(cor(datExpr, datTraits$Rloop_score, use = "p"))
colnames(geneTraitSignificance) <- "GS.Rloop"
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nrow(datExpr)))
colnames(GSPvalue) <- "GS.Pvalue"

geneInfo <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  GS.Rloop = geneTraitSignificance[,1],
  GS.Pvalue = GSPvalue[,1],
  stringsAsFactors = FALSE
)

for (i in 1:ncol(geneModuleMembership)) {
  geneInfo[[paste0("MM_", colnames(geneModuleMembership)[i])]] <- geneModuleMembership[, i]
  geneInfo[[paste0("MMp_", colnames(geneModuleMembership)[i])]] <- MMPvalue[, i]
}

write.csv(geneInfo,
          file.path(out_dir, "06_Data_Output", "03_GeneInfo_AllModules.csv"),
          row.names = FALSE)

module_genes <- moduleColors == bestModuleColor
module_column <- paste0("ME", bestModuleColor)

df_plot <- data.frame(
  MM = abs(geneModuleMembership[module_genes, module_column]),
  GS = abs(geneTraitSignificance[module_genes, 1])
)

cor_test <- cor.test(df_plot$MM, df_plot$GS, method = "pearson")

p_mmgs <- ggplot(df_plot, aes(MM, GS)) +
  geom_point(color = bestModuleColor, alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  labs(
    title = paste0("MM vs GS in ", bestModuleColor, " module"),
    x = paste0("Module membership in ", bestModuleColor, " module"),
    y = "Gene significance for R-loop score"
  ) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(cor_test$estimate, 3),
                          "\nP = ", signif(cor_test$p.value, 3)),
           hjust = 1.1, vjust = 1.1, size = 4) +
  theme_sci()

pdf(file.path(out_dir, "05_Gene_Analysis", paste0("MM_vs_GS_", bestModuleColor, ".pdf")),
    width = 6, height = 6)
print(p_mmgs)
dev.off()

best_genes <- colnames(datExpr)[moduleColors == bestModuleColor]

write.table(best_genes,
            file.path(out_dir, "06_Data_Output", paste0("04_BestModule_", bestModuleColor, "_genes.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

for (mc in unique(moduleColors)) {
  tmp_genes <- colnames(datExpr)[moduleColors == mc]
  write.table(tmp_genes,
              file.path(out_dir, "06_Data_Output", paste0("Module_", mc, "_genes.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)
})

workDir <- file.path(PROJECT_DIR, "results/bulk/discovery")
dir.create(workDir, recursive = TRUE, showWarnings = FALSE)
setwd(workDir)

expr_file <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_exp.csv")
clin_file <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_HRp_HERn_clinical.csv")
gene_file <- file.path(PROJECT_DIR, "results/bulk/discovery/WGCNA_Rloop_Only/06_Data_Output/Module_blue_genes.txt")

out_dir <- file.path(workDir, "BlueModule_Rloop_Filter")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rloop_score_col <- "Rloop_ssGSEA_Score"

cor_method <- "spearman"
cor_cutoff <- 0.40
fdr_cutoff <- 0.05

top_n_heatmap <- 50

expr_raw <- read.csv(expr_file, header = TRUE, row.names = 1, check.names = FALSE)
expr_mat <- as.matrix(expr_raw)
mode(expr_mat) <- "numeric"

if (max(expr_mat, na.rm = TRUE) > 100) {
  expr_mat <- log2(expr_mat + 1)
} else {
}

clin <- read.csv(clin_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

if (!"sample" %in% colnames(clin)) {
  stop("Clinical file lacks column: sample")
}
if (!rloop_score_col %in% colnames(clin)) {
  stop(paste("Clinical file lacks column:", rloop_score_col))
}

clin <- clin[, c("sample", rloop_score_col)]
colnames(clin) <- c("sample", "Rloop_score")
clin$Rloop_score <- suppressWarnings(as.numeric(clin$Rloop_score))

clin <- clin[!is.na(clin$Rloop_score), , drop = FALSE]

common_samples <- intersect(colnames(expr_mat), clin$sample)

if (length(common_samples) < 30) {
  stop("Too few common samples after matching expression and clinical data.")
}

expr_mat <- expr_mat[, common_samples, drop = FALSE]
clin <- clin[match(common_samples, clin$sample), , drop = FALSE]

blue_genes <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE,
                         sep = "\t", quote = "", fill = TRUE)[, 1]
blue_genes <- unique(trimws(as.character(blue_genes)))
blue_genes <- blue_genes[!is.na(blue_genes) & blue_genes != ""]

matched_genes <- intersect(blue_genes, rownames(expr_mat))
missing_genes <- setdiff(blue_genes, rownames(expr_mat))

write.table(data.frame(Missing_Genes = missing_genes),
            file = file.path(out_dir, "01_Missing_blue_module_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

blue_expr <- expr_mat[matched_genes, , drop = FALSE]

rloop_vec <- clin$Rloop_score
names(rloop_vec) <- clin$sample

cor_res <- lapply(rownames(blue_expr), function(g) {
  x <- as.numeric(blue_expr[g, ])
  ct <- suppressWarnings(cor.test(x, rloop_vec, method = cor_method, exact = FALSE))
  data.frame(
    GeneSymbol = g,
    Correlation = unname(ct$estimate),
    Pvalue = ct$p.value,
    stringsAsFactors = FALSE
  )
})

cor_df <- do.call(rbind, cor_res)
cor_df$FDR <- p.adjust(cor_df$Pvalue, method = "BH")
cor_df$AbsCorrelation <- abs(cor_df$Correlation)

selected_df <- subset(cor_df, AbsCorrelation >= cor_cutoff & FDR < fdr_cutoff)
selected_df <- selected_df[order(selected_df$AbsCorrelation, decreasing = TRUE), ]

write.csv(cor_df,
          file = file.path(out_dir, "02_BlueModule_RloopCorrelation_AllGenes.csv"),
          row.names = FALSE)

write.table(selected_df$GeneSymbol,
            file = file.path(out_dir, "03_BlueModule_RloopCorrelation_SelectedGenes.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

selected_expr <- blue_expr[selected_df$GeneSymbol, , drop = FALSE]

write.csv(
  data.frame(GeneSymbol = rownames(selected_expr), selected_expr, check.names = FALSE),
  file = file.path(out_dir, "04_BlueModule_RloopCorrelation_SelectedExpression.csv"),
  row.names = FALSE
)

plot_df <- cor_df
plot_df$Significant <- ifelse(plot_df$AbsCorrelation >= cor_cutoff & plot_df$FDR < fdr_cutoff,
                              "Selected", "Not selected")
plot_df$log10FDR <- -log10(plot_df$FDR + 1e-300)

p_cor <- ggplot(plot_df, aes(x = Correlation, y = log10FDR)) +
  geom_point(aes(color = Significant), alpha = 0.75, size = 2) +
  scale_color_manual(values = c("Not selected" = "grey70", "Selected" = "#b43665")) +
  geom_vline(xintercept = c(-cor_cutoff, cor_cutoff), linetype = "dashed", color = "#6fa6cf") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "#6fa6cf") +
  labs(
    title = "Blue Module Genes Correlated with R-loop Score",
    subtitle = paste0("Method: ", cor_method,
                      " | Selected genes = ", nrow(selected_df)),
    x = "Correlation with R-loop score",
    y = expression(-log[10](FDR))
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(color = "black")
  )

ggsave(file.path(out_dir, "05_BlueModule_RloopCorrelation_Volcano.pdf"),
       plot = p_cor, width = 8, height = 6)
ggsave(file.path(out_dir, "05_BlueModule_RloopCorrelation_Volcano.png"),
       plot = p_cor, width = 8, height = 6, dpi = 600)

if (nrow(selected_df) >= 2) {
  top_genes <- head(selected_df$GeneSymbol, min(top_n_heatmap, nrow(selected_df)))
  heat_expr <- blue_expr[top_genes, , drop = FALSE]

  heat_expr_scaled <- t(scale(t(heat_expr)))

  anno_col <- data.frame(Rloop_score = clin$Rloop_score)
  rownames(anno_col) <- clin$sample

  pdf(file.path(out_dir, "06_BlueModule_RloopCorrelation_TopGenes_Heatmap.pdf"),
      width = 10, height = max(6, 0.18 * nrow(heat_expr_scaled) + 3))
  pheatmap(
    heat_expr_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "none",
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = anno_col,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Top R-loop-correlated genes in blue module",
    fontsize_row = 8,
    fontsize_col = 8,
    border_color = NA
  )
  dev.off()
}

summary_df <- data.frame(
  Input_Blue_Genes = length(blue_genes),
  Matched_Blue_Genes = length(matched_genes),
  Missing_Blue_Genes = length(missing_genes),
  Selected_Genes = nrow(selected_df),
  Correlation_Method = cor_method,
  Correlation_Cutoff = cor_cutoff,
  FDR_Cutoff = fdr_cutoff,
  stringsAsFactors = FALSE
)

write.csv(summary_df,
          file = file.path(out_dir, "07_BlueModule_RloopCorrelation_Summary.csv"),
          row.names = FALSE)

gc()

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(statmod)
  library(ggplot2)
  library(ggrepel)
})

workDir <- file.path(PROJECT_DIR, "results/bulk/discovery")
dir.create(workDir, showWarnings = FALSE, recursive = TRUE)
setwd(workDir)

expr_file <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_expr_HighLow.csv")

gene_file <- file.path(PROJECT_DIR, "results/bulk/discovery/BlueModule_Rloop_Filter/03_BlueModule_RloopCorrelation_SelectedGenes.txt")

fitThr <- 1.0

minMeanForFit <- 1

label_top_n <- 15

col_point <- "#e58027"
col_fit   <- "#187d79"
col_thr   <- "#6fa6cf"
col_hvg   <- "#b43665"

getMostVarGenes <- function(data, fitThr = 1.0, minMeanForFit = 1, label_top_n = 15) {

  data_no0 <- as.matrix(data[rowSums(data, na.rm = TRUE) > 0, , drop = FALSE])

  if (nrow(data_no0) == 0) {
    stop("gene0, HVG.")
  }

  meanGeneExp <- rowMeans(data_no0, na.rm = TRUE)
  varGenes    <- rowVars(data_no0, na.rm = TRUE)

  cv2 <- varGenes / (meanGeneExp^2)

  useForFit <- meanGeneExp >= minMeanForFit

  if (sum(useForFit) < 3) {
    stop("gene3, ,  minMeanForFit.")
  }

  fit <- glmgam.fit(
    cbind(
      a0      = 1,
      a1tilde = 1 / meanGeneExp[useForFit]
    ),
    cv2[useForFit]
  )

  fitModel <- fit$fitted.values
  names(fitModel) <- names(meanGeneExp[useForFit])

  obs_cv2_use <- cv2[useForFit]
  hvg_flag <- obs_cv2_use > fitModel * fitThr
  HVGenes <- names(obs_cv2_use)[hvg_flag]

  stat_df <- data.frame(
    GeneSymbol = names(meanGeneExp),
    MeanExp = as.numeric(meanGeneExp),
    Variance = as.numeric(varGenes),
    CV2 = as.numeric(cv2),
    UsedForFit = names(meanGeneExp) %in% names(fitModel),
    FittedCV2 = NA_real_,
    ExcessRatio = NA_real_,
    IsHVG = FALSE,
    stringsAsFactors = FALSE
  )

  stat_df$FittedCV2[match(names(fitModel), stat_df$GeneSymbol)] <- fitModel
  stat_df$ExcessRatio <- stat_df$CV2 / stat_df$FittedCV2
  stat_df$IsHVG[stat_df$GeneSymbol %in% HVGenes] <- TRUE

  stat_df <- stat_df[order(stat_df$IsHVG, stat_df$ExcessRatio, decreasing = TRUE), ]

  plot_df <- stat_df[stat_df$UsedForFit, ]
  plot_df$log10Mean <- log10(plot_df$MeanExp)
  plot_df$log10CV2  <- log10(plot_df$CV2)
  plot_df$log10Fit  <- log10(plot_df$FittedCV2)
  plot_df$log10Thr  <- log10(plot_df$FittedCV2 * fitThr)

  label_df <- plot_df[plot_df$IsHVG, ]
  label_df <- label_df[order(label_df$ExcessRatio, decreasing = TRUE), ]
  if (nrow(label_df) > label_top_n) {
    label_df <- label_df[1:label_top_n, ]
  }

  line_df <- plot_df[order(plot_df$log10Mean), ]

  p <- ggplot(plot_df, aes(x = log10Mean, y = log10CV2)) +
    geom_point(color = col_point, size = 1.2, alpha = 0.65) +
    geom_point(
      data = subset(plot_df, IsHVG),
      color = col_hvg, size = 1.5, alpha = 0.9
    ) +
    geom_line(
      data = line_df,
      aes(y = log10Fit),
      color = col_fit, linewidth = 1.1
    ) +
    geom_line(
      data = line_df,
      aes(y = log10Thr),
      color = col_thr, linewidth = 1.1, linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = GeneSymbol),
      size = 3.4,
      color = "black",
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey50",
      max.overlaps = 100
    ) +
    labs(
      title = "Highly Variable Gene Selection in WGCNA Blue Module",
      subtitle = paste0(
        "Threshold = fitted CV² x ", fitThr,
        " | Mean expression cutoff = ", minMeanForFit,
        " | HVGs = ", sum(stat_df$IsHVG, na.rm = TRUE)
      ),
      x = "log10(Mean expression)",
      y = expression(log[10](CV^2))
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11.5, hjust = 0.5, color = "black"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8)
    )

  HVG_mat <- data_no0[rownames(data_no0) %in% HVGenes, , drop = FALSE]

  return(list(
    HVG_mat = HVG_mat,
    stat_df = stat_df,
    plot_df = plot_df,
    plot = p
  ))
}

expr_raw <- read.csv(
  expr_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

expr_mat <- as.matrix(expr_raw)
mode(expr_mat) <- "numeric"

blue_genes <- read.table(
  gene_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  sep = "\t",
  quote = "",
  fill = TRUE
)[, 1]

blue_genes <- unique(trimws(as.character(blue_genes)))
blue_genes <- blue_genes[!is.na(blue_genes) & blue_genes != ""]

matched_genes <- intersect(blue_genes, rownames(expr_mat))
missing_genes <- setdiff(blue_genes, rownames(expr_mat))

write.table(
  data.frame(Missing_Genes = missing_genes),
  file = "01_Missing_blue_module_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

blue_expr <- expr_mat[matched_genes, , drop = FALSE]

write.csv(
  data.frame(GeneSymbol = rownames(blue_expr), blue_expr, check.names = FALSE),
  file = "02_Blue_module_expression_matrix.csv",
  row.names = FALSE
)

hvg_res <- getMostVarGenes(
  data = blue_expr,
  fitThr = fitThr,
  minMeanForFit = minMeanForFit,
  label_top_n = label_top_n
)

hvg_mat <- hvg_res$HVG_mat
stat_df <- hvg_res$stat_df
p_hvg   <- hvg_res$plot

write.table(
  rownames(hvg_mat),
  file = "03_HVG_blue_module_gene_list.txt",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.csv(
  data.frame(GeneSymbol = rownames(hvg_mat), hvg_mat, check.names = FALSE),
  file = "04_HVG_blue_module_expression_matrix.csv",
  row.names = FALSE
)

write.csv(
  stat_df,
  file = "05_HVG_statistics_all_blue_module_genes.csv",
  row.names = FALSE
)

write.csv(
  subset(stat_df, IsHVG),
  file = "06_HVG_statistics_selected_genes.csv",
  row.names = FALSE
)

ggsave(
  filename = "07_HVG_selection_plot_blue_module.pdf",
  plot = p_hvg,
  width = 8.2,
  height = 6.8
)

ggsave(
  filename = "07_HVG_selection_plot_blue_module.png",
  plot = p_hvg,
  width = 8.2,
  height = 6.8,
  dpi = 600
)

summary_df <- data.frame(
  Input_Blue_Genes = length(blue_genes),
  Matched_Genes = length(matched_genes),
  Missing_Genes = length(missing_genes),
  HVG_Count = nrow(hvg_mat),
  fitThr = fitThr,
  minMeanForFit = minMeanForFit,
  stringsAsFactors = FALSE
)

write.csv(summary_df, "08_HVG_summary_report.csv", row.names = FALSE)
