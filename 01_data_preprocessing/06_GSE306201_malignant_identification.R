# GSE306201 malignant epithelial-cell identification

PROJECT_DIR <- "."

setwd(file.path(PROJECT_DIR, "data/raw/TCGA_BRCA"))
outdir <- file.path(PROJECT_DIR, "results/GSE306201/SCIPAC")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

pd <- read.delim("TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

pd2 <- pd[, c("sampleID", "ER_Status_nature2012")]
rownames(pd2) <- paste0(gsub("-", ".", pd2$sampleID), "A")

tpm <- read.table("tpm/tpm_BRCA_mRNA_symbolgenenames.txt",
                  sep = "\t", header = TRUE, row.names = 1,
                  stringsAsFactors = FALSE, check.names = FALSE)

er_pos_samples <- rownames(pd2)[pd2$ER_Status_nature2012 == "Positive"]
er_pos_tumor_01A <- er_pos_samples[grepl("01A$", er_pos_samples)]
er_pos_tumor_01A <- intersect(er_pos_tumor_01A, colnames(tpm))

normal_11A <- colnames(tpm)[grepl("11A$", colnames(tpm))]

final_samples <- c(er_pos_tumor_01A, normal_11A)
tpm_er_pos_all <- tpm[, final_samples, drop = FALSE]

out_file <- file.path(outdir, "ERpos_Tumor_Normal_TPM.csv")
write.csv(tpm_er_pos_all, file = out_file, quote = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SCIPAC)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(uwot)
  library(SCP)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

setwd(file.path(PROJECT_DIR, "results/GSE306201"))
ANNOTATED_RDS <- "GSE306201_Step4_the_end_Annotated.rds"
BULK_CSV      <- "SCIPAC/ERpos_Tumor_Normal_TPM.csv"
OUT_DIR       <- "SCIPAC/SCIPAC_allTCGA"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EPI_HVG   <- 2000L
EPI_PC    <- 20L
EPI_RES   <- 0.05

SCIPAC_HVG    <- 1000L
SCIPAC_PC     <- 30L
SCIPAC_RES    <- 1
ELA_ALPHA     <- 0.4
BT_SIZE       <- 50L
CI_ALPHA      <- 0.05
NFOLD         <- 10L
NUM_CORES     <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)

full_obj <- readRDS(ANNOTATED_RDS)
print(table(full_obj$manual_celltype))

epi_obj <- subset(full_obj, subset = manual_celltype == "Epithelial cells")
rm(full_obj); gc()

default_assay <- DefaultAssay(epi_obj)
layer_names <- Layers(epi_obj[[default_assay]])
if (!"data" %in% layer_names) {
  epi_obj <- NormalizeData(epi_obj, verbose = FALSE)} else {cat(" data layer,  NormalizeData()\n")}

epi_obj <- FindVariableFeatures(epi_obj, nfeatures = EPI_HVG, verbose = FALSE)
epi_obj <- ScaleData(epi_obj, verbose = FALSE)
epi_obj <- RunPCA(epi_obj, npcs = EPI_PC, verbose = FALSE)
epi_obj <- RunUMAP(epi_obj, dims = 1:EPI_PC, verbose = FALSE)
epi_obj <- FindNeighbors(epi_obj, dims = 1:EPI_PC, verbose = FALSE)
epi_obj <- FindClusters(epi_obj, resolution = EPI_RES, verbose = FALSE)

print(table(epi_obj$seurat_clusters))

plot_df <- cbind(
  as.data.frame(epi_obj@reductions$umap@cell.embeddings),
  epi_obj@meta.data
)
colnames(plot_df)[1:2] <- c("umap_1", "umap_2")
plot_df$cluster_label <- paste0("C", plot_df$seurat_clusters)

n_cluster <- length(unique(plot_df$cluster_label))
cluster_cols <- setNames(
  SCP::palette_scp(n = n_cluster),
  sort(unique(plot_df$cluster_label))
)

ct_pos <- plot_df %>%
  group_by(cluster_label) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

p_epi_umap <- ggplot(plot_df, aes(umap_1, umap_2)) +
  geom_point(aes(color = cluster_label), size = 0.01, alpha = 0.7) +
  geom_label_repel(
    data = ct_pos,
    aes(label = cluster_label, color = cluster_label),
    fontface = "bold", size = 3,
    box.padding = 0.4, segment.color = "grey50",
    fill = alpha("white", 0.7), show.legend = FALSE, max.overlaps = 20
  ) +
  scale_color_manual(values = cluster_cols,
                     guide = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  coord_fixed() +
  labs(color = "Cluster", title = "Epithelial Cells - Re-clustered") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    aspect.ratio = 1, plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(OUT_DIR, "Epi_01_Reclustered_UMAP.pdf"), p_epi_umap, width = 6.5, height = 6.5)

bulk_mat <- as.matrix(read.csv(BULK_CSV, row.names = 1, check.names = FALSE))
storage.mode(bulk_mat) <- "numeric"
if (anyDuplicated(rownames(bulk_mat))) bulk_mat <- rowsum(bulk_mat, rownames(bulk_mat))

bulk_group <- ifelse(grepl("01[A-Z]?$", colnames(bulk_mat)), "Tumor",
                     ifelse(grepl("11[A-Z]?$", colnames(bulk_mat)), "Normal", NA))
names(bulk_group) <- colnames(bulk_mat)
bulk_mat   <- bulk_mat[, !is.na(bulk_group), drop = FALSE]
bulk_group <- bulk_group[!is.na(bulk_group)]

print(table(bulk_group))

sc_mat <- tryCatch(
  as.matrix(GetAssayData(epi_obj, assay = "RNA", layer = "counts")),
  error = function(e) as.matrix(GetAssayData(epi_obj, assay = "RNA", slot = "counts"))
)
storage.mode(sc_mat) <- "numeric"
if (anyDuplicated(rownames(sc_mat))) sc_mat <- rowsum(sc_mat, rownames(sc_mat))

overlap_genes <- intersect(rownames(sc_mat), rownames(bulk_mat))
if (length(overlap_genes) < 500) warning("gene, checkgene")

sc_sub   <- sc_mat[overlap_genes, , drop = FALSE]
bulk_sub <- bulk_mat[overlap_genes, , drop = FALSE]
bulk_log <- log1p(bulk_sub)

seu_tmp <- CreateSeuratObject(sc_sub, project = "SCIPAC_EPI", min.cells = 3, min.features = 100)
seu_tmp <- NormalizeData(seu_tmp, verbose = FALSE)
seu_tmp <- FindVariableFeatures(seu_tmp, nfeatures = SCIPAC_HVG, verbose = FALSE)

sc_norm <- tryCatch(
  GetAssayData(seu_tmp, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(seu_tmp, assay = "RNA", slot = "data")
)
hvg_genes <- intersect(VariableFeatures(seu_tmp), rownames(bulk_log))

sc_prep   <- as.matrix(sc_norm[hvg_genes, , drop = FALSE])
bulk_prep <- as.matrix(bulk_log[hvg_genes, , drop = FALSE])
rm(seu_tmp); gc()

pca_res  <- SCIPAC::sc.bulk.pca(sc_prep, bulk_prep, do.pca.sc = FALSE, n.pc = SCIPAC_PC)
sc_rot   <- pca_res$sc.dat.rot
bulk_rot <- pca_res$bulk.dat.rot

ct_res <- SCIPAC::seurat.ct(sc_rot, res = SCIPAC_RES)
write.csv(ct_res$ct.assignment, file.path(OUT_DIR, "SCIPAC_Epi_cluster_assignment.csv"))

y_bin    <- factor(bulk_group, levels = c("Normal", "Tumor"))
bulk_use <- bulk_rot[names(y_bin), , drop = FALSE]

scipac_res <- SCIPAC::SCIPAC(
  bulk.dat      = bulk_use,
  y             = y_bin,
  family        = "binomial",
  ct.res        = ct_res,
  ela.net.alpha = ELA_ALPHA,
  bt.size       = BT_SIZE,
  numCores      = NUM_CORES,
  CI.alpha      = CI_ALPHA,
  nfold         = NFOLD
)

write.csv(scipac_res, file.path(OUT_DIR, "SCIPAC_Epi_Results.csv"))

common_cells <- intersect(colnames(epi_obj), rownames(scipac_res))

epi_obj@meta.data$scipac_lambda <- NA_real_
epi_obj@meta.data$scipac_sig <- NA_character_

epi_obj@meta.data[common_cells, "scipac_lambda"] <- scipac_res[common_cells, "Lambda.est"]
epi_obj@meta.data[common_cells, "scipac_sig"]    <- as.character(scipac_res[common_cells, "sig"])

epi_obj@meta.data$scipac_tumor_label <- dplyr::case_when(
  epi_obj@meta.data$scipac_sig == "Sig.pos" & epi_obj@meta.data$scipac_lambda > 0 ~ "High conf Tumor",
  epi_obj@meta.data$scipac_sig == "Sig.neg" & epi_obj@meta.data$scipac_lambda < 0 ~ "High conf Normal",
  TRUE ~ "Uncertain"
)

print(table(epi_obj@meta.data$scipac_sig, useNA = "ifany"))

print(table(epi_obj@meta.data$scipac_tumor_label, useNA = "ifany"))

umap_coords <- as.data.frame(Embeddings(epi_obj, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

plot_dat <- cbind(epi_obj@meta.data, umap_coords)
write.csv(plot_dat, file.path(OUT_DIR, "SCIPAC_Epi_PlotDat.csv"))

plot_dat$Lambda_clip <- pmax(pmin(plot_dat$scipac_lambda, 2), -2)

p_lambda <- ggplot(plot_dat, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = Lambda_clip), size = 0.01) +
  scale_color_gradient2(
    expression(Lambda), low = "#6fa6cf", mid = "grey80", high = "#b43665",
    midpoint = 0, limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2), labels = c("<=-2", "-1", "0", "1", ">=2")
  ) +
  coord_fixed() +
  labs(title = "SCIPAC Lambda - ER+ Tumor vs Normal") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "right",
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        aspect.ratio = 1)

ggsave(file.path(OUT_DIR, "SCIPAC_Epi_02_Lambda_UMAP.pdf"), p_lambda, width = 6.5, height = 6.5)

plot_dat$sig_label <- case_when(
  plot_dat$scipac_tumor_label == "High conf Tumor"  ~ "High conf Tumor (Sig.pos)",
  plot_dat$scipac_tumor_label == "High conf Normal" ~ "High conf Normal (Sig.neg)",
  TRUE ~ "Uncertain"
)
sig_cols <- c("High conf Tumor (Sig.pos)" = "#b43665", "High conf Normal (Sig.neg)" = "#6fa6cf", "Uncertain" = "grey85")

p_sig <- ggplot(plot_dat, aes(UMAP_1, UMAP_2)) +
  geom_point(data = subset(plot_dat, sig_label == "Uncertain"),
             color = "grey85", size = 0.01, alpha = 0.5) +
  geom_point(data = subset(plot_dat, sig_label != "Uncertain"),
             aes(color = sig_label), size = 0.2, alpha = 0.9) +
  scale_color_manual(values = sig_cols,
                     guide = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  coord_fixed() +
  labs(color = "SCIPAC Label", title = "SCIPAC Tumor/Normal Label") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "right",
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        aspect.ratio = 1)

ggsave(file.path(OUT_DIR, "SCIPAC_Epi_03_Sig_UMAP.pdf"), p_sig, width = 6.5, height = 6.5)

epi_obj@misc$step6_scipac <- list(
  bulk_group = bulk_group,
  scipac_res = scipac_res,
  ct_res     = ct_res,
  pca_res    = list(sc_rot = sc_rot, bulk_rot = bulk_rot),
  params     = list(SCIPAC_HVG = SCIPAC_HVG, SCIPAC_PC = SCIPAC_PC,
                    SCIPAC_RES = SCIPAC_RES, ELA_ALPHA = ELA_ALPHA,
                    BT_SIZE = BT_SIZE, EPI_RES = EPI_RES)
)

saveRDS(epi_obj, "Step6_Epithelial_SCIPAC.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

setwd(file.path(PROJECT_DIR, "results/GSE306201"))
STEP6_RDS <- "Step6_Epithelial_SCIPAC.rds"
OUT_DIR   <- "SCIPAC/MalignantClassify_SCIPACseed"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N_GENES   <- 50L
MAX_ITER      <- 20L
STABLE_RATIO  <- 0.95
LABEL_RATIO   <- 0.99
MIN_DEG_GENES <- 10L
set.seed(1234)

get_signature_from_scrna <- function(seu, mal_cells, nmal_cells, top_n, iter_label) {
  mal_cells  <- intersect(mal_cells, colnames(seu))
  nmal_cells <- intersect(nmal_cells, colnames(seu))

  if (length(mal_cells) < 10 || length(nmal_cells) < 10) {
    stop(paste0(" ", iter_label, " : cell < 10"))
  }

  seu$tmp_group <- NA_character_
  seu$tmp_group[colnames(seu) %in% mal_cells]  <- "Malignant"
  seu$tmp_group[colnames(seu) %in% nmal_cells] <- "NonMalignant"
  Idents(seu) <- "tmp_group"

  markers <- FindMarkers(
    seu, ident.1 = "Malignant", ident.2 = "NonMalignant",
    test.use = "wilcox", logfc.threshold = 0.1,
    min.pct = 0.1, only.pos = FALSE, verbose = FALSE
  )
  markers <- markers[markers$p_val_adj < 0.05, , drop = FALSE]

  if (nrow(markers) < MIN_DEG_GENES) {
    warning(paste0(" ", iter_label, " : gene < ", MIN_DEG_GENES))
  }

  tumor_sig <- rownames(markers[order(markers$avg_log2FC, decreasing = TRUE), , drop = FALSE])[1:top_n]
  nmal_m    <- markers[markers$avg_log2FC < 0, , drop = FALSE]
  normal_sig <- rownames(nmal_m[order(nmal_m$avg_log2FC, decreasing = FALSE), , drop = FALSE])[1:top_n]

  tumor_sig  <- tumor_sig[!is.na(tumor_sig)]
  normal_sig <- normal_sig[!is.na(normal_sig)]

  overlap <- intersect(tumor_sig, normal_sig)
  if (length(overlap) > 0) {
    tumor_sig  <- setdiff(tumor_sig, overlap)
    normal_sig <- setdiff(normal_sig, overlap)
  }

  list(tumor = tumor_sig, normal = normal_sig)
}

score_cells <- function(seu, sig_tumor, sig_normal, suffix) {
  sig_tumor  <- intersect(sig_tumor, rownames(seu))
  sig_normal <- intersect(sig_normal, rownames(seu))

  if (length(sig_tumor) < 5)  stop(paste0("Iter ", suffix, ": tumor sig gene <5"))
  if (length(sig_normal) < 5) stop(paste0("Iter ", suffix, ": normal sig gene <5"))

  seu <- AddModuleScore(seu, features = list(sig_tumor),  name = paste0("phenoiter_tumor_score_", suffix))
  seu <- AddModuleScore(seu, features = list(sig_normal), name = paste0("phenoiter_normal_score_", suffix))

  t_col <- paste0("phenoiter_tumor_score_", suffix, "1")
  n_col <- paste0("phenoiter_normal_score_", suffix, "1")
  seu[[paste0("phenoiter_diff_score_", suffix)]] <- as.numeric(scale(seu@meta.data[[t_col]] - seu@meta.data[[n_col]]))
  seu
}

kmeans_classify <- function(seu, suffix) {
  score_cols <- c(
    paste0("phenoiter_normal_score_", suffix, "1"),
    paste0("phenoiter_tumor_score_",  suffix, "1"),
    paste0("phenoiter_diff_score_",   suffix)
  )
  df <- seu@meta.data[, score_cols, drop = FALSE]
  colnames(df) <- c("normal_score", "tumor_score", "diff_score")

  set.seed(42)
  km <- kmeans(df[, "diff_score", drop = FALSE], centers = 2, nstart = 50)
  df$kmeans_class <- km$cluster

  c1_diff   <- mean(df$diff_score[df$kmeans_class == 1])
  c2_diff   <- mean(df$diff_score[df$kmeans_class == 2])
  malignant_cluster <- ifelse(c1_diff > c2_diff, 1, 2)
  df$classification <- ifelse(df$kmeans_class == malignant_cluster, "Malignant", "NonMalignant")
  df
}

calc_stability <- function(sig_new, sig_old) {
  mean(c(
    length(intersect(sig_new$tumor,  sig_old$tumor))  / max(length(sig_old$tumor), 1),
    length(intersect(sig_new$normal, sig_old$normal)) / max(length(sig_old$normal), 1)
  ))
}

calc_label_stability <- function(class_new, class_old) {
  common_cells <- intersect(names(class_new), names(class_old))
  if (length(common_cells) == 0) return(0)
  mean(class_new[common_cells] == class_old[common_cells])
}

plot_scatter <- function(df, iteration_label) {
  ggplot(df, aes(x = normal_score, y = tumor_score, color = classification)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_color_manual(values = c("Malignant" = "#D62728", "NonMalignant" = "#AEC7E8")) +
    labs(title = paste0("Iter ", iteration_label,
                        "  Malignant=", sum(df$classification == "Malignant"),
                        "  NonMalignant=", sum(df$classification == "NonMalignant")),
         x = "Normal Score", y = "Tumor Score", color = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          panel.border = element_rect(color = "black", fill = NA),
          aspect.ratio = 1, legend.position = "bottom")
}

epi_seu <- readRDS(STEP6_RDS)
DefaultAssay(epi_seu) <- "RNA"

if (!"scipac_tumor_label" %in% colnames(epi_seu@meta.data)) {
  stop("Step6 missing scipac_tumor_label")
}

seed_mal  <- colnames(epi_seu)[epi_seu$scipac_tumor_label == "High conf Tumor"]
seed_nmal <- colnames(epi_seu)[epi_seu$scipac_tumor_label == "High conf Normal"]

if (length(seed_mal) < 10 || length(seed_nmal) < 10) {
  stop("SCIPAC  Tumor/Normal cell, ")
}

sig_00 <- get_signature_from_scrna(epi_seu, seed_mal, seed_nmal, TOP_N_GENES, "00_SCIPAC")

write.csv(
  bind_rows(data.frame(iteration = "00_SCIPAC", type = "tumor",  gene = sig_00$tumor),
            data.frame(iteration = "00_SCIPAC", type = "normal", gene = sig_00$normal)),
  file.path(OUT_DIR, "Signature_iter00_SCIPACseed.csv"), row.names = FALSE
)

epi_seu     <- score_cells(epi_seu, sig_00$tumor, sig_00$normal, "01")
class_df_01 <- kmeans_classify(epi_seu, "01")
epi_seu$phenoiter_classification_01 <- class_df_01$classification

all_scatter_plots <- list(iter01 = plot_scatter(class_df_01, "01 (SCIPAC seed)"))
prev_sig          <- sig_00
prev_class        <- setNames(class_df_01$classification, rownames(class_df_01))
final_iter        <- "01"
stability_log     <- data.frame(iteration = character(), overlap = numeric(),
                                label_stability = numeric(),
                                n_malignant = integer(), n_nonmalignant = integer())

for (i in seq_len(MAX_ITER)) {
  iter_label <- sprintf("%02d", i + 1)

  mal_cells  <- names(prev_class)[prev_class == "Malignant"]
  nmal_cells <- names(prev_class)[prev_class == "NonMalignant"]

  if (length(mal_cells) < 10 || length(nmal_cells) < 10) {
break
  }

  new_sig <- tryCatch(
    get_signature_from_scrna(epi_seu, mal_cells, nmal_cells, TOP_N_GENES, iter_label),
    error = function(e) { cat("    FindMarkers :", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(new_sig)) break

  write.csv(
    bind_rows(data.frame(iteration = iter_label, type = "tumor",  gene = new_sig$tumor),
              data.frame(iteration = iter_label, type = "normal", gene = new_sig$normal)),
    file.path(OUT_DIR, paste0("Signature_iter", iter_label, "_scRNA.csv")),
    row.names = FALSE
  )

  epi_seu      <- score_cells(epi_seu, new_sig$tumor, new_sig$normal, iter_label)
  class_df_new <- kmeans_classify(epi_seu, iter_label)

  epi_seu[[paste0("phenoiter_classification_", iter_label)]] <- class_df_new$classification

  all_scatter_plots[[paste0("iter", iter_label)]] <-
    plot_scatter(class_df_new, paste0(iter_label, " (scRNA Wilcoxon)"))

  stability       <- calc_stability(new_sig, prev_sig)
  label_stability <- calc_label_stability(setNames(class_df_new$classification, rownames(class_df_new)), prev_class)
  n_mal           <- sum(class_df_new$classification == "Malignant")
  n_nmal          <- sum(class_df_new$classification == "NonMalignant")

  stability_log <- bind_rows(
    stability_log,
    data.frame(iteration = iter_label, overlap = round(stability, 4),
               label_stability = round(label_stability, 4),
               n_malignant = n_mal, n_nonmalignant = n_nmal)
  )

  prev_sig   <- new_sig
  prev_class <- setNames(class_df_new$classification, rownames(class_df_new))
  final_iter <- iter_label

  if (stability >= STABLE_RATIO && label_stability >= LABEL_RATIO) {
break
  }
}

epi_seu$phenoiter_tumor_label <- epi_seu@meta.data[[paste0("phenoiter_classification_", final_iter)]]
epi_seu$phenoiter_tumor_label <- ifelse(epi_seu$phenoiter_tumor_label == "Malignant", "Tumor", "Normal")

print(table(epi_seu$phenoiter_tumor_label, useNA = "ifany"))

umap_df <- as.data.frame(Embeddings(epi_seu, "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$phenoiter_label <- epi_seu$phenoiter_tumor_label

last_t_col <- paste0("phenoiter_tumor_score_", final_iter, "1")
last_n_col <- paste0("phenoiter_normal_score_", final_iter, "1")
last_d_col <- paste0("phenoiter_diff_score_", final_iter)

umap_df$tumor_score  <- epi_seu@meta.data[[last_t_col]]
umap_df$normal_score <- epi_seu@meta.data[[last_n_col]]
umap_df$diff_score   <- epi_seu@meta.data[[last_d_col]]

theme_umap <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(), aspect.ratio = 1)

p_class <- ggplot(umap_df[order(umap_df$phenoiter_label), ],
                  aes(UMAP_1, UMAP_2, color = phenoiter_label)) +
  geom_point(size = 0.01, alpha = 0.7) +
  scale_color_manual(values = c("Tumor" = "#D62728", "Normal" = "#AEC7E8"),
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = paste0("PhenoIter: Tumor vs Normal (iter ", final_iter, ")\n",
                      "Tumor=", sum(umap_df$phenoiter_label == "Tumor"),
                      "  Normal=", sum(umap_df$phenoiter_label == "Normal")),
       color = NULL) +
  theme_umap

p_diff <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = diff_score)) +
  geom_point(size = 0.01, alpha = 0.7) +
  scale_color_gradientn(
    colors = c("#6fa6cf", "#187d79", "#efb421", "#e58027", "#b43665"),
    name   = "Diff Score\n(z-scored)"
  ) +
  labs(title = "Differential Score (Tumor - Normal, z-scored)") +
  theme_umap +
  guides(color = guide_colorbar(barwidth = 0.8, barheight = 5))

p_scatter_final <- all_scatter_plots[[paste0("iter", final_iter)]] +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

if (nrow(stability_log) > 0) {
  p_stab <- ggplot(stability_log, aes(x = iteration, y = overlap, group = 1)) +
    geom_line(color = "#6fa6cf", linewidth = 0.8) +
    geom_point(color = "#b43665", size = 2) +
    geom_hline(yintercept = STABLE_RATIO, linetype = "dashed", color = "#8560af", linewidth = 0.6) +
    ylim(0, 1) +
    labs(title = "Signature Stability Convergence Curve",
         x = "Iteration",
         y = "Overlap Ratio") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA))
} else {
  p_stab <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Only 1 iteration") + theme_void()
}

ggsave(file.path(OUT_DIR, "PhenoIter_01_UMAP_Classification.pdf"), p_class, width = 6.5, height = 6.5)
ggsave(file.path(OUT_DIR, "PhenoIter_02_UMAP_DiffScore.pdf"), p_diff, width = 6.5, height = 6.5)
ggsave(file.path(OUT_DIR, "PhenoIter_03_Scatter_Final.pdf"), p_scatter_final, width = 6, height = 6)
ggsave(file.path(OUT_DIR, "PhenoIter_04_StabilityConvergence.pdf"), p_stab, width = 7, height = 4.5)
ggsave(file.path(OUT_DIR, "PhenoIter_00_Overview.pdf"),
       (p_class | p_diff) / (p_scatter_final | p_stab), width = 13, height = 12)

n_col_patch   <- min(3L, length(all_scatter_plots))
p_all_scatter <- wrap_plots(all_scatter_plots, ncol = n_col_patch)
ggsave(file.path(OUT_DIR, "PhenoIter_05_AllIter_Scatter.pdf"), p_all_scatter,
       width = n_col_patch * 5, height = ceiling(length(all_scatter_plots) / n_col_patch) * 5)

epi_seu@misc$step7_phenoiter <- list(
  init_method   = "SCIPAC_high_confidence_scRNA_DEG",
  final_iter    = final_iter,
  stability_log = stability_log,
  params        = list(TOP_N_GENES = TOP_N_GENES, MAX_ITER = MAX_ITER,
                       STABLE_RATIO = STABLE_RATIO, LABEL_RATIO = LABEL_RATIO)
)

saveRDS(epi_seu, "Step7_Epithelial_PhenoIter_SCIPACseed.rds")

write.csv(stability_log, file.path(OUT_DIR, "PhenoIter_StabilityLog.csv"), row.names = FALSE)
write.csv(
  data.frame(cell_barcode = colnames(epi_seu),
             scipac_tumor_label = epi_seu$scipac_tumor_label,
             phenoiter_tumor_label = epi_seu$phenoiter_tumor_label),
  file.path(OUT_DIR, "Step6_Step7_Labels.csv"), row.names = FALSE
)

suppressPackageStartupMessages({
  library(Seurat)
  library(fastCNV)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

setwd(file.path(PROJECT_DIR, "results/GSE306201"))
STEP7_RDS        <- "Step7_Epithelial_PhenoIter.rds"
FULL_RDS         <- "GSE306201_Step4_the_end_Annotated.rds"
OUT_DIR          <- "SCIPAC/FastCNV"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REFERENCE_LABELS <- c("T cells", "NK cells", "B cells")
TARGET_LABEL     <- "Epithelial cells"
MAX_REF_GLOBAL   <- 5000
SAMPLE_COL       <- "orig.ident"
set.seed(1234)

epi_seu  <- readRDS(STEP7_RDS)
full_obj <- readRDS(FULL_RDS)
DefaultAssay(epi_seu)  <- "RNA"
DefaultAssay(full_obj) <- "RNA"

if (packageVersion("Seurat") >= "5.0.0") {
  epi_seu  <- JoinLayers(epi_seu,  assay = "RNA")
  full_obj <- JoinLayers(full_obj, assay = "RNA")
}

print(table(epi_seu[[SAMPLE_COL]][, 1]))

print(table(full_obj$manual_celltype[full_obj$manual_celltype %in% REFERENCE_LABELS]))

ref_cells_sampled <- lapply(REFERENCE_LABELS, function(ct) {
  cells  <- colnames(full_obj)[full_obj$manual_celltype == ct]
  n_take <- min(length(cells), round(MAX_REF_GLOBAL / length(REFERENCE_LABELS)))
  if (length(cells) == 0) return(character(0))
  sample(cells, n_take)
}) %>% unlist()

common_genes <- intersect(rownames(epi_seu), rownames(full_obj))

ref_global_obj       <- full_obj[common_genes, ref_cells_sampled]
ref_global_obj$annot <- ref_global_obj$manual_celltype

ref_type_map <- setNames(ref_global_obj$manual_celltype,
                         colnames(ref_global_obj))

rm(full_obj); gc()

samples <- unique(epi_seu[[SAMPLE_COL]][, 1])

all_cnv_scores <- list()
all_ref_scores <- list()
all_ref_meta   <- list()

for (sid in samples) {

  epi_cells_orig <- colnames(epi_seu)[epi_seu[[SAMPLE_COL]][, 1] == sid]

  if (length(epi_cells_orig) < 10) {
next
  }

  epi_sub       <- epi_seu[common_genes, epi_cells_orig]
  epi_sub$annot <- TARGET_LABEL

  epi_sub_renamed <- RenameCells(epi_sub,
                                 new.names = paste0("EPI_", colnames(epi_sub)))
  ref_renamed     <- RenameCells(ref_global_obj,
                                 new.names = paste0("REF_", colnames(ref_global_obj)))

  cnv_input <- merge(epi_sub_renamed, ref_renamed)

  if (packageVersion("Seurat") >= "5.0.0") {
    cnv_input <- JoinLayers(cnv_input, assay = "RNA")
  }
  DefaultAssay(cnv_input) <- "RNA"

  counts_mat <- GetAssayData(cnv_input, assay = "RNA", layer = "counts")

  if (nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
    rm(epi_sub, epi_sub_renamed, ref_renamed, cnv_input, counts_mat); gc(); next
  }
  rm(counts_mat)

  sample_out_dir <- file.path(OUT_DIR, sid)
  dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(sample_out_dir)

  res <- tryCatch({
    fastCNV(
      seuratObj              = cnv_input,
      sampleName             = sid,
      referenceVar           = "annot",
      referenceLabel         = REFERENCE_LABELS,
      assay                  = "RNA",
      reClusterSeurat        = FALSE,
      getCNVPerChromosomeArm = FALSE,
      getCNVClusters         = FALSE,
      doPlot                 = TRUE,
      outputType             = "pdf"
    )
  }, error = function(e) {
NULL
  })

  setwd(old_wd)

  if (is.null(res)) {
    rm(epi_sub, epi_sub_renamed, ref_renamed, cnv_input); gc(); next
  }

  if (inherits(res, "Seurat")) {
    cnv_obj <- res
  } else if (is.list(res)) {
    seurat_idx <- which(sapply(res, inherits, "Seurat"))
    if (length(seurat_idx) == 0) {
print(names(res))
      rm(epi_sub, epi_sub_renamed, ref_renamed, cnv_input, res); gc(); next
    }
    cnv_obj <- res[[seurat_idx[1]]]
  }

  cnv_score_col <- grep("cnv_fraction|cnv.*score|score.*cnv",
                        colnames(cnv_obj@meta.data),
                        ignore.case = TRUE, value = TRUE)
  if (length(cnv_score_col) == 0) {
print(colnames(cnv_obj@meta.data))
    rm(epi_sub, epi_sub_renamed, ref_renamed, cnv_input, res, cnv_obj); gc(); next
  }
  cnv_score_col <- cnv_score_col[1]

  meta   <- cnv_obj@meta.data
  is_epi <- grepl("^EPI_", rownames(meta))
  is_ref <- grepl("^REF_", rownames(meta))

  meta_epi              <- meta[is_epi, , drop = FALSE]
  meta_epi$barcode_orig <- sub("^EPI_", "", rownames(meta_epi))
  score_epi             <- setNames(meta_epi[[cnv_score_col]],
                                    meta_epi$barcode_orig)

  n_match <- sum(meta_epi$barcode_orig %in% epi_cells_orig)
  if (n_match < length(epi_cells_orig) * 0.9)

  all_cnv_scores[[sid]] <- score_epi

  ref_score_vec         <- meta[is_ref, cnv_score_col]
  all_ref_scores[[sid]] <- ref_score_vec[is.finite(ref_score_vec)]

  meta_ref              <- meta[is_ref, , drop = FALSE]
  meta_ref$barcode_orig <- sub("^REF_", "", rownames(meta_ref))
  meta_ref$ref_type     <- ref_type_map[meta_ref$barcode_orig]

  all_ref_meta[[sid]] <- data.frame(
    cnv_score = meta_ref[[cnv_score_col]],
    ref_type  = meta_ref$ref_type,
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(cnv_score), !is.na(ref_type))

  thr99_sid <- quantile(ref_score_vec, 0.99, na.rm = TRUE)

  rm(epi_sub, epi_sub_renamed, ref_renamed, cnv_input, res, cnv_obj,
     meta, meta_epi, meta_ref); gc()
}

if (length(all_cnv_scores) == 0) stop("sample, check")

all_scores_vec <- do.call(c, all_cnv_scores)

ref_scores_pooled <- do.call(c, all_ref_scores)
ref_scores_pooled <- ref_scores_pooled[is.finite(ref_scores_pooled)]
print(summary(ref_scores_pooled))

thr95 <- quantile(ref_scores_pooled, 0.95)
thr99 <- quantile(ref_scores_pooled, 0.99)

thr_use <- thr99

scores_vec_fixed <- all_scores_vec
names(scores_vec_fixed) <- sub("^[^.]+\\.", "", names(all_scores_vec))

epi_seu$cnv_score       <- NA_real_
epi_seu$cnv_tumor_label <- NA_character_

matched <- intersect(names(scores_vec_fixed), colnames(epi_seu))

if (length(matched) < ncol(epi_seu) * 0.8)
  warning("!!  < 80%, check barcode ")

epi_seu$cnv_score[match(matched, colnames(epi_seu))] <-
  scores_vec_fixed[matched]

epi_seu$cnv_tumor_label <- ifelse(
  !is.na(epi_seu$cnv_score) & epi_seu$cnv_score > thr_use,
  "Tumor", "Normal"
)

print(table(epi_seu$cnv_tumor_label, useNA = "ifany"))

print(summary(epi_seu$cnv_score))

print(round(prop.table(
  table(epi_seu[[SAMPLE_COL]][, 1], epi_seu$cnv_tumor_label),
  margin = 1) * 100, 1))

if (!"umap" %in% names(epi_seu@reductions)) {
  stop("epi_seu  umap reduction,  RunUMAP")
}

p1 <- FeaturePlot(
  epi_seu,
  features   = "cnv_score",
  reduction  = "umap",
  pt.size    = 0.01,
  cols       = c("lightgrey", "#b43665")
) +
  labs(title = "CNV Score", color = "CNV Score") +
  theme_classic(base_size = 13)

pdf(file.path(OUT_DIR, "CNV_01_Score_UMAP.pdf"), width = 7, height = 6)
print(p1)
dev.off()

p2 <- DimPlot(
  epi_seu,
  group.by  = "cnv_tumor_label",
  reduction = "umap",
  pt.size   = 0.01,
  cols      = c("Tumor" = "#b43665", "Normal" = "#6fa6cf")
) +
  labs(title = "CNV-based Tumor / Normal") +
  theme_classic(base_size = 13)

pdf(file.path(OUT_DIR, "CNV_02_Label_UMAP.pdf"), width = 7, height = 6)
print(p2)
dev.off()

df_epi <- data.frame(
  cnv_score = epi_seu$cnv_score,
  group     = epi_seu$cnv_tumor_label,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(cnv_score), !is.na(group)) %>%
  filter(group %in% c("Tumor", "Normal"))

ref_scores_with_type <- bind_rows(all_ref_meta) %>%
  mutate(group = case_when(
    ref_type == "T cells"  ~ "T",
    ref_type == "NK cells" ~ "NK",
    ref_type == "B cells"  ~ "B",
    TRUE ~ ref_type
  )) %>%
  select(cnv_score, group) %>%
  filter(!is.na(cnv_score), !is.na(group))

df_vln <- bind_rows(df_epi, ref_scores_with_type) %>%
  mutate(group = factor(group, levels = c("Tumor", "Normal", "T", "NK", "B")))

print(table(df_vln$group, useNA = "ifany"))

p3 <- ggplot(df_vln, aes(x = group, y = cnv_score, fill = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
  geom_boxplot(
    width         = 0.12,
    outlier.shape = NA,
    fill          = "white",
    alpha         = 0.7
  ) +
  geom_hline(
    yintercept = thr_use,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.6
  ) +
  annotate(
    "text",
    x     = 5.3,
    y     = thr_use * 1.02,
    label = paste0("99% ref threshold\n(", round(thr_use, 3), ")"),
    hjust = 1,
    size  = 3.2,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "Tumor"  = "#b43665",
    "Normal" = "#6fa6cf",
    "T"      = "#187d79",
    "NK"     = "#efb421",
    "B"      = "#8560af"
  )) +
  labs(
    title = "CNV Score Distribution",
    x     = NULL,
    y     = "CNV Score"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

pdf(file.path(OUT_DIR, "CNV_03_Score_Violin.pdf"), width = 7, height = 6)
print(p3)
dev.off()

if (!"phenoiter_tumor_label" %in% colnames(epi_seu@meta.data)) {
  stop("epi_seu missing phenoiter_tumor_label,  Step7 results")
}

if (!"cnv_score" %in% colnames(epi_seu@meta.data)) {
  stop("epi_seu missing cnv_score, completed Step8 CNV score ")
}

df_epi_pheno <- data.frame(
  cnv_score = epi_seu$cnv_score,
  group     = epi_seu$phenoiter_tumor_label,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(cnv_score), !is.na(group)) %>%
  filter(group %in% c("Tumor", "Normal"))

ref_scores_with_type <- bind_rows(all_ref_meta) %>%
  mutate(group = case_when(
    ref_type == "T cells"  ~ "T",
    ref_type == "NK cells" ~ "NK",
    ref_type == "B cells"  ~ "B",
    TRUE ~ ref_type
  )) %>%
  select(cnv_score, group) %>%
  filter(!is.na(cnv_score), !is.na(group)) %>%
  filter(group %in% c("T", "NK", "B"))

df_vln_pheno <- bind_rows(df_epi_pheno, ref_scores_with_type) %>%
  mutate(group = factor(group, levels = c("Tumor", "Normal", "T", "NK", "B")))

print(table(df_vln_pheno$group, useNA = "ifany"))

print(
  df_vln_pheno %>%
    group_by(group) %>%
    summarise(
      n      = n(),
      min    = min(cnv_score, na.rm = TRUE),
      q25    = quantile(cnv_score, 0.25, na.rm = TRUE),
      median = median(cnv_score, na.rm = TRUE),
      mean   = mean(cnv_score, na.rm = TRUE),
      q75    = quantile(cnv_score, 0.75, na.rm = TRUE),
      max    = max(cnv_score, na.rm = TRUE)
    )
)

p4 <- ggplot(df_vln_pheno, aes(x = group, y = cnv_score, fill = group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = thr_use,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.6
  ) +
  annotate(
    "text",
    x     = 5.3,
    y     = thr_use * 1.02,
    label = paste0("99% ref threshold\n(", round(thr_use, 3), ")"),
    hjust = 1,
    size  = 3.2,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "Tumor"  = "#b43665",
    "Normal" = "#6fa6cf",
    "T"      = "#187d79",
    "NK"     = "#efb421",
    "B"      = "#8560af"
  )) +
  labs(
    title = "CNV Score Distribution by PhenoIter Label",
    subtitle = "Tumor/Normal labels from Step7 phenoiter_tumor_label; T/NK/B as reference cells",
    x = NULL,
    y = "CNV Score"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

pdf(file.path(OUT_DIR, "CNV_04_PhenoIterLabel_Violin.pdf"), width = 7.5, height = 6.2)
print(p4)
dev.off()

saveRDS(epi_seu, "Step8_Epithelial_FastCNV.rds")
