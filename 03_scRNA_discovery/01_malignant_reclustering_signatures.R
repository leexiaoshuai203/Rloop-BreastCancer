# GSE306201 malignant epithelial analysis

PROJECT_DIR <- "."

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

setwd(file.path(PROJECT_DIR, "results/GSE306201"))

STEP7_RDS <- "Step7_Epithelial_PhenoIter_SCIPACseed.rds"

OUT_DIR <- "Malignant"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MAL_HVG   <- 2000L
MAL_PC    <- 20L
MAL_RES   <- 0.08
BATCH_COL <- "orig.ident"

set.seed(1234)

epi_obj <- readRDS(STEP7_RDS)
DefaultAssay(epi_obj) <- "RNA"

if (!"phenoiter_tumor_label" %in% colnames(epi_obj@meta.data)) {
  stop("Step7 missing phenoiter_tumor_label")
}

print(table(epi_obj$phenoiter_tumor_label, useNA = "ifany"))

mal_obj <- subset(epi_obj, subset = phenoiter_tumor_label == "Tumor")
DefaultAssay(mal_obj) <- "RNA"

if (ncol(mal_obj) < 20) {
  stop("malignant_epithelialcell, ")
}

print(table(mal_obj[[BATCH_COL, drop = TRUE]]))

for (rr in Reductions(mal_obj))       mal_obj[[rr]] <- NULL
for (gg in names(mal_obj@graphs))     mal_obj@graphs[[gg]] <- NULL
for (nn in names(mal_obj@neighbors))  mal_obj@neighbors[[nn]] <- NULL

cluster_like_cols <- grep(
  "(^seurat_clusters$)|(^RNA_snn_res\\.)|(^SCT_snn_res\\.)|(^integrated_snn_res\\.)|(^wsnn_res\\.)",
  colnames(mal_obj@meta.data), value = TRUE
)
if (length(cluster_like_cols) > 0) mal_obj@meta.data[, cluster_like_cols] <- NULL

default_assay <- DefaultAssay(mal_obj)
layer_names   <- Layers(mal_obj[[default_assay]])

if (!"data" %in% layer_names) {
  mal_obj <- NormalizeData(mal_obj, verbose = FALSE)
} else {
}

mal_obj <- FindVariableFeatures(mal_obj, nfeatures = MAL_HVG, verbose = FALSE)
mal_obj <- ScaleData(mal_obj, verbose = FALSE)
mal_obj <- RunPCA(mal_obj, npcs = MAL_PC, verbose = FALSE)

mal_obj <- mal_obj %>%
  RunHarmony(
    group.by.vars = BATCH_COL,
    verbose       = FALSE
  )

mal_obj <- RunUMAP(mal_obj,    reduction = "harmony", dims = 1:MAL_PC, verbose = FALSE)
mal_obj <- FindNeighbors(mal_obj, reduction = "harmony", dims = 1:MAL_PC, verbose = FALSE)
mal_obj <- FindClusters(mal_obj,  resolution = MAL_RES, verbose = FALSE)

print(table(mal_obj$seurat_clusters))

plot_df <- cbind(
  as.data.frame(mal_obj@reductions$umap@cell.embeddings),
  mal_obj@meta.data
)
colnames(plot_df)[1:2] <- c("umap_1", "umap_2")

plot_df$cluster_label <- paste0("C", plot_df$seurat_clusters)
plot_df$cluster_label <- factor(plot_df$cluster_label, levels = sort(unique(plot_df$cluster_label)))

n_cluster <- length(levels(plot_df$cluster_label))

if (requireNamespace("SCP", quietly = TRUE)) {
  scp_cols <- SCP::palette_scp(n = n_cluster)
} else {
  scp_cols <- scales::hue_pal()(n_cluster)
}
cluster_cols <- setNames(scp_cols, levels(plot_df$cluster_label))
if ("C9" %in% names(cluster_cols)) cluster_cols["C9"] <- "#D95F02"

ct_pos <- plot_df %>%
  group_by(cluster_label) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), .groups = "drop")

xmin <- min(plot_df$umap_1); xmax <- max(plot_df$umap_1)
ymin <- min(plot_df$umap_2); ymax <- max(plot_df$umap_2)
xrange <- xmax - xmin;       yrange <- ymax - ymin

ax_len_x <- xrange * 0.22;   ax_len_y <- yrange * 0.22
ax_x0    <- xmin - xrange * 0.02
ax_y0    <- ymin - yrange * 0.02

p_mal_umap <- ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = cluster_label), size = 0.01, alpha = 0.6) +
  stat_density_2d(
    aes(color = cluster_label),
    geom = "density_2d", linewidth = 0.35, linetype = "dashed",
    contour_var = "ndensity", breaks = 0.15
  ) +
  geom_label_repel(
    data = ct_pos,
    aes(label = cluster_label, color = cluster_label),
    fontface = "bold", size = 3, box.padding = 0.35,
    segment.color = "grey45", fill = alpha("white", 0.75),
    show.legend = FALSE, max.overlaps = 20
  ) +
  scale_color_manual(
    values = cluster_cols,
    guide  = guide_legend(override.aes = list(size = 4, alpha = 1, linetype = 0), ncol = 2)
  ) +
  annotate("segment",
           x = ax_x0, xend = ax_x0 + ax_len_x, y = ax_y0, yend = ax_y0,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           linewidth = 0.7, color = "black"
  ) +
  annotate("segment",
           x = ax_x0, xend = ax_x0, y = ax_y0, yend = ax_y0 + ax_len_y,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           linewidth = 0.7, color = "black"
  ) +
  annotate("text",
           x = ax_x0 + ax_len_x / 2, y = ax_y0 - yrange * 0.04,
           label = "UMAP_1", fontface = "bold", size = 3.5
  ) +
  annotate("text",
           x = ax_x0 - xrange * 0.04, y = ax_y0 + ax_len_y / 2,
           label = "UMAP_2", angle = 90, fontface = "bold", size = 3.5
  ) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin      = margin(15, 15, 15, 20)
  ) +
  labs(
    color = "Cluster",
    title = paste0("Malignant Epithelial Cells - Harmony Re-clustered (n=", ncol(mal_obj), ")")
  )

ggsave(
  file.path(OUT_DIR, "Malignant_UMAP_contour_harmony_reclustered.pdf"),
  plot = p_mal_umap, width = 6.5, height = 5.8
)

saveRDS(mal_obj, file.path(OUT_DIR, "Malignant_Epithelial_Reclustered.rds"))

cluster_tab <- table(mal_obj$seurat_clusters)
print(cluster_tab)
for (cl in names(cluster_tab)) {
}

print(table(mal_obj$orig.ident[mal_obj$seurat_clusters == "0"]))

c8_markers <- FindMarkers(
  mal_obj,
  ident.1 = "0",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)
print(head(c8_markers, 10))

c_cells <- rownames(mal_obj@meta.data[mal_obj$seurat_clusters == "0", ])
other_cells <- rownames(mal_obj@meta.data[mal_obj$seurat_clusters != "0", ])
