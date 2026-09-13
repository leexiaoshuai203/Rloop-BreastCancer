# GSE245601 malignant epithelial reclustering

PROJECT_DIR <- "."

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(SCP)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(base_dir, "08.tumorselect",
                        "GSE245601_Malignant_Epithelial.rds")
out_dir <- file.path(base_dir, "10.malignat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MAL_HVG <- 2000L
MAL_PC <- 20L
MAL_RES <- 0.08

if (!file.exists(input_file)) stop("not_foundinput: ", input_file)

mal_obj <- readRDS(input_file)
DefaultAssay(mal_obj) <- "RNA"

batch_col <- if ("sample_gsm" %in% colnames(mal_obj[[]])) {
  "sample_gsm"
} else {
  "orig.ident"
}
if (!batch_col %in% colnames(mal_obj[[]])) stop("missingsample")

print(table(mal_obj[[batch_col, drop = TRUE]]))

if (length(Layers(mal_obj[["RNA"]])) > 1) {
  mal_obj <- JoinLayers(mal_obj, assay = "RNA")
}

for (rr in Reductions(mal_obj)) mal_obj[[rr]] <- NULL
mal_obj@graphs <- list()
mal_obj@neighbors <- list()

old_cols <- grep(
  "(^seurat_clusters$)|(_snn_res\\.)|(^wsnn_res\\.)",
  colnames(mal_obj[[]]),
  value = TRUE
)
if (length(old_cols) > 0) {
  mal_obj@meta.data <- mal_obj@meta.data[
    , setdiff(colnames(mal_obj[[]]), old_cols), drop = FALSE
  ]
}

if (!"data" %in% Layers(mal_obj[["RNA"]])) {
  mal_obj <- NormalizeData(mal_obj, verbose = FALSE)
}

mal_obj <- FindVariableFeatures(
  mal_obj, nfeatures = MAL_HVG, verbose = FALSE
)
mal_obj <- ScaleData(
  mal_obj, features = VariableFeatures(mal_obj), verbose = FALSE
)
mal_obj <- RunPCA(
  mal_obj, features = VariableFeatures(mal_obj),
  npcs = MAL_PC, verbose = FALSE
)

mal_obj <- mal_obj %>%
  RunHarmony(
    group.by.vars = batch_col,
    verbose = FALSE
  )

mal_obj <- FindNeighbors(
  mal_obj, reduction = "harmony",
  dims = 1:MAL_PC, verbose = FALSE
)
mal_obj <- FindClusters(
  mal_obj, resolution = MAL_RES,
  algorithm = 1, random.seed = 1234,
  verbose = FALSE
)
mal_obj <- RunUMAP(
  mal_obj, reduction = "harmony",
  dims = 1:MAL_PC,
  n.neighbors = 30,
  min.dist = 0.30,
  seed.use = 1234,
  verbose = FALSE
)

cluster_ids <- sort(unique(as.integer(
  as.character(mal_obj$seurat_clusters)
)))
mal_obj$C_cluster <- factor(
  paste0("C", mal_obj$seurat_clusters),
  levels = paste0("C", cluster_ids)
)
Idents(mal_obj) <- "seurat_clusters"

print(table(mal_obj$C_cluster))

plot_df <- cbind(
  as.data.frame(Embeddings(mal_obj, "umap")),
  mal_obj[[]]
)
colnames(plot_df)[1:2] <- c("umap_1", "umap_2")
plot_df$Cluster <- mal_obj$C_cluster

cluster_cols <- setNames(
  SCP::palette_scp(n = length(cluster_ids)),
  levels(mal_obj$C_cluster)
)

label_pos <- plot_df %>%
  group_by(Cluster) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

p_cluster <- ggplot(plot_df, aes(umap_1, umap_2)) +
  geom_point(aes(color = Cluster), size = 0.01, alpha = 0.70) +
  stat_density_2d(
    aes(color = Cluster),
    geom = "density_2d",
    linewidth = 0.30,
    linetype = "dashed",
    contour_var = "ndensity",
    breaks = 0.15
  ) +
  geom_label_repel(
    data = label_pos,
    aes(label = Cluster, color = Cluster),
    fontface = "bold",
    size = 3.3,
    segment.color = "grey50",
    fill = scales::alpha("white", 0.75),
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = cluster_cols) +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = paste0(
      "Malignant epithelial cells (n = ",
      format(ncol(mal_obj), big.mark = ","), ")"
    )
  )

ggsave(
  file.path(out_dir, "01_malignant_epithelialUMAP.pdf"),
  p_cluster, width = 7, height = 6
)

p_sample <- SCP::CellDimPlot(
  srt = mal_obj,
  group.by = batch_col,
  reduction = "umap",
  pt.size = 0.01,
  alpha = 0.70,
  label = FALSE,
  theme_use = "theme_blank"
) +
  ggtitle("Malignant epithelial cells by sample") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    aspect.ratio = 1
  )

ggsave(
  file.path(out_dir, "02_malignant_epithelialsampleUMAP.pdf"),
  p_sample, width = 8, height = 6.5
)

meta_df <- data.frame(
  Cluster = mal_obj$C_cluster,
  Sample = mal_obj[[batch_col, drop = TRUE]]
)

cluster_count <- meta_df %>%
  count(Cluster, name = "CellNumber") %>%
  mutate(Percentage = 100 * CellNumber / sum(CellNumber))

cluster_sample <- meta_df %>%
  count(Cluster, Sample, name = "CellNumber") %>%
  group_by(Cluster) %>%
  mutate(Percentage = 100 * CellNumber / sum(CellNumber)) %>%
  ungroup()

write.csv(
  cluster_count,
  file.path(out_dir, "03_clustercell.csv"),
  row.names = FALSE
)
write.csv(
  cluster_sample,
  file.path(out_dir, "04_clustersample.csv"),
  row.names = FALSE
)

output_file <- file.path(
  out_dir,
  "GSE245601_Malignant_Epithelial_Reclustered.rds"
)
saveRDS(mal_obj, output_file)
