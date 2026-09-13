# CytoTRACE2 and Slingshot analysis

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(scop)
})

# Paths
project_dir <- "."
input_file <- file.path(
  project_dir, "results", "GSE306201", "Malignant",
  "Malignant_Epithelial_Reclustered.rds"
)
out_dir <- file.path(project_dir, "results", "GSE306201", "trajectory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_col <- "seurat_clusters"

# Load data
mal_obj <- readRDS(input_file)
DefaultAssay(mal_obj) <- "RNA"

# CytoTRACE2
mal_obj <- RunCytoTRACE(
  mal_obj,
  species = "Homo_sapiens"
)

p_cytotrace <- CytoTRACEPlot(
  mal_obj,
  group.by = cluster_col,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

ggsave(
  file.path(out_dir, "CytoTRACE2.pdf"),
  p_cytotrace,
  width = 12,
  height = 6
)

# Slingshot
mal_obj <- RunSlingshot(
  mal_obj,
  group.by = cluster_col,
  reduction = "umap"
)

n_lineages <- sum(grepl("^Lineage", colnames(mal_obj@meta.data)))
lineage_names <- paste0("Lineage", seq_len(n_lineages))

p_trajectory <- CellDimPlot(
  mal_obj,
  group.by = cluster_col,
  lineages = lineage_names,
  reduction = "umap",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
) +
  ggtitle("Slingshot trajectory") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_pseudotime <- FeatureDimPlot(
  mal_obj,
  features = lineage_names,
  reduction = "umap",
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

ggsave(
  file.path(out_dir, "Slingshot_trajectory.pdf"),
  p_trajectory,
  width = 8,
  height = 7
)

ggsave(
  file.path(out_dir, "Slingshot_pseudotime.pdf"),
  p_pseudotime,
  width = 12,
  height = 6
)

saveRDS(
  mal_obj,
  file.path(out_dir, "GSE306201_CytoTRACE2_Slingshot.rds")
)
