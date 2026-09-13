# hdWGCNA network construction

PROJECT_DIR <- "."

library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(UCell)
library(ggradar)

WORK_DIR <- file.path(PROJECT_DIR, "results", "GSE306201", "hdWGCNA")
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

scRNA <- file.path(PROJECT_DIR, "results", "GSE306201", "Malignant", "Malignant_Epithelial_Reclustered.rds")

dat_sub <- readRDS(scRNA)

dat_sub$C5_cluster_label <- paste0("C", as.character(dat_sub$seurat_clusters))

seurat_obj <- SetupForWGCNA(
  dat_sub,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "Malignant_C5"
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,

  group.by = c("orig.ident", "C5_cluster_label"),
  reduction = 'harmony',
  k = 15,
  max_shared = 5,
  ident.group = 'C5_cluster_label',
  min_cells = 20
)

seurat_obj <- NormalizeMetacells(seurat_obj)

all_clusters <- sort(unique(dat_sub$C5_cluster_label))

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = all_clusters,
  group.by = 'C5_cluster_label',
  assay = 'RNA',
  slot = 'data'
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  networkType = 'signed'
)

plot_list <- PlotSoftPowers(seurat_obj)

combined_plot <- wrap_plots(plot_list, ncol = 2)

ggsave(filename = "1.pdf", plot = combined_plot, width = 10, height = 8, units = "in", device = "pdf")

seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = NULL,
  minModuleSize = 30,
  deepSplit = 2,
  mergeCutHeight = 0.25,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  tom_outdir = "TOM_Malignant",
  tom_name = "Malignant_network",
  randomSeed = 2026,
  overwrite_tom = TRUE,
  useDiskCache = TRUE
)

pdf("2plot.pdf", width = 10, height = 8)
PlotDendrogram(seurat_obj, main = 'Malignant hdWGCNA Dendrogram')
dev.off()

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

seurat_obj <- ModuleEigengenes(seurat_obj)

hMEs <- GetMEs(seurat_obj)

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = NULL,
  group_name = NULL
)

module_plot <- PlotKMEs(seurat_obj, ncol = 4)

ggsave(filename = "3KMEplot.pdf", plot = module_plot, width = 12, height = 10, units = "in", device = "pdf")

modules <- GetModules(seurat_obj) %>%
  subset(module != 'grey')

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  wgcna_name = NULL,
  method = 'UCell'
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'hMEs',
  order = TRUE
)

combined_plot <- wrap_plots(plot_list, ncol = 4)

ggsave(filename = "4hMEsplot.pdf", plot = combined_plot, width = 16, height = 12, units = "in", device = "pdf", dpi = 300)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'scores',
  order = 'shuffle',
  ucell = TRUE
)

hub_plot <- wrap_plots(plot_list, ncol = 4)

ggsave(
  filename = "5Hubgeneplot.pdf",
  plot = hub_plot,
  width = 16,
  height = 12,
  units = "in",
  device = "pdf"
)

radar_plot <- ModuleRadarPlot(
  seurat_obj,
  group.by = 'C5_cluster_label',
  features = "hMEs",
  barcodes = NULL,
  wgcna_name = NULL,
  fill = TRUE,
  draw.points = FALSE,
  grid.label.size = 4
)

ggsave(
  filename = "6plot.pdf",
  plot = radar_plot,
  width = 20,
  height = 20,
  units = "in",
  device = "pdf",
  dpi = 300
)

pdf("7plot.pdf", width = 10, height = 10)
ModuleCorrelogram(seurat_obj)
dev.off()

MEs <- GetMEs(seurat_obj, harmonized = TRUE)

modules <- GetModules(seurat_obj)

mods <- levels(modules$module)
mods <- mods[mods != 'grey']

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- DotPlot(
  seurat_obj,
  features = mods,
  group.by = 'C5_cluster_label'
)

p <- p +
  RotatedAxis() +
  scale_color_gradient2(
    high = '#8B0000',
    mid = 'grey95',
    low = '#008B8B'
  ) +
  labs(
    x = "Modules",
    y = "Malignant cluster",
    title = "Module Eigengene Expression",
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

ggsave(
  filename = "8MEsplot.pdf",
  plot = p,
  width = 12,
  height = 8,
  units = "in",
  device = "pdf",
  dpi = 300
)

if (!dir.exists('ModuleNetworks')) {
  dir.create('ModuleNetworks', recursive = TRUE)
}

ModuleNetworkPlot(
  seurat_obj,
  outdir = 'ModuleNetworks',
  n_inner = 20,
  n_outer = 30,
  n_conns = Inf,
  plot_size = c(10, 10),
  vertex.label.cex = 1
)

modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

options(future.globals.maxSize = 10 * 1024^3)
pdf("9A_Hubgenenetwork_.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2,
  n_other = 2,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

selected_mods <- mods[1:min(5, length(mods))]

pdf("9B_Hubgenenetwork_.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10,
  n_other = 10,
  edge_prop = 0.75,
  mods = selected_mods
)
dev.off()

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  n_neighbors = 15,
  min_dist = 0.1
)

umap_df <- GetModuleUMAP(seurat_obj)

p_ggplot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(
    color = umap_df$color,
    size = umap_df$kME * 2
  ) +
  umap_theme() +
  labs(
    title = "Module UMAP (ggplot2)",
    subtitle = "Point size = connectivity (kME), Color = module"
  )

ggsave(
  filename = "10A_UMAP_ggplot2.pdf",
  plot = p_ggplot,
  width = 10,
  height = 10,
  units = "in",
  device = "pdf",
  dpi = 300
)

pdf("10B_module_UMAP.pdf", width = 10, height = 10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha = 0.25,
  sample_edges = TRUE,
  edge_prop = 0.1,
  label_hubs = 2,
  keep_grey_edges = FALSE
)
dev.off()

saveRDS(seurat_obj, "hdWGCNA_Malignant.rds")

library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(UCell)
library(ggradar)

WORK_DIR <- file.path(PROJECT_DIR, "results", "GSE306201", "hdWGCNA")
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

seurat_obj <- readRDS("hdWGCNA_Malignant.rds")

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj$C5_cluster_label <- paste0("C", as.character(seurat_obj$seurat_clusters))

module_genes <- GetModules(seurat_obj) %>%
  pull(gene_name) %>%
  unique()

seurat_obj <- ScaleData(
  seurat_obj,
  features = module_genes,
  verbose = FALSE
)

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars = "orig.ident",
  assay = "RNA",
  verbose = FALSE
)

hMEs <- GetMEs(seurat_obj)

all_clusters <- sort(unique(seurat_obj$C5_cluster_label))

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "C5_cluster_label",
  group_name = all_clusters,
  harmonized = TRUE,
  assay = "RNA"
)

module_plot <- PlotKMEs(seurat_obj, ncol = 4)
ggsave(filename = "3KMEplot.pdf", plot = module_plot,
       width = 12, height = 10, units = "in", device = "pdf")

modules <- GetModules(seurat_obj) %>%
  subset(module != 'grey')

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  wgcna_name = NULL,
  method = 'UCell'
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'hMEs',
  order = TRUE
)
combined_plot <- wrap_plots(plot_list, ncol = 4)
ggsave(filename = "4hMEsplot.pdf", plot = combined_plot,
       width = 16, height = 12, units = "in", device = "pdf", dpi = 300)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features = 'scores',
  order = 'shuffle',
  ucell = TRUE
)
hub_plot <- wrap_plots(plot_list, ncol = 4)
ggsave(filename = "5Hubgeneplot.pdf", plot = hub_plot,
       width = 16, height = 12, units = "in", device = "pdf")

radar_plot <- ModuleRadarPlot(
  seurat_obj,
  group.by = 'C5_cluster_label',
  features = "hMEs",
  barcodes = NULL,
  wgcna_name = NULL,
  fill = TRUE,
  draw.points = FALSE,
  grid.label.size = 4
)
ggsave(filename = "6plot.pdf", plot = radar_plot,
       width = 20, height = 20, units = "in", device = "pdf", dpi = 300)

pdf("7plot.pdf", width = 10, height = 10)
ModuleCorrelogram(seurat_obj)
dev.off()

MEs <- GetMEs(seurat_obj, harmonized = TRUE)

modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

seurat_obj@meta.data <- seurat_obj@meta.data[
  , !(colnames(seurat_obj@meta.data) %in% colnames(MEs)), drop = FALSE
]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- DotPlot(
  seurat_obj,
  features = mods,
  group.by = 'C5_cluster_label'
)

p <- p +
  RotatedAxis() +
  scale_color_gradient2(
    high = '#8B0000',
    mid = 'grey95',
    low = '#008B8B'
  ) +
  labs(
    x = "Modules",
    y = "Malignant cluster",
    title = "Module Eigengene Expression",
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

ggsave(
  filename = "8MEsplot.pdf",
  plot = p,
  width = 12,
  height = 8,
  units = "in",
  device = "pdf",
  dpi = 300
)

if (!dir.exists('ModuleNetworks')) {
  dir.create('ModuleNetworks', recursive = TRUE)
}

ModuleNetworkPlot(
  seurat_obj,
  outdir = 'ModuleNetworks',
  n_inner = 20,
  n_outer = 30,
  n_conns = Inf,
  plot_size = c(10, 10),
  vertex.label.cex = 1
)

modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

options(future.globals.maxSize = 10 * 1024^3)

pdf("9A_Hubgenenetwork_.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 2,
  n_other = 2,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

selected_mods <- mods[1:min(5, length(mods))]

pdf("9B_Hubgenenetwork_.pdf", width = 10, height = 10)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10,
  n_other = 10,
  edge_prop = 0.75,
  mods = selected_mods
)
dev.off()

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  n_neighbors = 15,
  min_dist = 0.1
)

umap_df <- GetModuleUMAP(seurat_obj)

p_ggplot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(
    color = umap_df$color,
    size = umap_df$kME * 2
  ) +
  umap_theme() +
  labs(
    title = "Module UMAP (ggplot2)",
    subtitle = "Point size = connectivity (kME), Color = module"
  )

ggsave(
  filename = "10A_UMAP_ggplot2.pdf",
  plot = p_ggplot,
  width = 10,
  height = 10,
  units = "in",
  device = "pdf",
  dpi = 300
)

pdf("10B_module_UMAP.pdf", width = 10, height = 10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha = 0.25,
  sample_edges = TRUE,
  edge_prop = 0.1,
  label_hubs = 2,
  keep_grey_edges = FALSE
)
dev.off()

saveRDS(seurat_obj, "hdWGCNA_Malignant.rds")
