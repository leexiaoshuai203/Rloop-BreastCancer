# Malignant epithelial markers and GO analysis

PROJECT_DIR <- "."

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

setwd(file.path(PROJECT_DIR, "results/GSE306201"))

STEP7_RDS      <- "Step7_Epithelial_PhenoIter_SCIPACseed.rds"
FULL_OBJ_RDS   <- "GSE306201_Step4_the_end_Annotated.rds"

OUT_P7  <- "Malignant/SCOP/Part7_DEG_GO"
dir.create(OUT_P7, showWarnings = FALSE, recursive = TRUE)

MAL_HVG   <- 2000L
MAL_PC    <- 20L
MAL_RES   <- 0.1
BATCH_COL <- "orig.ident"
CLUSTER_COL <- "seurat_clusters"

set.seed(1234)

mal_obj_path <- file.path(file.path(PROJECT_DIR, "results/GSE306201/Malignant"), "Malignant_Epithelial_Reclustered.rds")

if (!file.exists(mal_obj_path)) {
  stop(paste0(": ", mal_obj_path,
              ", ."))
}

mal_obj <- readRDS(mal_obj_path)
DefaultAssay(mal_obj) <- "RNA"

print(table(mal_obj$seurat_clusters))

Idents(mal_obj) <- CLUSTER_COL
mal_markers <- FindAllMarkers(
  mal_obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox",
  verbose         = FALSE
)
write.csv(mal_markers, file.path(OUT_P7, "All_Markers.csv"), row.names = FALSE)

marker_count <- mal_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(MarkerNumber = n(), .groups = "drop") %>%
  mutate(CellType = paste0("C", cluster))

go_results_list <- list()

for (cl in sort(unique(mal_markers$cluster))) {

  genes <- mal_markers %>%
    dplyr::filter(cluster == cl, p_val_adj < 0.05) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
    dplyr::pull(gene)

  genes <- unique(genes)

  if (length(genes) < 10) {
    next
  }

  eg <- tryCatch(
    bitr(
      genes,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db,
      drop     = TRUE
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(eg) || nrow(eg) == 0) {
    next
  }

  go_res <- tryCatch(
    enrichGO(
      gene          = unique(eg$ENTREZID),
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(go_res)) {
    next
  }

  go_df <- as.data.frame(go_res)

  if (nrow(go_df) == 0) {
    next
  }

  top3 <- go_df %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 3) %>%
    dplyr::mutate(
      Description = ifelse(
        nchar(Description) > 45,
        paste0(substr(Description, 1, 42), "..."),
        Description
      )
    )

  if (nrow(top3) > 0) {
    go_results_list[[paste0("C", cl)]] <- data.frame(
      CellType = paste0("C", cl),
      GOItem   = top3$Description,
      Pvalue   = top3$p.adjust,
      Rank     = seq_len(nrow(top3)),
      stringsAsFactors = FALSE
    )
  } else {
  }
}

go_top1 <- bind_rows(go_results_list) %>% filter(Rank == 1)

go_all <- bind_rows(go_results_list)
write.csv(go_all, file.path(OUT_P7, "GO_BP_Top3_per_cluster.csv"), row.names = FALSE)

plot_data <- left_join(marker_count, go_top1, by = "CellType") %>% drop_na()
plot_data <- plot_data %>% mutate(logP = -log10(Pvalue))

ct_order  <- plot_data$CellType
go_order  <- plot_data$GOItem
plot_data$CellType <- factor(plot_data$CellType, levels = rev(ct_order))
plot_data$GOItem   <- factor(plot_data$GOItem,   levels = rev(go_order))

n_cluster <- nrow(plot_data)
color_pool <- c("#D65656","#5FAB5F","#DDB370","#E1A99A","#CFE6A1",
                "#72C3E3","#D43B63","#D796C3","#9B3683","#294B2E",
                "#B49267","#E69F2F","#66A4B8","#5BC5BA","#EB2B1E")
Color <- rev(color_pool[1:n_cluster])

p1 <- ggplot(plot_data, aes(x = MarkerNumber, y = CellType)) +
  geom_segment(aes(y = CellType, yend = CellType, x = 0, xend = MarkerNumber),
               linewidth = 0.8, color = "grey80") +
  geom_point(aes(color = CellType), size = 6) +
  scale_color_manual(values = Color) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, max(plot_data$MarkerNumber) * 1.25)) +
  labs(title = "Marker gene numbers per cluster", x = NULL, y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text     = element_text(color = "black"),
        panel.grid    = element_blank(),
        legend.position = "none")

p2 <- ggplot(plot_data, aes(y = GOItem)) +
  geom_bar(aes(x = logP, fill = GOItem),
           stat = "identity", width = 0.55,
           color = "transparent", alpha = 0.8) +
  geom_text(aes(x = 0.15, label = GOItem), hjust = 0, size = 4.2) +
  labs(title = "Top GO term (BP) per cluster",
       x = "-log10(adj. P)", y = NULL) +
  scale_fill_manual(values = Color) +
  theme_bw(base_size = 14) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x   = element_text(color = "black"),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.grid    = element_blank(),
        legend.position = "none")

p_v1 <- p1 + p2 + plot_layout(widths = c(1, 1.8))
ggsave(file.path(OUT_P7, "Version1_Lollipop_GO.pdf"),
       plot = p_v1, height = 6, width = 6)

p3 <- ggplot(plot_data, aes(x = 0, y = CellType)) +
  geom_point(aes(color = CellType, size = MarkerNumber)) +
  scale_color_manual(values = Color, guide = "none") +
  scale_size_continuous(name = "Marker\nCount",
                        range = c(3, 10),
                        breaks = pretty(plot_data$MarkerNumber, n = 4)) +
  scale_x_continuous(expand = expansion(mult = 0.8)) +
  labs(title = "", x = NULL, y = NULL) +
  theme_bw(base_size = 14) +
  theme(panel.border   = element_blank(),
        axis.text.x    = element_blank(),
        axis.text.y    = element_text(hjust = 1, color = "black"),
        axis.ticks     = element_blank(),
        panel.grid     = element_blank(),
        legend.position = "right")

p_v2 <- p3 + p2 + plot_layout(widths = c(0.18, 1.3), guides = "collect")
ggsave(file.path(OUT_P7, "Version2_Bubble_GO.pdf"),
       plot = p_v2, height = 6, width = 6)

suppressPackageStartupMessages({
  library(Seurat)
  library(ClusterGVis)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
})

setwd(file.path(PROJECT_DIR, "results/GSE306201"))

MAL_OBJ_RDS <- "Malignant/Malignant_Epithelial_Reclustered.rds"
OUT_P7      <- "Malignant/SCOP/Part7_DEG_GO"
MARKERS_CSV <- file.path(OUT_P7, "All_Markers.csv")
dir.create(OUT_P7, showWarnings = FALSE, recursive = TRUE)

set.seed(1234)

mal_obj <- readRDS(MAL_OBJ_RDS)
DefaultAssay(mal_obj) <- "RNA"
Idents(mal_obj) <- "seurat_clusters"

if (file.exists(MARKERS_CSV)) {
  mal_markers <- read.csv(MARKERS_CSV, stringsAsFactors = FALSE)
} else {
  mal_markers <- FindAllMarkers(
    mal_obj,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    test.use        = "wilcox",
    verbose         = FALSE
  )
  write.csv(mal_markers, MARKERS_CSV, row.names = FALSE)
}

mal_markers$cluster <- factor(mal_markers$cluster,
                              levels = sort(unique(as.numeric(as.character(mal_markers$cluster)))))

top_markers <- mal_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC) %>%
  dplyr::ungroup()

print(table(top_markers$cluster))

st.data <- prepareDataFromscRNA(
  object       = mal_obj,
  diffData     = top_markers,
  showAverage  = TRUE,
  keep.uniqGene = TRUE
)

str(st.data, max.level = 1)

enrich <- enrichCluster(
  object        = st.data,
  OrgDb         = org.Hs.eg.db,
  type          = "BP",
  organism      = "hsa",
  pvalueCutoff  = 0.5,
  topn          = 5,
  seed          = 5201314
)

head(enrich)
write.csv(enrich, file.path(OUT_P7, "ClusterGVis_GO_BP_Top5.csv"), row.names = FALSE)

markGenes <- top_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 4, wt = avg_log2FC) %>%
  dplyr::pull(gene) %>%
  unique()

print(markGenes)

n_cluster <- length(unique(top_markers$cluster))
cluster_order <- seq_len(n_cluster)

pdf(file.path(OUT_P7, "ClusterGVis_Fig1_LinePlot.pdf"),
    width = 10, height = 7, onefile = FALSE)
visCluster(
  object    = st.data,
  plot.type = "line"
)
dev.off()

pdf(file.path(OUT_P7, "ClusterGVis_Fig2_Heatmap.pdf"),
    width = 7, height = 11, onefile = FALSE)
visCluster(
  object           = st.data,
  plot.type        = "heatmap",
  column_names_rot = 45,
  markGenes        = markGenes,
  cluster.order    = cluster_order,
  ctAnno.col       = ggsci::pal_npg()(n_cluster)
)
dev.off()

go_col <- rep(ggsci::pal_d3()(n_cluster), each = 5)

n_per_cluster <- enrich %>% dplyr::count(group) %>% dplyr::pull(n)
go_col <- unlist(mapply(function(col, n) rep(col, n),
                        ggsci::pal_d3()(n_cluster)[seq_along(n_per_cluster)],
                        n_per_cluster, SIMPLIFY = FALSE))

pdf(file.path(OUT_P7, "ClusterGVis_Fig3_Heatmap_with_GO.pdf"),
    width = 16, height = 11, onefile = FALSE)
visCluster(
  object          = st.data,
  plot.type       = "both",
  column_names_rot = 45,
  show_row_dend   = FALSE,
  markGenes       = markGenes,
  markGenes.side  = "left",
  annoTerm.data   = enrich,
  line.side       = "left",
  cluster.order   = cluster_order,
  go.col          = go_col,
  textbar.pos     = c(0.8, 0.2),
  add.bar         = TRUE,
  ctAnno.col      = ggsci::pal_npg()(n_cluster)
)
dev.off()

pdf(file.path(OUT_P7, "ClusterGVis_Fig4_Heatmap_with_GOterm.pdf"),
    width = 10, height = 11, onefile = FALSE)
visCluster(
  object          = st.data,
  plot.type       = "both",
  column_names_rot = 45,
  show_row_dend   = FALSE,
  markGenes       = markGenes,
  markGenes.side  = "left",
  annoTerm.data   = enrich,
  line.side       = "left",
  cluster.order   = cluster_order,
  go.col          = go_col,
  add.bar         = FALSE,
  ctAnno.col      = ggsci::pal_npg()(n_cluster)
)
dev.off()
