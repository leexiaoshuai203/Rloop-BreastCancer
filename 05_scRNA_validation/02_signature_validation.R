# RL-Sig65, CIN70 and C5-program validation

PROJECT_DIR <- "."

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(UCell)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(
  base_dir, "10.malignat",
  "GSE245601_Malignant_Epithelial_Reclustered.rds"
)
rl65_file <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")
out_dir <- file.path(base_dir, "11.RL65&CIN70")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS <- 20L
UCELL_MAXRANK <- 5000L

for (f in c(input_file, rl65_file)) {
  if (!file.exists(f)) stop("not_found: ", f)
}

cin70 <- c(
  "TPX2","PRC1","FOXM1","CDK1","TGIF2","MCM2","H2AFZ","TOP2A","PCNA","UBE2C",
  "MELK","TRIP13","NCAPD2","MCM7","RNASEH2A","RAD51AP1","KIF20A","CDC45",
  "MAD2L1","ESPL1","CCNB2","FEN1","TTK","CCT5","RFC4","ATAD2","CKAP5",
  "NUP205","CDC20","CKS2","RRM2","ELAVL1","CCNB1","RRM1","AURKB","MSH6",
  "EZH2","CTPS1","DKC1","OIP5","CDCA8","PTTG1","CEP55","H2AFX","CMAS",
  "NCAPH","MCM10","LSM4","NCAPG2","ASF1B","ZWINT","PBK","CDCA3","ECT2",
  "CDC6","UNG","MTCH2","RAD21","ACTL6A","GPI","SRSF2","HDGF","NXT1",
  "NEK2","DHCR7","AURKA","NDUFAB1","MIIP","KIF4A"
)

mal_obj <- readRDS(input_file)
DefaultAssay(mal_obj) <- "RNA"

sample_col <- if ("sample_gsm" %in% colnames(mal_obj[[]])) {
  "sample_gsm"
} else {
  "orig.ident"
}

if (!sample_col %in% colnames(mal_obj[[]])) stop("missingsample")
if (!"C_cluster" %in% colnames(mal_obj[[]])) stop("missingC_cluster")

if (length(Layers(mal_obj[["RNA"]])) > 1) {
  mal_obj <- JoinLayers(mal_obj, assay = "RNA")
}

cluster_levels <- paste0(
  "C",
  sort(unique(as.integer(gsub("C", "", as.character(mal_obj$C_cluster)))))
)
mal_obj$C_cluster <- factor(
  as.character(mal_obj$C_cluster),
  levels = cluster_levels
)

print(table(mal_obj$C_cluster))

meta_df <- data.frame(
  Cell = colnames(mal_obj),
  Sample = as.character(mal_obj[[sample_col, drop = TRUE]]),
  Cluster = as.character(mal_obj$C_cluster),
  stringsAsFactors = FALSE
)

pb_map <- meta_df %>%
  count(Sample, Cluster, name = "CellNumber") %>%
  arrange(Sample, factor(Cluster, levels = cluster_levels)) %>%
  mutate(
    Valid = CellNumber >= MIN_CELLS,
    PB_id = sprintf("PB%03d", row_number())
  )

write.csv(
  pb_map,
  file.path(out_dir, "00_bulkgroupcell.csv"),
  row.names = FALSE
)

valid_map <- pb_map %>% filter(Valid)
id_key <- setNames(
  valid_map$PB_id,
  paste(valid_map$Sample, valid_map$Cluster, sep = "||")
)

cell_key <- paste(meta_df$Sample, meta_df$Cluster, sep = "||")
mal_obj$PB_id <- unname(id_key[cell_key])
mal_obj <- subset(mal_obj, subset = !is.na(PB_id))
mal_obj$PB_id <- factor(mal_obj$PB_id, levels = valid_map$PB_id)

pb_counts <- AggregateExpression(
  mal_obj,
  assays = "RNA",
  group.by = "PB_id",
  slot = "counts",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

if (!setequal(colnames(pb_counts), valid_map$PB_id)) {
  stop("PBAggregateExpressionresults")
}

lib_size <- Matrix::colSums(pb_counts)
pb_expr <- log1p(t(t(pb_counts) / lib_size) * 10000)

rl65_raw <- unique(na.omit(trimws(readLines(rl65_file, warn = FALSE))))
rl65_raw <- rl65_raw[rl65_raw != ""]

rl65_use <- intersect(rl65_raw, rownames(pb_expr))
cin70_use <- intersect(cin70, rownames(pb_expr))

if (length(rl65_use) < 5) stop("RL65gene5")
if (length(cin70_use) < 5) stop("CIN70gene5")

write.csv(
  data.frame(Gene = rl65_use),
  file.path(out_dir, "01_RL65gene.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Gene = cin70_use),
  file.path(out_dir, "02_CIN70gene.csv"),
  row.names = FALSE
)

score_raw <- UCell::ScoreSignatures_UCell(
  matrix = pb_expr,
  features = list(RL65 = rl65_use, CIN70 = cin70_use),
  maxRank = min(UCELL_MAXRANK, nrow(pb_expr))
)
score_raw <- as.data.frame(score_raw, check.names = FALSE)

rl_col <- grep("RL65", colnames(score_raw), ignore.case = TRUE, value = TRUE)[1]
cin_col <- grep("CIN70", colnames(score_raw), ignore.case = TRUE, value = TRUE)[1]

if (is.na(rl_col) || is.na(cin_col)) {
  stop("not_foundUCell: ", paste(colnames(score_raw), collapse = ", "))
}

score_df <- data.frame(
  PB_id = rownames(score_raw),
  RL65 = as.numeric(score_raw[, rl_col, drop = TRUE]),
  CIN70 = as.numeric(score_raw[, cin_col, drop = TRUE])
)

pb_df <- valid_map %>%
  select(PB_id, Sample, Cluster, CellNumber) %>%
  left_join(score_df, by = "PB_id")

if (anyNA(pb_df)) stop("PBNA, checkPB_id")
pb_df$Cluster <- factor(pb_df$Cluster, levels = cluster_levels)

write.csv(
  pb_df,
  file.path(out_dir, "03_Pseudobulk_RL65_CIN70_Scores.csv"),
  row.names = FALSE
)

cluster_summary <- pb_df %>%
  group_by(Cluster) %>%
  summarise(
    N_samples = n(),
    Total_cells = sum(CellNumber),
    RL65_mean = mean(RL65),
    RL65_median = median(RL65),
    CIN70_mean = mean(CIN70),
    CIN70_median = median(CIN70),
    .groups = "drop"
  )

write.csv(
  cluster_summary,
  file.path(out_dir, "04_clustersummary.csv"),
  row.names = FALSE
)

cor_res <- cor.test(
  pb_df$RL65,
  pb_df$CIN70,
  method = "spearman",
  exact = FALSE
)

write.csv(
  data.frame(
    Analysis = "Overall_Spearman",
    Spearman_rho = unname(cor_res$estimate),
    P_value = cor_res$p.value,
    N = nrow(pb_df)
  ),
  file.path(out_dir, "05_statistics.csv"),
  row.names = FALSE
)

cluster_cols <- setNames(
  scales::hue_pal()(length(cluster_levels)),
  cluster_levels
)

long_df <- pb_df %>%
  pivot_longer(
    c(RL65, CIN70),
    names_to = "Signature",
    values_to = "Score"
  )

p_box <- ggplot(long_df, aes(Cluster, Score)) +
  geom_boxplot(
    aes(fill = Cluster),
    width = 0.55,
    alpha = 0.75,
    outlier.shape = NA
  ) +
  geom_point(
    aes(fill = Cluster),
    shape = 21,
    color = "white",
    size = 2.6,
    position = position_jitter(width = 0.08)
  ) +
  facet_wrap(~Signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = cluster_cols) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "RL65 and CIN70 across malignant epithelial clusters",
    x = "Cluster",
    y = "Pseudo-bulk UCell score"
  )

p_scatter <- ggplot(pb_df, aes(RL65, CIN70, color = Cluster)) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "grey40",
    fill = "grey85",
    linetype = "dashed"
  ) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = paste0(Sample, "_", Cluster)),
    size = 2.2,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_cols) +
  annotate(
    "text",
    x = min(pb_df$RL65),
    y = max(pb_df$CIN70),
    hjust = 0,
    vjust = 1,
    label = paste0(
      "Spearman rho = ",
      round(unname(cor_res$estimate), 3),
      "\np = ",
      signif(cor_res$p.value, 3)
    ),
    size = 4
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "RL65 versus CIN70",
    x = "RL65 UCell score",
    y = "CIN70 UCell score",
    color = "Cluster"
  )

ggsave(
  file.path(out_dir, "06_RL65_CIN70resultsplot.pdf"),
  p_box | p_scatter,
  width = 15,
  height = 5.5
)

p_heatmap <- ggplot(
  long_df,
  aes(Cluster, Sample, fill = Score)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_wrap(~Signature, scales = "free", nrow = 1) +
  scale_fill_gradientn(colors = c("#6fa6cf", "white", "#b43665")) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(title = "Sample-cluster pseudo-bulk scores")

ggsave(
  file.path(out_dir, "07_RL65_CIN70plot.pdf"),
  p_heatmap,
  width = 10,
  height = 6
)

saveRDS(
  list(
    scores = pb_df,
    cluster_summary = cluster_summary,
    correlation = cor_res,
    RL65_genes = rl65_use,
    CIN70_genes = cin70_use
  ),
  file.path(out_dir, "Pseudobulk_RL65_CIN70_Result.rds")
)

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ClusterGVis)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")

input_file <- file.path(
  base_dir, "10.malignat",
  "GSE245601_Malignant_Epithelial_Reclustered.rds"
)

out_dir <- file.path(
  base_dir, "11.RL65&CIN70",
  "11B_DEG_GO"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TOP_HEATMAP <- 20L
LABEL_N <- 4L
ENRICH_TOP_GENES <- 300L
GO_TOP_PLOT <- 5L
GO_TOP_TXT <- 10L

if (!file.exists(input_file)) {
  stop("not_foundinput: ", input_file)
}

mal_obj <- readRDS(input_file)
DefaultAssay(mal_obj) <- "RNA"

if (!"seurat_clusters" %in% colnames(mal_obj[[]])) {
  stop("missingseurat_clusters")
}

if (length(Layers(mal_obj[["RNA"]])) > 1) {
  mal_obj <- JoinLayers(mal_obj, assay = "RNA")
}

if (!"data" %in% Layers(mal_obj[["RNA"]])) {
  mal_obj <- NormalizeData(mal_obj, verbose = FALSE)
}

cluster_ids <- sort(unique(as.integer(
  as.character(mal_obj$seurat_clusters)
)))

cluster_levels <- as.character(cluster_ids)

mal_obj$seurat_clusters <- factor(
  as.character(mal_obj$seurat_clusters),
  levels = cluster_levels
)

Idents(mal_obj) <- "seurat_clusters"

print(table(mal_obj$seurat_clusters))

markers <- FindAllMarkers(
  mal_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = FALSE
)

markers <- markers %>%
  dplyr::filter(
    p_val_adj < 0.05,
    !is.na(gene),
    gene != ""
  )

if (nrow(markers) == 0) {
  stop("Marker")
}

markers$cluster <- factor(
  as.character(markers$cluster),
  levels = cluster_levels
)

marker_count <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    MarkerNumber = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    CellType = as.character(cluster)
  )

top_markers_plot <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = TOP_HEATMAP,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

top_markers_enrich <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = ENRICH_TOP_GENES,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

st_enrich <- prepareDataFromscRNA(
  object = mal_obj,
  diffData = top_markers_enrich,
  showAverage = TRUE,
  keep.uniqGene = TRUE
)

enrich_internal <- enrichCluster(
  object = st_enrich,
  OrgDb = org.Hs.eg.db,
  type = "BP",
  organism = "hsa",
  pvalueCutoff = 0.05,
  topn = GO_TOP_TXT,
  seed = 1234
)

if (is.null(enrich_internal) || nrow(enrich_internal) == 0) {
  stop("GO BPenrichmentresults")
}

required_cols <- c("group", "Description", "pvalue")

if (!all(required_cols %in% colnames(enrich_internal))) {
  stop(
    "GOresultsmissing, : ",
    paste(colnames(enrich_internal), collapse = ", ")
  )
}

enrich_internal$group <- as.character(enrich_internal$group)

internal_groups <- unique(enrich_internal$group)

internal_groups <- internal_groups[
  order(as.integer(sub("^C", "", internal_groups)))
]

if (length(internal_groups) != length(cluster_levels)) {
  stop(
    "ClusterGVisgroupcluster.: ",
    paste(internal_groups, collapse = ", "),
    "; : ",
    paste(cluster_levels, collapse = ", ")
  )
}

group_map <- setNames(
  cluster_levels,
  internal_groups
)

enrich_mapped <- enrich_internal
enrich_mapped$group <- unname(
  group_map[enrich_mapped$group]
)

if (anyNA(enrich_mapped$group)) {
  stop("ClusterGVisgroupNA")
}

enrich_internal <- enrich_internal %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(pvalue, .by_group = TRUE) %>%
  dplyr::slice_head(n = GO_TOP_TXT) %>%
  dplyr::ungroup()

enrich_mapped <- enrich_mapped %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(pvalue, .by_group = TRUE) %>%
  dplyr::slice_head(n = GO_TOP_TXT) %>%
  dplyr::ungroup() %>%
  dplyr::rename(cluster = group) %>%
  dplyr::arrange(
    factor(cluster, levels = cluster_levels),
    pvalue
  )

write.table(
  enrich_mapped,
  file.path(out_dir, "11B_GO_BP_Top10_per_cluster.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

go_top1 <- enrich_mapped %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_min(
    order_by = pvalue,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    CellType = as.character(cluster),
    GOItem = Description,
    Pvalue = pvalue
  )

plot_data <- marker_count %>%
  dplyr::left_join(
    go_top1,
    by = "CellType"
  ) %>%
  tidyr::drop_na(GOItem, Pvalue) %>%
  dplyr::mutate(
    logP = -log10(Pvalue)
  )

if (nrow(plot_data) == 0) {
  stop("plot1GOresults")
}

cell_order <- rev(
  cluster_levels[cluster_levels %in% plot_data$CellType]
)

plot_data$CellType <- factor(
  plot_data$CellType,
  levels = cell_order
)

cluster_cols <- setNames(
  ggsci::pal_npg()(length(cluster_levels)),
  cluster_levels
)

p1 <- ggplot(
  plot_data,
  aes(MarkerNumber, CellType)
) +
  geom_segment(
    aes(
      x = 0,
      xend = MarkerNumber,
      yend = CellType
    ),
    color = "grey80",
    linewidth = 0.8
  ) +
  geom_point(
    aes(color = CellType),
    size = 5
  ) +
  scale_color_manual(values = cluster_cols) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(
      0,
      max(plot_data$MarkerNumber) * 1.20
    )
  ) +
  labs(
    title = "Marker genes",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank(),
    legend.position = "none"
  )

p2 <- ggplot(
  plot_data,
  aes(logP, CellType)
) +
  geom_col(
    aes(fill = CellType),
    width = 0.55,
    alpha = 0.82
  ) +
  geom_text(
    aes(
      x = 0.10,
      label = GOItem
    ),
    hjust = 0,
    size = 3.3
  ) +
  scale_fill_manual(values = cluster_cols) +
  labs(
    title = "Top GO term",
    x = "-log10(P)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "11B_Fig1_Marker_GO.pdf"),
  p1 + p2 + plot_layout(widths = c(1, 1.8)),
  width = 4,
  height = 4
)

st_plot <- prepareDataFromscRNA(
  object = mal_obj,
  diffData = top_markers_plot,
  showAverage = TRUE,
  keep.uniqGene = TRUE
)

enrich_top5_internal <- enrich_internal %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(pvalue, .by_group = TRUE) %>%
  dplyr::slice_head(n = GO_TOP_PLOT) %>%
  dplyr::ungroup()

mark_genes <- top_markers_plot %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = LABEL_N,
    with_ties = FALSE
  ) %>%
  dplyr::pull(gene) %>%
  unique()

group_values <- unique(
  as.character(enrich_top5_internal$group)
)

group_map_color <- setNames(
  ggsci::pal_d3()(length(group_values)),
  group_values
)

go_cols <- unname(
  group_map_color[
    as.character(enrich_top5_internal$group)
  ]
)

pdf(
  file.path(
    out_dir,
    "11B_Fig2_ClusterGVis_Heatmap_GO.pdf"
  ),
  width = 11,
  height = 7,
  onefile = FALSE
)

visCluster(
  object = st_plot,
  plot.type = "both",
  column_names_rot = 45,
  show_row_dend = FALSE,
  markGenes = mark_genes,
  markGenes.side = "left",
  annoTerm.data = enrich_top5_internal,
  line.side = "left",
  cluster.order = seq_along(cluster_levels),
  go.col = go_cols,
  textbar.pos = c(0.8, 0.2),
  add.bar = TRUE,
  ctAnno.col = ggsci::pal_npg()(length(cluster_levels))
)

dev.off()

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
discovery_dir <- file.path(PROJECT_DIR, "results/GSE306201/supplementary")

input_file <- file.path(base_dir, "10.malignat",
                        "GSE245601_Malignant_Epithelial_Reclustered.rds")
c5_file <- file.path(discovery_dir, "02", "C5_stable_up_genes.txt")
rlcin_file <- file.path(base_dir, "11.RL65&CIN70",
                        "03_Pseudobulk_RL65_CIN70_Scores.csv")
out_dir <- file.path(base_dir, "11.RL65&CIN70", "11C_C5_Projection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS <- 20L
MAX_RANK <- 5000L

for (f in c(input_file, c5_file, rlcin_file)) {
  if (!file.exists(f)) stop("not_found: ", f)
}

safe_vec <- function(x) {
  if (is.data.frame(x)) x <- x[[1]]
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  as.character(x)
}

c5_genes <- unique(trimws(readLines(c5_file, warn = FALSE)))
c5_genes <- c5_genes[nzchar(c5_genes)]
if (length(c5_genes) < 5) stop("C5 signaturegene5")

obj <- readRDS(input_file)
DefaultAssay(obj) <- "RNA"

sample_col <- if ("sample_gsm" %in% colnames(obj[[]])) "sample_gsm" else "orig.ident"
if (!sample_col %in% colnames(obj[[]])) stop("missingsample")
if (!"C_cluster" %in% colnames(obj[[]])) stop("missingC_cluster")

if (length(Layers(obj[["RNA"]])) > 1) {
  obj <- JoinLayers(obj, assay = "RNA")
}

sample_vec <- safe_vec(obj@meta.data[[sample_col]])
cluster_vec <- safe_vec(obj@meta.data[["C_cluster"]])

if (length(sample_vec) != ncol(obj) || length(cluster_vec) != ncol(obj)) {
  stop("samplecluster")
}

cluster_levels <- paste0("C", sort(unique(as.integer(gsub("C", "", cluster_vec)))))
obj$C_cluster <- factor(cluster_vec, levels = cluster_levels)

print(table(obj$C_cluster))

pb_map <- as.data.frame(table(sample_vec, cluster_vec), stringsAsFactors = FALSE)
colnames(pb_map) <- c("Sample", "Cluster", "CellNumber")
pb_map <- pb_map[pb_map$CellNumber > 0, , drop = FALSE]
pb_map <- pb_map[order(pb_map$Sample, match(pb_map$Cluster, cluster_levels)), ]
pb_map$Valid <- pb_map$CellNumber >= MIN_CELLS
pb_map$PB_id <- sprintf("PB%03d", seq_len(nrow(pb_map)))
valid_map <- pb_map[pb_map$Valid, , drop = FALSE]

write.csv(pb_map, file.path(out_dir, "00_bulkgroupcell.csv"), row.names = FALSE)

id_map <- setNames(valid_map$PB_id,
                   paste(valid_map$Sample, valid_map$Cluster, sep = "||"))
obj$PB_id <- unname(id_map[paste(sample_vec, cluster_vec, sep = "||")])

obj <- subset(obj, cells = colnames(obj)[!is.na(obj$PB_id)])
obj$PB_id <- factor(obj$PB_id, levels = valid_map$PB_id)

pb_counts <- AggregateExpression(
  obj, assays = "RNA", group.by = "PB_id",
  slot = "counts", return.seurat = FALSE, verbose = FALSE
)$RNA

if (!setequal(colnames(pb_counts), valid_map$PB_id)) {
  stop("AggregateExpressionPB")
}

pb_counts <- pb_counts[, valid_map$PB_id, drop = FALSE]
lib_size <- Matrix::colSums(pb_counts)
if (any(lib_size <= 0)) stop("0bulk")
pb_expr <- log1p(t(t(pb_counts) / lib_size) * 10000)

c5_use <- intersect(c5_genes, rownames(pb_expr))
if (length(c5_use) < 5) stop("C5 signaturegene5")

write.csv(data.frame(Gene = c5_use),
          file.path(out_dir, "01_C5_signaturegene.csv"),
          row.names = FALSE)

score_raw <- UCell::ScoreSignatures_UCell(
  matrix = pb_expr,
  features = list(C5_SIGNATURE = c5_use),
  maxRank = min(MAX_RANK, nrow(pb_expr))
)
score_raw <- as.data.frame(score_raw, check.names = FALSE)

c5_col <- grep("C5_SIGNATURE", colnames(score_raw),
               ignore.case = TRUE, value = TRUE)[1]
if (is.na(c5_col)) stop("not_foundC5 signature")

c5_df <- data.frame(
  PB_id = rownames(score_raw),
  C5_SIGNATURE = as.numeric(score_raw[, c5_col, drop = TRUE])
)

projection_df <- dplyr::left_join(
  valid_map[, c("PB_id", "Sample", "Cluster", "CellNumber")],
  c5_df, by = "PB_id"
)
if (anyNA(projection_df)) stop("C5 signatureNA")

rlcin_df <- read.csv(rlcin_file, stringsAsFactors = FALSE, check.names = FALSE)
need_cols <- c("Sample", "Cluster", "RL65", "CIN70")
missing_cols <- setdiff(need_cols, colnames(rlcin_df))

if (length(missing_cols) > 0) {
  stop("RL65/CIN70missing: ", paste(missing_cols, collapse = ", "))
}

rlcin_df$Sample <- safe_vec(rlcin_df$Sample)
rlcin_df$Cluster <- safe_vec(rlcin_df$Cluster)

combined_df <- dplyr::left_join(
  projection_df,
  rlcin_df[, need_cols],
  by = c("Sample", "Cluster")
)

score_cols <- c("C5_SIGNATURE", "RL65", "CIN70")
if (anyNA(combined_df[, score_cols])) {
  stop("NA, 11A11CMIN_CELLS=20")
}

combined_df$Cluster <- factor(combined_df$Cluster, levels = cluster_levels)
write.csv(combined_df,
          file.path(out_dir, "02_signature_pseudobulk_results.csv"),
          row.names = FALSE)

summary_df <- combined_df %>%
  pivot_longer(all_of(score_cols), names_to = "Signature", values_to = "Score") %>%
  group_by(Cluster, Signature) %>%
  summarise(N_samples = n(), Mean = mean(Score), Median = median(Score),
            .groups = "drop") %>%
  group_by(Signature) %>%
  arrange(desc(Median), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

top_df <- summary_df %>%
  filter(Rank == 1) %>%
  select(Signature, TopCluster = Cluster, TopMedian = Median)

get_top <- function(x) {
  y <- as.character(top_df$TopCluster[top_df$Signature == x])
  if (length(y) != 1) stop("", x, "cluster")
  y
}

c5_top <- get_top("C5_SIGNATURE")
rl65_top <- get_top("RL65")
cin70_top <- get_top("CIN70")

decision <- if (all(c(c5_top, rl65_top, cin70_top) == "C4")) {
  "C4 showed the highest C5-signature, RL-Sig65 and CIN70 scores."
} else {
  paste0("C4: C5 signature=", c5_top,
         "; RL65=", rl65_top, "; CIN70=", cin70_top)
}

write.csv(summary_df,
          file.path(out_dir, "03_clustersummary.csv"),
          row.names = FALSE)
write.csv(top_df,
          file.path(out_dir, "04_cluster.csv"),
          row.names = FALSE)
writeLines(c(
  paste0("Top cluster for C5 signature: ", c5_top),
  paste0("Top cluster for RL-Sig65: ", rl65_top),
  paste0("Top cluster for CIN70: ", cin70_top),
  paste0("Conclusion: ", decision)
), file.path(out_dir, "05_C4_C5_similarity_summary.txt"))

plot_long <- combined_df %>%
  pivot_longer(all_of(score_cols), names_to = "Signature", values_to = "Score")

p_box <- ggplot(plot_long, aes(Cluster, Score, group = Sample)) +
  geom_line(color = "grey75", linewidth = 0.3, alpha = 0.5) +
  geom_boxplot(aes(fill = Cluster, group = Cluster),
               width = 0.55, alpha = 0.78, outlier.shape = NA) +
  geom_point(aes(fill = Cluster), shape = 21, color = "white", size = 2.3) +
  facet_wrap(~Signature, scales = "free_y", nrow = 1) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  labs(title = "Validation of the discovery C5-related state",
       x = "Validation cluster",
       y = "Pseudo-bulk UCell score")

ggsave(file.path(out_dir, "06_clustercomparison.pdf"),
       p_box, width = 11, height = 4.5)

heat_df <- summary_df %>%
  group_by(Signature) %>%
  mutate(Score_z = ifelse(sd(Median) > 0, as.numeric(scale(Median)), 0)) %>%
  ungroup()

heat_df$Signature <- factor(heat_df$Signature, levels = rev(score_cols))

p_heat <- ggplot(heat_df, aes(Cluster, Signature, fill = Score_z)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Median)), size = 3.6) +
  scale_fill_gradient2(low = "#6fa6cf", mid = "white", high = "#b43665",
                       midpoint = 0, name = "Row Z-score") +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "Median score by validation cluster")

ggsave(file.path(out_dir, "07_plot.pdf"),
       p_heat, width = 7, height = 4)

saveRDS(list(
  C5_genes = c5_use,
  scores = combined_df,
  summary = summary_df,
  decision = decision
), file.path(out_dir, "11C_C5_Projection_Result.rds"))
