# GSE245601 QC, integration and annotation

PROJECT_DIR <- "."

options(stringsAsFactors = FALSE)
set.seed(123)

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)

if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("hdf5r: install.packages('hdf5r')")
}

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
raw_root <- file.path(PROJECT_DIR, "data/raw/GSE245601")
raw_dir  <- file.path(raw_root, "GSE245601_RAW")
tar_file <- file.path(raw_root, "GSE245601_RAW.tar")
out_dir  <- file.path(base_dir, "01.QC")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

h5_files <- list.files(raw_dir, pattern = "\\.h5$", full.names = TRUE,
                       recursive = TRUE, ignore.case = TRUE)

if (length(h5_files) == 0) {
  if (!file.exists(tar_file)) stop("not_foundH5GSE245601_RAW.tar")
  utils::untar(tar_file, exdir = raw_dir)
  h5_files <- list.files(raw_dir, pattern = "\\.h5$", full.names = TRUE,
                         recursive = TRUE, ignore.case = TRUE)
}

if (length(h5_files) == 0) stop("not_foundH5")

sample_info <- data.frame(
  sample_gsm = paste0("GSM", seq(7845552, 7845570, by = 2)),
  patient_id = sprintf("Tumor_%02d", 1:10),
  sample_label = paste0(sprintf("Tumor_%02d", 1:10), "_Control"),
  stringsAsFactors = FALSE
)

print(sample_info)

match_h5_file <- function(gsm, sample_label, all_files) {
  file_names <- basename(all_files)
  hit_gsm <- grepl(gsm, file_names, fixed = TRUE)
  hit_label <- grepl(sample_label, file_names, fixed = TRUE, ignore.case = TRUE)
  matched <- all_files[hit_gsm | hit_label]

  if (length(matched) == 0) stop("not_foundsample: ", gsm, " / ", sample_label)
  if (length(matched) > 1) {
    print(matched)
    stop("checkH5")
  }
  matched
}

sample_info$h5_file <- mapply(
  FUN = match_h5_file,
  gsm = sample_info$sample_gsm,
  sample_label = sample_info$sample_label,
  MoreArgs = list(all_files = h5_files),
  USE.NAMES = FALSE
)

print(sample_info[, c("sample_gsm", "patient_id", "h5_file")])

read_one_h5 <- function(h5_file) {
  mat <- Read10X_h5(h5_file, use.names = TRUE, unique.features = TRUE)

  if (is.list(mat)) {
    if ("Gene Expression" %in% names(mat)) mat <- mat[["Gene Expression"]]
    else mat <- mat[[1]]
  }

  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  mat
}

seu_list <- vector("list", nrow(sample_info))
names(seu_list) <- sample_info$sample_gsm

for (i in seq_len(nrow(sample_info))) {
  sid <- sample_info$sample_gsm[i]

  counts <- read_one_h5(sample_info$h5_file[i])

  obj <- CreateSeuratObject(
    counts = counts, project = "GSE245601_ControlTumor",
    min.cells = 3, min.features = 200
  )

  obj$orig.ident  <- sid
  obj$sample_gsm  <- sid
  obj$patient_id  <- sample_info$patient_id[i]
  obj$sample_name <- sample_info$sample_label[i]
  obj$sample_type <- "Primary_tumor"
  obj$treatment   <- "Control"
  obj$dataset     <- "GSE245601"

  seu_list[[i]] <- obj
  rm(counts, obj)
  gc()
}

print(sapply(seu_list, ncol))

seurat_obj <- merge(
  x = seu_list[[1]], y = seu_list[-1],
  add.cell.ids = names(seu_list),
  project = "GSE245601_ControlTumor"
)

rm(seu_list)
gc()

print(table(seurat_obj$sample_gsm))

rna_layers <- Layers(seurat_obj[["RNA"]])

if (length(rna_layers) > 1) {
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
}

DefaultAssay(seurat_obj) <- "RNA"
gene_names <- rownames(seurat_obj)

n_mt_upper <- sum(grepl("^MT-", gene_names))
n_mt_lower <- sum(grepl("^mt-", gene_names))
mt_pattern <- if (n_mt_upper >= n_mt_lower) "^MT-" else "^mt-"

if (max(n_mt_upper, n_mt_lower) == 0) {
  stop("gene, checkH5gene")
}

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(RPL|RPS)")

p_before <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  group.by = "patient_id", pt.size = 0, ncol = 2
)

pdf(file.path(out_dir, "01_pre_QCQCplot.pdf"), width = 14, height = 9)
print(p_before)
dev.off()
rm(p_before)
gc()

min_features <- 200
max_features <- 6000
min_counts   <- 500
max_counts   <- 50000
mt_max       <- 10
rb_max       <- 40

seurat_obj$qc_pass <- (
  seurat_obj$nFeature_RNA >= min_features &
    seurat_obj$nFeature_RNA <= max_features &
    seurat_obj$nCount_RNA >= min_counts &
    seurat_obj$nCount_RNA <= max_counts &
    seurat_obj$percent.mt <= mt_max &
    seurat_obj$percent.rb <= rb_max
)
seurat_obj$qc_pass[is.na(seurat_obj$qc_pass)] <- FALSE

print(table(seurat_obj$sample_gsm, seurat_obj$qc_pass))

seurat_qc <- subset(seurat_obj, subset = qc_pass)
rm(seurat_obj)
gc()

rna_layers <- Layers(seurat_qc[["RNA"]])
if (length(rna_layers) > 1) seurat_qc <- JoinLayers(seurat_qc, assay = "RNA")

DefaultAssay(seurat_qc) <- "RNA"
sce <- as.SingleCellExperiment(seurat_qc)
sce <- scDblFinder(sce, samples = "sample_gsm", BPPARAM = SerialParam())

seurat_qc$scDblFinder.class <- sce$scDblFinder.class
seurat_qc$scDblFinder.score <- sce$scDblFinder.score

print(table(sce$sample_gsm, sce$scDblFinder.class))
rm(sce)
gc()

meta_df <- as.data.frame(seurat_qc[[]], stringsAsFactors = FALSE)
required_cols <- c("patient_id", "sample_gsm", "scDblFinder.class", "scDblFinder.score")
missing_cols <- setdiff(required_cols, colnames(meta_df))
if (length(missing_cols) > 0) {
  stop("seurat_qcmissing: ", paste(missing_cols, collapse = ", "))
}

if (is.list(meta_df$patient_id)) meta_df$patient_id <- unlist(meta_df$patient_id, use.names = FALSE)
if (is.list(meta_df$sample_gsm)) meta_df$sample_gsm <- unlist(meta_df$sample_gsm, use.names = FALSE)
if (is.list(meta_df$scDblFinder.class)) {
  meta_df$scDblFinder.class <- unlist(meta_df$scDblFinder.class, use.names = FALSE)
}

patient_levels <- sprintf("Tumor_%02d", 1:10)
meta_df$patient_id <- factor(as.character(meta_df$patient_id), levels = patient_levels)
meta_df$sample_gsm <- as.character(meta_df$sample_gsm)
meta_df$scDblFinder.class <- as.character(meta_df$scDblFinder.class)

doublet_df <- meta_df %>%
  dplyr::group_by(patient_id, scDblFinder.class) %>%
  dplyr::summarise(cell_number = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(patient_id) %>%
  dplyr::mutate(total_cells = sum(cell_number),
                percentage = 100 * cell_number / total_cells) %>%
  dplyr::ungroup()

p_doublet1 <- ggplot(doublet_df, aes(patient_id, percentage, fill = scDblFinder.class)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c(singlet = "#6fa6cf7f", doublet = "#b436658f")) +
  labs(title = "doubletfraction", x = "tumor_sample", y = "fraction(%)", fill = "class") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

p_doublet2 <- ggplot(meta_df, aes(patient_id, scDblFinder.score,
                                  fill = scDblFinder.class)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = c(singlet = "#6fa6cf7f", doublet = "#b436658f")) +
  labs(title = "doublet", x = "tumor_sample",
       y = "doublet", fill = "class") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

pdf(file.path(out_dir, "02_doubletplot.pdf"), width = 14, height = 6)
print(p_doublet1 + p_doublet2)
dev.off()

singlet_cells <- rownames(meta_df)[
  !is.na(meta_df$scDblFinder.class) & meta_df$scDblFinder.class == "singlet"
]
if (length(singlet_cells) == 0) stop("singletcell, checkscDblFinderresults")

seurat_singlet <- subset(seurat_qc, cells = singlet_cells)

print(table(seurat_singlet$patient_id))

rm(seurat_qc, singlet_cells, p_doublet1, p_doublet2)
gc()

p_after <- VlnPlot(
  seurat_singlet,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  group.by = "patient_id", pt.size = 0, ncol = 2
)

pdf(file.path(out_dir, "03_post_QCQCplot.pdf"), width = 14, height = 9)
print(p_after)
dev.off()

p_scatter1 <- FeatureScatter(
  seurat_singlet, feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA", group.by = "patient_id"
)
pdf(file.path(out_dir, "04_UMIgeneplot.pdf"), width = 8, height = 6)
print(p_scatter1)
dev.off()

p_scatter2 <- FeatureScatter(
  seurat_singlet, feature1 = "nCount_RNA",
  feature2 = "percent.mt", group.by = "patient_id"
)
pdf(file.path(out_dir, "05_UMIfractionplot.pdf"), width = 8, height = 6)
print(p_scatter2)
dev.off()

rm(p_after, p_scatter1, p_scatter2)
gc()

final_meta <- as.data.frame(seurat_singlet[[]], stringsAsFactors = FALSE)
final_cell_table <- table(factor(as.character(final_meta$patient_id),
                                 levels = patient_levels))

cell_summary <- data.frame(
  patient_id = patient_levels,
  sample_gsm = paste0("GSM", seq(7845552, 7845570, by = 2)),
  final_singlet_cells = as.integer(final_cell_table),
  stringsAsFactors = FALSE
)

write.csv(cell_summary, file.path(out_dir, "cellnumberstatistics.csv"), row.names = FALSE)
write.csv(doublet_df, file.path(out_dir, "doubletstatistics.csv"), row.names = FALSE)

print(cell_summary)

output_file <- file.path(out_dir, "GSE245601__QC_Singlet.rds")
saveRDS(seurat_singlet, output_file)

options(stringsAsFactors = FALSE)
set.seed(123)

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(celda)
library(SingleCellExperiment)
library(future)

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(base_dir, "01.QC", "GSE245601__QC_Singlet.rds")
out_dir <- file.path(base_dir, "02_Integration")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = 10000 * 1024^2)
plan(multisession, workers = 6)

if (!file.exists(input_file)) stop("not_foundinput: ", input_file)

seurat_obj <- readRDS(input_file)

required_meta <- c("sample_gsm", "patient_id", "percent.mt")
missing_meta <- setdiff(required_meta, colnames(seurat_obj[[]]))
if (length(missing_meta) > 0) {
  stop("missing: ", paste(missing_meta, collapse = ", "))
}

print(table(seurat_obj$patient_id))

plan(sequential)

rna_layers <- tryCatch(Layers(seurat_obj[["RNA"]]), error = function(e) character(0))
if (length(rna_layers) > 1) {
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
}

counts_raw <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

sce_decont <- SingleCellExperiment(
  assays = list(counts = counts_raw),
  colData = data.frame(
    sample_gsm = as.character(seurat_obj$sample_gsm),
    row.names = colnames(seurat_obj)
  )
)

sce_decont <- decontX(sce_decont, batch = sce_decont$sample_gsm)
seurat_obj$decontX_contamination <- colData(sce_decont)$decontX_contamination

print(quantile(seurat_obj$decontX_contamination, probs = seq(0, 1, 0.1)))

p_contam <- ggplot(
  data.frame(contamination = seurat_obj$decontX_contamination),
  aes(contamination)
) +
  geom_histogram(bins = 60, fill = "#6fa6cf8f", color = "white") +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "#b436659f") +
  labs(title = "decontX contamination", x = "Contamination proportion",
       y = "Cell count") +
  theme_bw(base_size = 12)

pdf(file.path(out_dir, "01_decontX.pdf"), width = 6, height = 4)
print(p_contam)
dev.off()

decontX_counts <- round(assay(sce_decont, "decontXcounts"))
seurat_obj <- SetAssayData(
  object = seurat_obj, assay = "RNA",
  layer = "counts", new.data = decontX_counts
)

n_before <- ncol(seurat_obj)
keep_cells <- colnames(seurat_obj)[seurat_obj$decontX_contamination < 0.2]
seurat_obj <- subset(seurat_obj, cells = keep_cells)
n_after <- ncol(seurat_obj)

rm(sce_decont, counts_raw, decontX_counts, keep_cells, p_contam)
gc()
plan(multisession, workers = 6)

seurat_obj <- NormalizeData(
  seurat_obj, normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
)

seurat_obj <- FindVariableFeatures(
  seurat_obj, selection.method = "vst",
  nfeatures = 2000, verbose = FALSE
)

top10 <- head(VariableFeatures(seurat_obj), 10)
p_hvg <- VariableFeaturePlot(seurat_obj)
p_hvg <- LabelPoints(p_hvg, points = top10, repel = TRUE)

pdf(file.path(out_dir, "02_geneplot.pdf"), width = 8, height = 5)
print(p_hvg)
dev.off()

seurat_obj <- CellCycleScoring(
  seurat_obj,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)
seurat_obj$Phase <- factor(seurat_obj$Phase, levels = c("G1", "S", "G2M"))

p_cc <- VlnPlot(
  seurat_obj, features = c("S.Score", "G2M.Score"),
  group.by = "patient_id", pt.size = 0, ncol = 2
)

pdf(file.path(out_dir, "03_cell.pdf"), width = 12, height = 5)
print(p_cc)
dev.off()

print(table(seurat_obj$Phase))

plan(sequential)

hvg_features <- VariableFeatures(seurat_obj)
seurat_obj <- ScaleData(
  seurat_obj,
  features = hvg_features,
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
  verbose = FALSE
)

plan(multisession, workers = 6)

seurat_obj <- RunPCA(
  seurat_obj, features = hvg_features,
  npcs = 50, verbose = FALSE
)

p_elbow <- ElbowPlot(seurat_obj, ndims = 50)
pdf(file.path(out_dir, "04_PCAplot.pdf"), width = 6, height = 4)
print(p_elbow)
dev.off()

pdf(file.path(out_dir, "05_PCAplot.pdf"), width = 14, height = 10)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

p_pca <- DimPlot(
  seurat_obj, reduction = "pca",
  group.by = "patient_id", pt.size = 0.1
) + ggtitle("PCA before Harmony") + theme_bw()

pdf(file.path(out_dir, "06_HarmonyPCAplot.pdf"), width = 9, height = 6)
print(p_pca)
dev.off()

seurat_obj <- RunHarmony(
  object = seurat_obj,
  group.by.vars = "sample_gsm",
  reduction = "pca",
  reduction.save = "harmony",
  plot_convergence = TRUE,
  max.iter.harmony = 30
)

p_harmony <- DimPlot(
  seurat_obj, reduction = "harmony",
  group.by = "patient_id", pt.size = 0.1
) + ggtitle("Harmony embedding") + theme_bw()

pdf(file.path(out_dir, "07_Harmonyresults.pdf"), width = 9, height = 6)
print(p_harmony)
dev.off()

output_file <- file.path(out_dir, "GSE245601_Step2_Harmony.rds")
saveRDS(seurat_obj, output_file)

print(table(seurat_obj$patient_id))

options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(SeuratObject)
library(SCP)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(future)
library(RColorBrewer)

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(base_dir, "02_Integration", "GSE245601_Step2_Harmony.rds")
out_dir <- file.path(base_dir, "03.UMAP")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = 10000 * 1024^2)
plan(multisession, workers = 6)

dims_use <- 1:30
selected_res <- 0.20
res_seq <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00)

if (!file.exists(input_file)) stop("not_foundinput: ", input_file)

seurat_obj <- readRDS(input_file)

required_reductions <- c("pca", "harmony")
missing_reductions <- setdiff(required_reductions, names(seurat_obj@reductions))
if (length(missing_reductions) > 0) {
  stop("missingresults: ", paste(missing_reductions, collapse = ", "))
}

required_meta <- c("sample_gsm", "patient_id", "percent.mt",
                   "nFeature_RNA", "nCount_RNA", "decontX_contamination")
missing_meta <- setdiff(required_meta, colnames(seurat_obj[[]]))
if (length(missing_meta) > 0) {
  stop("missing: ", paste(missing_meta, collapse = ", "))
}

print(table(seurat_obj$patient_id))

seurat_obj <- FindNeighbors(
  seurat_obj, reduction = "harmony",
  dims = dims_use, k.param = 20, verbose = FALSE
)

plan(sequential)

for (res in res_seq) {
  seurat_obj <- FindClusters(
    seurat_obj, resolution = res,
    algorithm = 1, random.seed = 42, verbose = FALSE
  )
  n_clust <- length(unique(seurat_obj$seurat_clusters))
}

pdf(file.path(out_dir, "01_clustering_resolution_tree.pdf"), width = 12, height = 9)
print(clustree(seurat_obj@meta.data, prefix = "RNA_snn_res."))
dev.off()

seurat_obj <- FindClusters(
  seurat_obj, resolution = selected_res,
  algorithm = 1, random.seed = 42, verbose = FALSE
)

seurat_obj$seurat_clusters <- factor(
  seurat_obj$seurat_clusters,
  levels = sort(unique(as.integer(as.character(seurat_obj$seurat_clusters))))
)
Idents(seurat_obj) <- "seurat_clusters"

cluster_tab <- table(seurat_obj$seurat_clusters)
print(cluster_tab)

res_cols <- grep("^RNA_snn_res\\.", colnames(seurat_obj[[]]), value = TRUE)
keep_col <- paste0("RNA_snn_res.", selected_res)
remove_cols <- setdiff(res_cols, keep_col)
for (nm in remove_cols) seurat_obj[[nm]] <- NULL

plan(multisession, workers = 6)

seurat_obj <- RunUMAP(
  seurat_obj, reduction = "harmony", dims =  1:20,
  n.neighbors = 30, min.dist = 0.30, spread = 1,
  seed.use = 42, reduction.name = "umap",
  reduction.key = "UMAP_", verbose = FALSE
)

seurat_obj <- RunUMAP(
  seurat_obj, reduction = "pca", dims =  1:20,
  n.neighbors = 30, min.dist = 0.30, spread = 1,
  seed.use = 42, reduction.name = "umap_pca",
  reduction.key = "pUMAP_", verbose = FALSE
)

n_clusters <- length(levels(seurat_obj$seurat_clusters))
n_samples <- length(unique(seurat_obj$patient_id))

cluster_palette <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"),
                     brewer.pal(8, "Dark2"))
clust_cols <- colorRampPalette(cluster_palette)(n_clusters)
names(clust_cols) <- levels(seurat_obj$seurat_clusters)

sample_palette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))
sample_cols <- colorRampPalette(sample_palette)(n_samples)
names(sample_cols) <- sort(unique(as.character(seurat_obj$patient_id)))

umap_df <- as.data.frame(Embeddings(seurat_obj, "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cluster <- seurat_obj$seurat_clusters

centers <- umap_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

theme_umap <- theme_scp() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    aspect.ratio = 1
  )

p_cluster1 <- CellDimPlot(
  seurat_obj, group.by = "seurat_clusters", reduction = "umap",
  cols = clust_cols, pt.size = 0.01, alpha = 0.75, label = FALSE
) + ggtitle("UMAP clusters") + theme_umap

p_cluster2 <- p_cluster1 +
  geom_text(
    data = centers, aes(UMAP_1, UMAP_2, label = cluster),
    color = "black", size = 3.5, fontface = "bold",
    inherit.aes = FALSE
  ) + ggtitle("UMAP clusters with labels")

pdf(file.path(out_dir, "02_UMAPplot.pdf"), width = 13, height = 6)
print((p_cluster1 | p_cluster2) + plot_layout(guides = "collect") &
        theme(legend.position = "right"))
dev.off()

p_sample <- CellDimPlot(
  seurat_obj, group.by = "patient_id", reduction = "umap",
  cols = sample_cols, pt.size = 0.01, alpha = 0.65, label = FALSE
) + ggtitle("UMAP by sample") + theme_umap

pdf(file.path(out_dir, "03_sampleUMAPplot.pdf"), width = 9, height = 7)
print(p_sample)
dev.off()

p_split <- DimPlot(
  seurat_obj, reduction = "umap", group.by = "seurat_clusters",
  split.by = "patient_id", cols = clust_cols,
  pt.size = 0.01, ncol = 5
) + NoLegend()

pdf(file.path(out_dir, "04_sampleUMAPplot.pdf"), width = 15, height = 8)
print(p_split)
dev.off()

qc_features <- c("percent.mt", "nFeature_RNA",
                 "nCount_RNA", "decontX_contamination")

p_qc <- FeaturePlot(
  seurat_obj, features = qc_features, reduction = "umap",
  pt.size = 0.01, ncol = 2, order = TRUE
) & theme(aspect.ratio = 1)

pdf(file.path(out_dir, "05_QCUMAPplot.pdf"), width = 11, height = 10)
print(p_qc)
dev.off()

p_pre <- DimPlot(
  seurat_obj, reduction = "umap_pca",
  group.by = "patient_id", cols = sample_cols, pt.size = 0.01
) + ggtitle("PCA-based UMAP") + theme_umap + NoLegend()

p_post <- DimPlot(
  seurat_obj, reduction = "umap",
  group.by = "patient_id", cols = sample_cols, pt.size = 0.01
) + ggtitle("Harmony-based UMAP") + theme_umap

pdf(file.path(out_dir, "06_PCAHarmonyUMAP.pdf"), width = 13, height = 6)
print(p_pre | p_post)
dev.off()

meta_df <- as.data.frame(seurat_obj[[]], stringsAsFactors = FALSE)

cluster_summary <- meta_df %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_samples = dplyr::n_distinct(patient_id),
    .groups = "drop"
  ) %>%
  dplyr::arrange(as.integer(as.character(seurat_clusters)))

cluster_sample <- meta_df %>%
  dplyr::group_by(seurat_clusters, patient_id) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(pct_cluster = 100 * n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

p_comp <- ggplot(
  cluster_sample,
  aes(x = seurat_clusters, y = pct_cluster, fill = patient_id)
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = sample_cols) +
  labs(x = "Cluster", y = "Cell proportion (%)", fill = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(out_dir, "07_sample.pdf"), width = 11, height = 6)
print(p_comp)
dev.off()

write.csv(cluster_summary, file.path(out_dir, "cellstatistics.csv"),
          row.names = FALSE)
write.csv(cluster_sample, file.path(out_dir, "sample.csv"),
          row.names = FALSE)

output_file <- file.path(out_dir, "GSE245601_Step3_Clustered_UMAP.rds")
saveRDS(seurat_obj, output_file)

options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(SeuratObject)
library(SingleR)
library(celldex)
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(base_dir, "03.UMAP", "GSE245601_Step3_Clustered_UMAP.rds")
out_dir <- file.path(base_dir, "05.singleR")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("not_foundinput: ", input_file)

seurat_obj <- readRDS(input_file)

if (!"RNA" %in% names(seurat_obj@assays)) stop("RNA assay")
if (!"seurat_clusters" %in% colnames(seurat_obj[[]])) {
  stop("seurat_clusters")
}
if (!"umap" %in% names(seurat_obj@reductions)) stop("UMAPresults")

DefaultAssay(seurat_obj) <- "RNA"

rna_layers <- Layers(seurat_obj[["RNA"]])
if (length(rna_layers) > 1) {
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
}

if (!"data" %in% Layers(seurat_obj[["RNA"]])) {
  seurat_obj <- NormalizeData(
    seurat_obj, normalization.method = "LogNormalize",
    scale.factor = 10000, verbose = FALSE
  )
}

cluster_num <- sort(unique(as.integer(as.character(seurat_obj$seurat_clusters))))
seurat_obj$seurat_clusters <- factor(
  seurat_obj$seurat_clusters,
  levels = cluster_num
)
Idents(seurat_obj) <- "seurat_clusters"

expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
cluster_ids <- seurat_obj$seurat_clusters

print(table(cluster_ids))

ref_bpe <- celldex::BlueprintEncodeData(ensembl = FALSE)

singler_cell <- SingleR(
  test = expr_matrix,
  ref = ref_bpe,
  labels = ref_bpe$label.main,
  assay.type.ref = "logcounts",
  BPPARAM = SerialParam()
)

cell_match <- match(colnames(seurat_obj), rownames(singler_cell))
if (any(is.na(cell_match))) stop("cellSingleRresults")

raw_label <- as.character(singler_cell$labels[cell_match])

if ("pruned.labels" %in% colnames(singler_cell)) {
  pruned_label <- as.character(singler_cell$pruned.labels[cell_match])
} else {
  pruned_label <- rep(NA_character_, length(raw_label))
}

final_cell_label <- ifelse(
  is.na(pruned_label) | pruned_label == "",
  raw_label,
  pruned_label
)

seurat_obj$SingleR_raw_label <- setNames(raw_label, colnames(seurat_obj))
seurat_obj$SingleR_pruned_label <- setNames(pruned_label, colnames(seurat_obj))
seurat_obj$SingleR_cell_label <- setNames(final_cell_label, colnames(seurat_obj))
seurat_obj$SingleR_pruned_status <- setNames(
  ifelse(is.na(pruned_label) | pruned_label == "", "Not_pruned", "Pruned"),
  colnames(seurat_obj)
)

print(sort(table(seurat_obj$SingleR_cell_label), decreasing = TRUE))

vote_table <- data.frame(
  Cluster = as.character(cluster_ids),
  CellLabel = as.character(seurat_obj$SingleR_cell_label),
  stringsAsFactors = FALSE
) %>%
  dplyr::group_by(Cluster, CellLabel) %>%
  dplyr::summarise(CellNumber = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(
    ClusterTotal = sum(CellNumber),
    Percentage = 100 * CellNumber / ClusterTotal
  ) %>%
  dplyr::arrange(dplyr::desc(CellNumber), .by_group = TRUE) %>%
  dplyr::mutate(Rank = dplyr::row_number()) %>%
  dplyr::ungroup()

top1 <- vote_table %>%
  dplyr::filter(Rank == 1) %>%
  dplyr::transmute(
    Cluster,
    SingleR_label = CellLabel,
    Dominant_cells = CellNumber,
    Cluster_cells = ClusterTotal,
    Dominant_pct = round(Percentage, 2)
  )

top2 <- vote_table %>%
  dplyr::filter(Rank == 2) %>%
  dplyr::transmute(
    Cluster,
    Second_label = CellLabel,
    Second_pct = round(Percentage, 2)
  )

annotation_table <- top1 %>%
  dplyr::left_join(top2, by = "Cluster") %>%
  dplyr::mutate(
    Vote_margin = round(Dominant_pct - ifelse(is.na(Second_pct), 0, Second_pct), 2)
  ) %>%
  dplyr::arrange(as.integer(Cluster))

print(annotation_table)

cluster_to_type <- setNames(
  annotation_table$SingleR_label,
  annotation_table$Cluster
)

seurat_obj$SingleR_cluster_label <- setNames(
  cluster_to_type[as.character(seurat_obj$seurat_clusters)],
  colnames(seurat_obj)
)

if (any(is.na(seurat_obj$SingleR_cluster_label))) {
  stop("cluster")
}

seurat_obj$SingleR_cluster_label <- factor(
  seurat_obj$SingleR_cluster_label,
  levels = sort(unique(seurat_obj$SingleR_cluster_label))
)

write.csv(
  annotation_table,
  file.path(out_dir, "01_SingleRannotation.csv"),
  row.names = FALSE
)

write.csv(
  vote_table,
  file.path(out_dir, "02_SingleR.csv"),
  row.names = FALSE
)

celltype_summary <- as.data.frame(table(seurat_obj$SingleR_cluster_label))
colnames(celltype_summary) <- c("CellType", "CellNumber")

write.csv(
  celltype_summary,
  file.path(out_dir, "03_SingleRcellnumber.csv"),
  row.names = FALSE
)

singler_cluster <- SingleR(
  test = expr_matrix,
  ref = ref_bpe,
  labels = ref_bpe$label.main,
  clusters = cluster_ids,
  assay.type.ref = "logcounts",
  BPPARAM = SerialParam()
)

pdf(file.path(out_dir, "04_SingleRplot.pdf"), width = 10, height = 8)
plotScoreHeatmap(
  singler_cluster,
  show.pruned = TRUE,
  cluster_cols = TRUE,
  max.labels = 50
)
dev.off()

celltypes <- levels(seurat_obj$SingleR_cluster_label)
n_types <- length(celltypes)

base_cols <- c(
  brewer.pal(12, "Set3"),
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Dark2")
)
anno_cols <- colorRampPalette(base_cols)(n_types)
names(anno_cols) <- celltypes

p_cluster <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "SingleR_cluster_label",
  cols = anno_cols,
  pt.size = 0.08,
  label = TRUE,
  repel = TRUE,
  label.size = 3.5
) +
  ggtitle("SingleR cluster annotation") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    aspect.ratio = 1
  )

pdf(file.path(out_dir, "05_SingleRannotationUMAP.pdf"), width = 10, height = 8)
print(p_cluster)
dev.off()

cell_labels <- sort(unique(seurat_obj$SingleR_cell_label))
cell_cols <- colorRampPalette(base_cols)(length(cell_labels))
names(cell_cols) <- cell_labels

p_cell <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "SingleR_cell_label",
  cols = cell_cols,
  pt.size = 0.05
) +
  ggtitle("SingleR cell-level labels") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    aspect.ratio = 1
  )

pdf(file.path(out_dir, "06_SingleRcellUMAP.pdf"), width = 10, height = 8)
print(p_cell)
dev.off()

p_vote <- ggplot(
  vote_table,
  aes(x = factor(Cluster, levels = as.character(cluster_num)),
      y = Percentage, fill = CellLabel)
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cell_cols) +
  labs(
    title = "SingleR labels within each cluster",
    x = "Cluster", y = "Cell proportion (%)", fill = "SingleR label"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 7)
  )

pdf(file.path(out_dir, "07_SingleR.pdf"), width = 12, height = 7)
print(p_vote)
dev.off()

output_file <- file.path(out_dir, "GSE245601_Step4_SingleR_Annotated.rds")
saveRDS(seurat_obj, output_file)

save(
  singler_cell,
  singler_cluster,
  annotation_table,
  vote_table,
  file = file.path(out_dir, "SingleRanalysisresults.RData")
)

options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(
  base_dir, "05.singleR",
  "GSE245601_Step4_SingleR_Annotated.rds"
)
out_dir <- file.path(base_dir, "06.marker")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("not_foundinput: ", input_file)

seurat_obj <- readRDS(input_file)

if (!"RNA" %in% names(seurat_obj@assays)) stop("RNA assay")
if (!"seurat_clusters" %in% colnames(seurat_obj[[]])) {
  stop("seurat_clusters")
}
if (!"umap" %in% names(seurat_obj@reductions)) {
  stop("UMAPresults")
}

DefaultAssay(seurat_obj) <- "RNA"

rna_layers <- Layers(seurat_obj[["RNA"]])
if (length(rna_layers) > 1) {
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
}

if (!"data" %in% Layers(seurat_obj[["RNA"]])) {
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

Idents(seurat_obj) <- "seurat_clusters"

print(table(seurat_obj$seurat_clusters))

marker_panel <- list(
  Epi = c("EPCAM", "KRT8", "KRT18", "KRT19", "FOXA1"),

  CD4T = c("IL7R", "CCR7", "TCF7", "LTB"),

  CD8T = c("CD8A", "CD8B", "CCL5", "GZMK"),

  B = c("MS4A1", "CD79A", "CD74", "CD37"),

  Macro = c("LST1", "FCER1G", "CD68", "C1QC"),

  DC = c("CD1C", "FCER1A", "CLEC10A", "CLEC9A"),

  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT"),

  Fib = c("COL1A1", "COL1A2", "DCN", "LUM"),

  Endo = c("PECAM1", "VWF", "CLDN5", "CDH5")
)

feature_info <- do.call(
  rbind,
  lapply(names(marker_panel), function(cell_type) {
    data.frame(
      gene = marker_panel[[cell_type]],
      cell_type = cell_type,
      stringsAsFactors = FALSE
    )
  })
)

feature_info <- feature_info[!duplicated(feature_info$gene), ]

feature_info <- feature_info[
  feature_info$gene %in% rownames(seurat_obj),
]

if (nrow(feature_info) == 0) {
  stop("Markergene, checkgene")
}

feature_info$display <- paste0(
  feature_info$gene,
  " [",
  feature_info$cell_type,
  "]"
)

features_use <- feature_info$gene
display_labels <- setNames(
  feature_info$display,
  feature_info$gene
)

p_dot <- DotPlot(
  seurat_obj,
  features = features_use,
  assay = "RNA",
  group.by = "seurat_clusters",
  cols = c("grey90", "#b43665"),
  dot.scale = 6
) +
  scale_x_discrete(
    limits = features_use,
    labels = display_labels
  ) +
  RotatedAxis() +
  labs(
    title = "Canonical markers across Seurat clusters",
    x = NULL,
    y = "Cluster",
    color = "Average expression",
    size = "Percent expressed"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1,
      face = "italic"
    ),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_line(
      color = "grey92",
      linewidth = 0.25
    )
  )

pdf(
  file.path(out_dir, "01_cellMarkerplot.pdf"),
  width = 20,
  height = 8
)
print(p_dot)
dev.off()

feature_plots <- lapply(seq_len(nrow(feature_info)), function(i) {
  gene_i <- feature_info$gene[i]
  title_i <- feature_info$display[i]

  FeaturePlot(
    seurat_obj,
    features = gene_i,
    reduction = "umap",
    pt.size = 0.01,
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95"
  ) +
    ggtitle(title_i) +
    theme_classic() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "italic",
        size = 10
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      aspect.ratio = 1
    )
})

ncol_plot <- 5
nrow_plot <- ceiling(length(feature_plots) / ncol_plot)

p_umap <- wrap_plots(
  feature_plots,
  ncol = ncol_plot
)

pdf(
  file.path(out_dir, "02_cellMarker_UMAPplot.pdf"),
  width = ncol_plot * 3.6,
  height = nrow_plot * 3.3
)
print(p_umap)
dev.off()

options(stringsAsFactors = FALSE)
set.seed(42)

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(
  base_dir, "05.singleR",
  "GSE245601_Step4_SingleR_Annotated.rds"
)
out_dir <- file.path(base_dir, "07.celltype_UMAP")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("not_foundinput: ", input_file)

seurat_obj <- readRDS(input_file)

if (!"seurat_clusters" %in% colnames(seurat_obj[[]])) {
  stop("seurat_clusters")
}
if (!"umap" %in% names(seurat_obj@reductions)) {
  stop("umapresults")
}

cluster_ids <- as.character(sort(unique(
  as.integer(as.character(seurat_obj$seurat_clusters))
)))

cluster_to_type <- c(
  "0"  = "T cells",
  "1"  = "Epithelial cells",
  "2"  = "Fibroblasts",
  "3"  = "Epithelial cells",
  "4"  = "Macrophages",
  "5"  = "B cells",
  "6"  = "Mast cells",
  "7"  = "Epithelial cells",
  "8"  = "Endothelial cells",
  "9"  = "Fibroblasts",
  "10" = "Epithelial cells",
  "11" = "B cells",
  "12" = "Epithelial cells",
  "13" = "DCs"
)

missing_cluster <- setdiff(cluster_ids, names(cluster_to_type))
extra_cluster <- setdiff(names(cluster_to_type), cluster_ids)

if (length(missing_cluster) > 0) {
  stop("missingcluster: ", paste(missing_cluster, collapse = ", "))
}
if (length(extra_cluster) > 0) {
  stop("cluster: ",
       paste(extra_cluster, collapse = ", "))
}

blank_cluster <- names(cluster_to_type)[trimws(cluster_to_type) == ""]

if (length(blank_cluster) > 0) {
  stop(
    "clustercell: ",
    paste(blank_cluster, collapse = ", ")
  )
}

celltype_vec <- unname(
  cluster_to_type[as.character(seurat_obj$seurat_clusters)]
)

celltype_levels <- unique(unname(cluster_to_type[cluster_ids]))

seurat_obj$manual_celltype <- factor(
  celltype_vec,
  levels = celltype_levels
)

Idents(seurat_obj) <- "manual_celltype"

print(table(seurat_obj$manual_celltype))

default_colors <- c(
  "Epithelial cells" = "#D95F5F",
  "T cells" = "#3C8DAD",
  "B cells" = "#7E57C2",
  "Macrophages" = "#D98C3F",
  "DCs" = "#9C6B30",
  "Mast cells" = "#E6B84A",
  "Fibroblasts" = "#66A86C",
  "Endothelial cells" = "#C85A9E"
)

celltypes <- levels(seurat_obj$manual_celltype)
celltype_cols <- default_colors[celltypes]

missing_color <- is.na(celltype_cols)
if (any(missing_color)) {
  celltype_cols[missing_color] <- grDevices::hcl.colors(
    sum(missing_color), palette = "Dark 3"
  )
}
names(celltype_cols) <- celltypes

cluster_cols <- setNames(
  grDevices::hcl.colors(length(cluster_ids), palette = "Dark 3"),
  cluster_ids
)

p_cluster <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = cluster_cols,
  pt.size = 0.08,
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = FALSE
) +
  ggtitle("Seurat clusters") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1
  )

p_celltype <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "manual_celltype",
  cols = celltype_cols,
  pt.size = 0.08,
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  raster = FALSE
) +
  ggtitle("Manual cell-type annotation") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    aspect.ratio = 1
  )

pdf(
  file.path(out_dir, "01_cellUMAP.pdf"),
  width = 14, height = 6
)
print(p_cluster | p_celltype)
dev.off()

pdf(
  file.path(out_dir, "02_cellUMAP.pdf"),
  width = 8, height = 7
)
print(p_celltype)
dev.off()

mapping_table <- data.frame(
  Cluster = cluster_ids,
  CellType = unname(cluster_to_type[cluster_ids]),
  stringsAsFactors = FALSE
)

celltype_count <- as.data.frame(table(seurat_obj$manual_celltype))
colnames(celltype_count) <- c("CellType", "CellNumber")

write.csv(
  mapping_table,
  file.path(out_dir, "00_Clustercell.csv"),
  row.names = FALSE
)

write.csv(
  celltype_count,
  file.path(out_dir, "03_cellnumber.csv"),
  row.names = FALSE
)

output_file <- file.path(
  out_dir,
  "GSE245601_Step5_Manual_Annotated.rds"
)
saveRDS(seurat_obj, output_file)

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(SCIPAC)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_rds <- file.path(base_dir, "07.celltype_UMAP",
                       "GSE245601_Step5_Manual_Annotated.rds")
out_dir <- file.path(base_dir, "08.tumorselect")

tcga_dir <- file.path(PROJECT_DIR, "data/raw/TCGA_BRCA")
clinical_file <- file.path(tcga_dir, "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix")
tpm_file <- file.path(tcga_dir, "tpm", "tpm_BRCA_mRNA_symbolgenenames.txt")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c(input_rds, clinical_file, tpm_file)) {
  if (!file.exists(f)) stop("not_found: ", f)
}

EPI_HVG <- 2000L
EPI_PC <- 20L
EPI_RES <- 0.10
SCIPAC_HVG <- 1000L
SCIPAC_PC <- 30L
SCIPAC_RES <- 1
NUM_CORES <- if (.Platform$OS.type == "windows") 1L else {
  max(1L, parallel::detectCores() - 1L)
}

find_col <- function(df, candidates, label) {
  idx <- match(tolower(candidates), tolower(colnames(df)))
  idx <- idx[!is.na(idx)]
  if (!length(idx)) stop("not_found", label, ": ",
                         paste(candidates, collapse = ", "))
  colnames(df)[idx[1]]
}

patient_id <- function(x) {
  substr(gsub("\\.", "-", as.character(x)), 1, 12)
}

sample_code <- function(x) {
  substr(gsub("\\.", "-", as.character(x)), 14, 15)
}

clinical <- read.delim(
  clinical_file, sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, check.names = FALSE
)

sample_col <- find_col(clinical,
                       c("sampleID", "sample_id", "SampleID"), "sample")
er_col <- find_col(clinical,
                   c("ER_Status_nature2012", "ER_Status",
                     "ER status", "er_status_by_ihc"), "ER")
her2_col <- find_col(clinical,
                     c("HER2_Final_Status_nature2012",
                       "HER2_Status_nature2012", "HER2_Status",
                       "HER2 status", "her2_status_by_ihc"), "HER2")

clinical$patient_key <- patient_id(clinical[[sample_col]])
clinical$ER_clean <- toupper(trimws(as.character(clinical[[er_col]])))
clinical$HER2_clean <- toupper(trimws(as.character(clinical[[her2_col]])))

selected_patients <- unique(clinical$patient_key[
  grepl("^POS", clinical$ER_clean) &
    grepl("^NEG", clinical$HER2_clean)
])
if (!length(selected_patients)) {
  stop("ER+/HER2-, checkER/HER2")
}

tpm <- as.matrix(read.table(
  tpm_file, sep = "\t", header = TRUE, row.names = 1,
  stringsAsFactors = FALSE, check.names = FALSE
))
storage.mode(tpm) <- "numeric"
if (anyDuplicated(rownames(tpm))) tpm <- rowsum(tpm, rownames(tpm))

tumor_samples <- colnames(tpm)[
  sample_code(colnames(tpm)) == "01" &
    patient_id(colnames(tpm)) %in% selected_patients
]
normal_samples <- colnames(tpm)[sample_code(colnames(tpm)) == "11"]

if (!length(tumor_samples)) stop("TPMER+/HER2-")
if (!length(normal_samples)) stop("TPMsample")

bulk_samples <- c(tumor_samples, normal_samples)
bulk_mat <- tpm[, bulk_samples, drop = FALSE]
bulk_group <- setNames(
  c(rep("Tumor", length(tumor_samples)),
    rep("Normal", length(normal_samples))),
  bulk_samples
)

write.csv(
  bulk_mat,
  file.path(out_dir, "00_TCGA_ERpos_HER2neg_Tumor_Normal_TPM.csv"),
  quote = FALSE
)
write.csv(
  data.frame(
    Sample = bulk_samples,
    Patient = patient_id(bulk_samples),
    Group = unname(bulk_group[bulk_samples])
  ),
  file.path(out_dir, "00_TCGAsample.csv"),
  row.names = FALSE
)

full_obj <- readRDS(input_rds)
if (!"manual_celltype" %in% colnames(full_obj[[]])) {
  stop("manual_celltype")
}

epi_old <- subset(full_obj, subset = manual_celltype == "Epithelial cells")
rm(full_obj); gc()

DefaultAssay(epi_old) <- "RNA"
if (length(Layers(epi_old[["RNA"]])) > 1) {
  epi_old <- JoinLayers(epi_old, assay = "RNA")
}

epi_obj <- CreateSeuratObject(
  counts = GetAssayData(epi_old, assay = "RNA", layer = "counts"),
  meta.data = epi_old[[]],
  project = "GSE245601_Epithelial"
)
rm(epi_old); gc()

batch_col <- if ("sample_gsm" %in% colnames(epi_obj[[]])) {
  "sample_gsm"
} else {
  "orig.ident"
}
if (!batch_col %in% colnames(epi_obj[[]])) stop("missingsample")

epi_obj <- NormalizeData(epi_obj, verbose = FALSE)
epi_obj <- FindVariableFeatures(epi_obj, nfeatures = EPI_HVG, verbose = FALSE)

regress_vars <- intersect(c("percent.mt", "S.Score", "G2M.Score"),
                          colnames(epi_obj[[]]))
epi_obj <- ScaleData(epi_obj, vars.to.regress = regress_vars, verbose = FALSE)
epi_obj <- RunPCA(epi_obj, npcs = EPI_PC, verbose = FALSE)
epi_obj <- epi_obj %>%
  RunHarmony(
    group.by.vars = batch_col,
    verbose = FALSE
  )
epi_obj <- FindNeighbors(epi_obj, reduction = "harmony",
                         dims = 1:EPI_PC, verbose = FALSE)
epi_obj <- FindClusters(epi_obj, resolution = EPI_RES, verbose = FALSE)
epi_obj <- RunUMAP(epi_obj, reduction = "harmony",
                   dims = 1:EPI_PC, verbose = FALSE)

p_cluster <- DimPlot(
  epi_obj, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE, pt.size = 0.05
) +
  ggtitle("Reclustered epithelial cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1)

ggsave(file.path(out_dir, "01_cellUMAP.pdf"),
       p_cluster, width = 7, height = 6)
saveRDS(epi_obj, file.path(out_dir, "01_Epithelial_Reclustered.rds"))

sc_counts <- GetAssayData(epi_obj, assay = "RNA", layer = "counts")
common_genes <- intersect(rownames(sc_counts), rownames(bulk_mat))
if (length(common_genes) < 500) stop("cellTCGAgene500")

tmp <- CreateSeuratObject(sc_counts[common_genes, , drop = FALSE],
                          min.cells = 3, min.features = 100)
tmp <- NormalizeData(tmp, verbose = FALSE)
tmp <- FindVariableFeatures(tmp, nfeatures = SCIPAC_HVG, verbose = FALSE)

hvg <- intersect(VariableFeatures(tmp), rownames(bulk_mat))
if (length(hvg) < 200) stop("SCIPACHVG200")

sc_prep <- as.matrix(
  GetAssayData(tmp, assay = "RNA", layer = "data")[hvg, , drop = FALSE]
)
bulk_prep <- log1p(as.matrix(bulk_mat[hvg, , drop = FALSE]))
rm(tmp, sc_counts); gc()

pca_res <- SCIPAC::sc.bulk.pca(
  sc_prep, bulk_prep, do.pca.sc = FALSE, n.pc = SCIPAC_PC
)
sc_rot <- pca_res$sc.dat.rot
bulk_rot <- pca_res$bulk.dat.rot
ct_res <- SCIPAC::seurat.ct(sc_rot, res = SCIPAC_RES)

y_bin <- factor(
  bulk_group[rownames(bulk_rot)],
  levels = c("Normal", "Tumor")
)

scipac_res <- SCIPAC::SCIPAC(
  bulk.dat = bulk_rot,
  y = y_bin,
  family = "binomial",
  ct.res = ct_res,
  ela.net.alpha = 0.4,
  bt.size = 50L,
  numCores = NUM_CORES,
  CI.alpha = 0.05,
  nfold = 10L
)

write.csv(scipac_res, file.path(out_dir, "02_SCIPACcellresults.csv"))

epi_obj$scipac_lambda <- NA_real_
epi_obj$scipac_sig <- NA_character_
common_cells <- intersect(colnames(epi_obj), rownames(scipac_res))

epi_obj@meta.data[common_cells, "scipac_lambda"] <-
  scipac_res[common_cells, "Lambda.est"]
epi_obj@meta.data[common_cells, "scipac_sig"] <-
  as.character(scipac_res[common_cells, "sig"])

epi_obj$scipac_seed <- case_when(
  epi_obj$scipac_sig == "Sig.pos" & epi_obj$scipac_lambda > 0 ~
    "High_conf_Tumor",
  epi_obj$scipac_sig == "Sig.neg" & epi_obj$scipac_lambda < 0 ~
    "High_conf_Normal",
  TRUE ~ "Uncertain"
)

p_lambda <- FeaturePlot(
  epi_obj, features = "scipac_lambda",
  reduction = "umap", pt.size = 0.02, order = TRUE
) +
  scale_color_gradient2(low = "#6fa6cf", mid = "grey85",
                        high = "#b43665", midpoint = 0) +
  ggtitle("SCIPAC lambda: ER+/HER2- tumor vs normal")

p_seed <- DimPlot(
  epi_obj, reduction = "umap", group.by = "scipac_seed",
  cols = c("High_conf_Tumor" = "#b43665",
           "High_conf_Normal" = "#6fa6cf",
           "Uncertain" = "grey85"),
  pt.size = 0.05
) +
  ggtitle("SCIPAC high-confidence seeds")

ggsave(file.path(out_dir, "02_SCIPACresultsUMAP.pdf"),
       p_lambda | p_seed, width = 13, height = 6)

epi_obj@misc$SCIPAC <- list(
  selected_group = "TCGA ER+/HER2- primary tumor vs normal breast",
  batch_col = batch_col,
  ct_res = ct_res,
  parameters = list(HVG = SCIPAC_HVG, PC = SCIPAC_PC,
                    resolution = SCIPAC_RES)
)

saveRDS(epi_obj, file.path(out_dir, "02_Epithelial_SCIPAC.rds"))

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
out_dir <- file.path(base_dir, "08.tumorselect")
input_rds <- file.path(out_dir, "02_Epithelial_SCIPAC.rds")

if (!file.exists(input_rds)) stop("not_found: ", input_rds)

TOP_N <- 50L
MAX_ITER <- 20L
SIG_STABLE <- 0.95
LABEL_STABLE <- 0.99

get_signature <- function(seu, tumor_cells, normal_cells) {
  if (length(tumor_cells) < 10 || length(normal_cells) < 10) {
    stop("TumorNormal10")
  }

  seu$tmp_group <- NA_character_
  seu$tmp_group[colnames(seu) %in% tumor_cells] <- "Tumor"
  seu$tmp_group[colnames(seu) %in% normal_cells] <- "Normal"
  Idents(seu) <- "tmp_group"

  deg <- FindMarkers(
    seu, ident.1 = "Tumor", ident.2 = "Normal",
    test.use = "wilcox", logfc.threshold = 0.1,
    min.pct = 0.1, only.pos = FALSE, verbose = FALSE
  )
  deg <- deg[deg$p_val_adj < 0.05, , drop = FALSE]
  if (nrow(deg) < 10) stop("gene10")

  up <- deg[deg$avg_log2FC > 0, , drop = FALSE]
  down <- deg[deg$avg_log2FC < 0, , drop = FALSE]

  tumor_sig <- head(
    rownames(up[order(up$avg_log2FC, decreasing = TRUE), , drop = FALSE]),
    TOP_N
  )
  normal_sig <- head(
    rownames(down[order(down$avg_log2FC), , drop = FALSE]),
    TOP_N
  )

  list(tumor = tumor_sig, normal = normal_sig, deg = deg)
}

score_classify <- function(seu, sig, suffix) {
  tumor_sig <- intersect(sig$tumor, rownames(seu))
  normal_sig <- intersect(sig$normal, rownames(seu))
  if (length(tumor_sig) < 5 || length(normal_sig) < 5) {
    stop("TumorNormal signaturegene5")
  }

  seu <- AddModuleScore(
    seu, features = list(tumor_sig),
    name = paste0("tumor_score_", suffix)
  )
  seu <- AddModuleScore(
    seu, features = list(normal_sig),
    name = paste0("normal_score_", suffix)
  )

  t_col <- paste0("tumor_score_", suffix, "1")
  n_col <- paste0("normal_score_", suffix, "1")
  d_col <- paste0("diff_score_", suffix)

  seu[[d_col]] <- as.numeric(scale(
    seu[[t_col, drop = TRUE]] - seu[[n_col, drop = TRUE]]
  ))

  km <- kmeans(seu[[d_col, drop = TRUE]], centers = 2, nstart = 50)
  centers <- tapply(seu[[d_col, drop = TRUE]], km$cluster, mean)
  tumor_k <- as.integer(names(which.max(centers)))

  cls <- ifelse(km$cluster == tumor_k, "Tumor", "Normal")
  names(cls) <- colnames(seu)

  list(seu = seu, class = cls,
       tumor_col = t_col, normal_col = n_col, diff_col = d_col)
}

sig_stability <- function(new_sig, old_sig) {
  mean(c(
    length(intersect(new_sig$tumor, old_sig$tumor)) /
      max(length(old_sig$tumor), 1),
    length(intersect(new_sig$normal, old_sig$normal)) /
      max(length(old_sig$normal), 1)
  ))
}

epi_obj <- readRDS(input_rds)
DefaultAssay(epi_obj) <- "RNA"

if (!"scipac_seed" %in% colnames(epi_obj[[]])) {
  stop("scipac_seed")
}

seed_tumor <- colnames(epi_obj)[
  epi_obj$scipac_seed == "High_conf_Tumor"
]
seed_normal <- colnames(epi_obj)[
  epi_obj$scipac_seed == "High_conf_Normal"
]

current_sig <- get_signature(epi_obj, seed_tumor, seed_normal)
current_fit <- score_classify(epi_obj, current_sig, "01")
epi_obj <- current_fit$seu
current_class <- current_fit$class
final_iter <- "01"

stability_log <- data.frame(
  Iteration = character(),
  Signature_stability = numeric(),
  Label_stability = numeric(),
  Tumor = integer(),
  Normal = integer()
)

for (i in seq_len(MAX_ITER)) {
  iter <- sprintf("%02d", i + 1)

  new_sig <- get_signature(
    epi_obj,
    names(current_class)[current_class == "Tumor"],
    names(current_class)[current_class == "Normal"]
  )
  new_fit <- score_classify(epi_obj, new_sig, iter)
  epi_obj <- new_fit$seu
  new_class <- new_fit$class

  s_stab <- sig_stability(new_sig, current_sig)
  l_stab <- mean(new_class[colnames(epi_obj)] ==
                   current_class[colnames(epi_obj)])

  stability_log <- rbind(
    stability_log,
    data.frame(
      Iteration = iter,
      Signature_stability = round(s_stab, 4),
      Label_stability = round(l_stab, 4),
      Tumor = sum(new_class == "Tumor"),
      Normal = sum(new_class == "Normal")
    )
  )

  current_sig <- new_sig
  current_fit <- new_fit
  current_class <- new_class
  final_iter <- iter

  if (s_stab >= SIG_STABLE && l_stab >= LABEL_STABLE) break
}

epi_obj$tumorselect_label <- current_class[colnames(epi_obj)]
epi_obj$tumorselect_diff_score <- epi_obj[[
  current_fit$diff_col, drop = TRUE
]]

write.csv(
  stability_log,
  file.path(out_dir, "03_phenotype_iteration_stability.csv"),
  row.names = FALSE
)

write.csv(
  rbind(
    data.frame(Type = "Tumor", Gene = current_sig$tumor),
    data.frame(Type = "Normal", Gene = current_sig$normal)
  ),
  file.path(out_dir, "03_finalTumor_Normal_Signature.csv"),
  row.names = FALSE
)

p_class <- DimPlot(
  epi_obj, reduction = "umap", group.by = "tumorselect_label",
  cols = c("Tumor" = "#b43665", "Normal" = "#6fa6cf"),
  pt.size = 0.05
) +
  ggtitle(paste0("Tumor vs Normal epithelial cells, iter ", final_iter)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        aspect.ratio = 1)

p_score <- FeaturePlot(
  epi_obj, features = "tumorselect_diff_score",
  reduction = "umap", pt.size = 0.02, order = TRUE
) +
  scale_color_gradient2(low = "#6fa6cf", mid = "grey85",
                        high = "#b43665", midpoint = 0) +
  ggtitle("Tumor - Normal phenotype score")

ggsave(
  file.path(out_dir, "03_finalTumor_Normal_UMAP.pdf"),
  p_class | p_score,
  width = 13,
  height = 6
)

if (nrow(stability_log) > 0) {
  p_stability <- ggplot(
    stability_log,
    aes(x = Iteration, y = Label_stability, group = 1)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = LABEL_STABLE, linetype = "dashed") +
    ylim(0, 1) +
    theme_classic() +
    labs(y = "Label stability", title = "Phenotype iteration convergence") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    file.path(out_dir, "03_phenotype_iteration_stability.pdf"),
    p_stability,
    width = 7,
    height = 4.5
  )
}

malignant_obj <- subset(epi_obj, subset = tumorselect_label == "Tumor")
normal_obj <- subset(epi_obj, subset = tumorselect_label == "Normal")

epi_obj@misc$TumorSelect <- list(
  method = "SCIPAC ER+/HER2- seed plus phenotype iteration",
  final_iteration = final_iter,
  final_tumor_signature = current_sig$tumor,
  final_normal_signature = current_sig$normal,
  stability_log = stability_log
)

saveRDS(
  epi_obj,
  file.path(out_dir, "GSE245601_Epithelial_TumorSelect.rds")
)
saveRDS(
  malignant_obj,
  file.path(out_dir, "GSE245601_Malignant_Epithelial.rds")
)
saveRDS(
  normal_obj,
  file.path(out_dir, "GSE245601_NonMalignant_Epithelial.rds")
)

write.csv(
  data.frame(
    Cell = colnames(epi_obj),
    Sample = if ("sample_gsm" %in% colnames(epi_obj[[]])) {
      epi_obj$sample_gsm
    } else {
      epi_obj$orig.ident
    },
    Cluster = epi_obj$seurat_clusters,
    SCIPAC_seed = epi_obj$scipac_seed,
    Final_label = epi_obj$tumorselect_label,
    Final_score = epi_obj$tumorselect_diff_score
  ),
  file.path(out_dir, "04_cellfinal.csv"),
  row.names = FALSE
)

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(fastCNV)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
epi_rds <- file.path(base_dir, "08.tumorselect",
                     "GSE245601_Epithelial_TumorSelect.rds")
full_rds <- file.path(base_dir, "07.celltype_UMAP",
                      "GSE245601_Step5_Manual_Annotated.rds")
out_dir <- file.path(base_dir, "09.fastcnv")
run_dir <- file.path(out_dir, "FastCNV_by_sample")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(epi_rds, full_rds)) {
  if (!file.exists(f)) stop("not_foundinput: ", f)
}

ref_candidates <- c("T cells", "B cells", "Macrophages",
                    "DCs", "Mast cells", "Endothelial cells")
max_ref <- 5000L
threshold_q <- 0.99

epi_seu <- readRDS(epi_rds)
full_obj <- readRDS(full_rds)
DefaultAssay(epi_seu) <- "RNA"
DefaultAssay(full_obj) <- "RNA"

if (!"tumorselect_label" %in% colnames(epi_seu[[]])) {
  stop("missingtumorselect_label")
}
if (!"manual_celltype" %in% colnames(full_obj[[]])) {
  stop("missingmanual_celltype")
}
if (!"umap" %in% names(epi_seu@reductions)) stop("missingumap")

if (length(Layers(epi_seu[["RNA"]])) > 1) {
  epi_seu <- JoinLayers(epi_seu, assay = "RNA")
}
if (length(Layers(full_obj[["RNA"]])) > 1) {
  full_obj <- JoinLayers(full_obj, assay = "RNA")
}

sample_col <- if ("sample_gsm" %in% colnames(epi_seu[[]])) {
  "sample_gsm"
} else {
  "orig.ident"
}
if (!sample_col %in% colnames(epi_seu[[]])) stop("missingsample")

print(table(epi_seu$tumorselect_label, useNA = "ifany"))

ref_types <- intersect(ref_candidates,
                       unique(as.character(full_obj$manual_celltype)))
if (length(ref_types) < 2) stop("2")

max_each <- floor(max_ref / length(ref_types))
ref_cells <- unlist(lapply(ref_types, function(ct) {
  cells <- colnames(full_obj)[full_obj$manual_celltype == ct]
  sample(cells, min(length(cells), max_each))
}), use.names = FALSE)

common_genes <- intersect(rownames(epi_seu), rownames(full_obj))
if (length(common_genes) < 500) stop("gene500")

ref_obj <- subset(full_obj, features = common_genes, cells = ref_cells)
ref_obj$fastcnv_ref_type <- as.character(ref_obj$manual_celltype)
ref_map <- setNames(ref_obj$fastcnv_ref_type, colnames(ref_obj))

write.csv(
  data.frame(Cell = colnames(ref_obj),
             ReferenceType = ref_obj$fastcnv_ref_type),
  file.path(out_dir, "00_FastCNVcell.csv"),
  row.names = FALSE
)

rm(full_obj)
invisible(gc())
