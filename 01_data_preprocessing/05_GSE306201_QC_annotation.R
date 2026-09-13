# GSE306201 QC, integration and annotation

PROJECT_DIR <- "."

options(stringsAsFactors = FALSE)

library(Seurat)
library(Matrix)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(org.Hs.eg.db)
library(AnnotationDbi)

base_dir      <- file.path(PROJECT_DIR, "data/raw/GSE306201")
out_dir       <- file.path(PROJECT_DIR, "results/GSE306201")
features_file <- file.path(base_dir, "GSE306201_Features.tsv.gz")
raw_dir       <- file.path(base_dir, "GSE306201_RAW")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feat <- read_tsv(features_file, col_names = FALSE, show_col_types = FALSE)

if (any(grepl("feature", paste(feat[1, ], collapse = ""))) ||
    any(grepl("gene",    paste(feat[1, ], collapse = ""))) ||
    all(feat[1, ] %in% 0:2)) {
  feat <- feat[-1, ]
}
feat <- feat[complete.cases(feat), ]

gene_ids        <- feat[[1]]
gene_sym        <- feat[[2]]
rownames_to_use <- gene_sym

sample_ids_tmp <- list.dirs(raw_dir, full.names = FALSE, recursive = FALSE)
sample_ids_tmp <- sample_ids_tmp[sample_ids_tmp != ""]
if (length(sample_ids_tmp) == 0) stop("RAW directory GSM directory")

set.seed(123)
check_samples <- sample(sample_ids_tmp, min(3, length(sample_ids_tmp)))
same_nrow <- TRUE
for (sid in check_samples) {
  mtx_file <- file.path(raw_dir, sid, "matrix.mtx.gz")
  if (!file.exists(mtx_file)) { cat("  [] sample", sid, "missing matrix.mtx.gz\n"); next }
  test_mtx <- readMM(mtx_file)
  if (nrow(test_mtx) != nrow(feat)) {
    same_nrow <- FALSE
  }
}
if (same_nrow) {
} else {
}

feat_sym     <- gene_sym
mt_in_feat   <- feat_sym[grepl("^MT-", feat_sym)]
rb_in_feat   <- feat_sym[grepl("^RPL", feat_sym) | grepl("^RPS", feat_sym)]

sample_ids <- list.dirs(raw_dir, full.names = FALSE, recursive = FALSE)
sample_ids <- sample_ids[sample_ids != ""]
if (length(sample_ids) == 0) stop("RAW directory GSM directory")
print(sample_ids)

read_one_sample <- function(sample_id) {
  sample_dir <- file.path(raw_dir, sample_id)
  mtx_file   <- file.path(sample_dir, "matrix.mtx.gz")
  bc_file    <- file.path(sample_dir, "barcodes.tsv.gz")
  if (!file.exists(mtx_file)) stop("missing matrix : ", mtx_file)
  if (!file.exists(bc_file))  stop("missing barcodes : ", bc_file)
  mat <- readMM(mtx_file)
  bcs <- read_tsv(bc_file, col_names = FALSE, show_col_types = FALSE)[[1]]
  if (nrow(mat) != length(rownames_to_use)) {
    stop("sample ", sample_id, "  (", nrow(mat),
         ")  features  (", length(rownames_to_use), ") ")
  }
  rownames(mat) <- rownames_to_use
  colnames(mat) <- bcs
  mat
}

counts_list <- lapply(sample_ids, read_one_sample)
names(counts_list) <- sample_ids

seu_list_raw <- lapply(names(counts_list), function(nm) {
  obj <- CreateSeuratObject(
    counts     = counts_list[[nm]],
    project    = "ERplus",
    min.cells  = 3,
    min.features = 200
  )
  obj$orig.ident <- nm
  obj$sample_gsm <- nm
  return(obj)
})
names(seu_list_raw) <- sample_ids

print(sapply(seu_list_raw, ncol))

seurat_obj <- merge(
  x            = seu_list_raw[[1]],
  y            = seu_list_raw[-1],
  add.cell.ids = names(seu_list_raw)
)
print(table(seurat_obj$orig.ident))

rn <- rownames(seurat_obj)

mt_pat <- if (sum(grepl("^MT-", rn)) >= sum(grepl("^mt-", rn))) "^MT-" else "^mt-"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pat)
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^(RPL|RPS)")

mt_symbols_ref <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6",
                    "MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8","MT-CYB","MT-RNR1","MT-RNR2")
mt_genes_list  <- intersect(mt_symbols_ref, rn)
rb_genes_list  <- rn[grepl("^RPL", rn) | grepl("^RPS", rn)]
seurat_obj[["percent.mt.list"]] <- PercentageFeatureSet(seurat_obj, features = mt_genes_list)
seurat_obj[["percent.rb.list"]] <- PercentageFeatureSet(seurat_obj, features = rb_genes_list)

sym2ens     <- setNames(gene_ids, gene_sym)
ens_for_rn  <- sym2ens[rn]
ens_valid   <- unique(ens_for_rn[!is.na(ens_for_rn) & grepl("^ENSG", ens_for_rn)])
annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = ens_valid,
  columns = c("SYMBOL", "CHR"),
  keytype = "ENSEMBL"
)
mt_ens      <- unique(annot$ENSEMBL[grepl("^MT-", annot$SYMBOL)])
rb_ens      <- unique(annot$ENSEMBL[grepl("^RPL", annot$SYMBOL) | grepl("^RPS", annot$SYMBOL)])
mt_genes_ens <- rn[ens_for_rn %in% mt_ens]
rb_genes_ens <- rn[ens_for_rn %in% rb_ens]
seurat_obj[["percent.mt.ens"]] <- PercentageFeatureSet(seurat_obj, features = mt_genes_ens)
seurat_obj[["percent.rb.ens"]] <- PercentageFeatureSet(seurat_obj, features = rb_genes_ens)

safe_cor <- function(x, y) {
  if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
  if (sd(x, na.rm=TRUE)==0 || sd(y, na.rm=TRUE)==0) return(NA_real_)
  suppressWarnings(cor(x, y, use = "complete.obs"))
}

group_map <- c(
  "GSM9194698_P" = "Primary",
  "GSM9194699_P" = "Primary",
  "GSM9194700_P" = "Primary",
  "GSM9194701_P" = "Primary",
  "GSM9194702_P" = "Primary",
  "GSM9194703_P" = "Primary",
  "GSM9194714_P" = "Primary",
  "GSM9194716_P10" = "Primary",
  "GSM9194716_P8" = "Primary",
  "GSM9194718_P11" = "Primary",
  "GSM9194718_P12" = "Primary",
  "GSM9194718_P9" = "Primary",
  "GSM9194704_M" = "Metastatic",
  "GSM9194705_M" = "Metastatic",
  "GSM9194706_M" = "Metastatic",
  "GSM9194707_M" = "Metastatic",
  "GSM9194708_M" = "Metastatic",
  "GSM9194709_M" = "Metastatic",
  "GSM9194710_M" = "Metastatic",
  "GSM9194711_M" = "Metastatic",
  "GSM9194712_M" = "Metastatic",
  "GSM9194713_M" = "Metastatic",
  "GSM9194714_M" = "Metastatic"
)

group_map <- group_map[names(group_map) != "placeholder"]

seurat_obj$Group <- ifelse(
  seurat_obj$orig.ident %in% names(group_map),
  group_map[seurat_obj$orig.ident],
  "Unknown"
)

print(table(seurat_obj$Group))

min_features <- 200
max_features <- 6000
min_counts   <- 500
max_counts   <- 50000
mt_max       <- 10
rb_max       <- 40

qc_mask <- (seurat_obj$nFeature_RNA >= min_features) &
  (seurat_obj$nFeature_RNA <= max_features) &
  (seurat_obj$nCount_RNA   >= min_counts)   &
  (seurat_obj$nCount_RNA   <= max_counts)   &
  (seurat_obj$percent.mt   <= mt_max)        &
  (seurat_obj$percent.rb   <= rb_max)
qc_mask[is.na(qc_mask)] <- FALSE
seurat_obj$qc_pass <- qc_mask

print(table(seurat_obj$qc_pass, seurat_obj$orig.ident))

seurat_qc <- seurat_obj[, seurat_obj$qc_pass]

rm(seurat_obj, seu_list_raw, counts_list)
gc()

rna_layers <- tryCatch(Layers(seurat_qc[["RNA"]]), error = function(e) character(0))
if (length(rna_layers) > 1) {
  seurat_qc <- SeuratObject::JoinLayers(seurat_qc, assay = "RNA")
}

DefaultAssay(seurat_qc) <- "RNA"
sce    <- as.SingleCellExperiment(seurat_qc)
sce    <- scDblFinder(sce, samples = "sample_gsm", BPPARAM = SerialParam())

print(table(sce$scDblFinder.class, sce$sample_gsm))

n_doublet <- sum(sce$scDblFinder.class == "doublet")
n_total   <- length(sce$scDblFinder.class)

seurat_qc$scDblFinder.class <- sce$scDblFinder.class
seurat_qc$scDblFinder.score <- sce$scDblFinder.score
rm(sce); gc()

meta_df <- seurat_qc@meta.data

doublet_summary <- meta_df %>%
  dplyr::group_by(sample_gsm, scDblFinder.class) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sample_gsm) %>%
  dplyr::mutate(total = sum(count), pct = 100 * count / total)

p_db1 <- ggplot(doublet_summary, aes(x=sample_gsm, y=pct, fill=scDblFinder.class)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values = c("singlet"="#6fa6cf7f","doublet"="#b436658f")) +
  labs(title="Doublet Detection Summary", x="Sample (GSM)",
       y="Percentage (%)", fill="Class") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=6), legend.position="top")

p_db2 <- ggplot(meta_df, aes(x=sample_gsm, y=scDblFinder.score, fill=scDblFinder.class)) +
  geom_violin(scale="width", trim=TRUE) +
  scale_fill_manual(values = c("singlet"="#6fa6cf7f","doublet"="#b436658f")) +
  labs(title="Doublet Score Distribution", x="Sample (GSM)",
       y="Doublet Score", fill="Class") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=6), legend.position="top")

pdf(file.path(out_dir, "QC_01_Doublet_Detection.pdf"), width=14, height=6)
print(p_db1 + p_db2)
dev.off()

seurat_singlet <- seurat_qc[, seurat_qc$scDblFinder.class == "singlet"]
rm(seurat_qc, meta_df, doublet_summary, p_db1, p_db2); gc()

p_vln <- VlnPlot(
  seurat_singlet,
  features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),
  group.by = "sample_gsm",
  pt.size  = 0,
  ncol     = 2
)
pdf(file.path(out_dir, "QC_02_VlnPlot_afterFilter.pdf"), width=14, height=10)
print(p_vln)
dev.off()

p_scat1 <- FeatureScatter(seurat_singlet, feature1="nCount_RNA",
                          feature2="nFeature_RNA", group.by="sample_gsm")
pdf(file.path(out_dir, "QC_03_Scatter_Count_Feature.pdf"), width=8, height=6)
print(p_scat1)
dev.off()

p_scat2 <- FeatureScatter(seurat_singlet, feature1="nCount_RNA",
                          feature2="percent.mt", group.by="sample_gsm")
pdf(file.path(out_dir, "QC_04_Scatter_Count_MT.pdf"), width=8, height=6)
print(p_scat2)
dev.off()

saveRDS(seurat_singlet,
        file = file.path(out_dir, "GSE306201_Step1_QC_Singlet.rds"))

options(stringsAsFactors = FALSE)

library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(celda)
library(SingleCellExperiment)
library(future)
library(future.apply)

out_dir  <- file.path(PROJECT_DIR, "results/GSE306201")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan("multisession", workers = 6)
options(future.globals.maxSize = 10000 * 1024^2)

seurat_obj <- readRDS(file.path(out_dir, "GSE306201_Step1_QC_Singlet.rds"))
print(table(seurat_obj$orig.ident))
print(table(seurat_obj$Group))

plan(sequential)

rna_layers <- tryCatch(Layers(seurat_obj[["RNA"]]), error=function(e) character(0))
if (length(rna_layers) > 1) {
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay="RNA")
}

counts_raw <- GetAssayData(seurat_obj, assay="RNA", layer="counts")

sce_decont <- SingleCellExperiment(
  assays  = list(counts = counts_raw),
  colData = data.frame(sample_gsm = seurat_obj$sample_gsm,
                       row.names  = colnames(seurat_obj))
)
sce_decont <- decontX(sce_decont, batch = sce_decont$sample_gsm)

seurat_obj$decontX_contamination <- colData(sce_decont)$decontX_contamination
print(quantile(seurat_obj$decontX_contamination, probs=seq(0,1,0.1)))

p_contam <- ggplot(
  data.frame(contamination = seurat_obj$decontX_contamination),
  aes(x = contamination)
) +
  geom_histogram(bins=60, fill="#6fa6cf8f", color="white", alpha=0.85) +
  geom_vline(xintercept=0.2, linetype="dashed", color="#b436659f", linewidth=0.8) +
  annotate("text", x=0.22, y=Inf, vjust=2, hjust=0,
           label="Threshold = 0.2", color="#b436659f", size=3.5) +
  labs(title  = "DecontX Contamination Distribution",
       subtitle = paste0("n = ", ncol(seurat_obj), " cells (after QC & doublet removal)"),
       x = "Contamination Proportion",
       y = "Cell Count") +
  theme_bw(base_size=12) +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, color="grey50"))

pdf(file.path(out_dir, "INT_01_decontX_Contamination.pdf"), width=6, height=4)
print(p_contam)
dev.off()

decontX_counts <- round(assay(sce_decont, "decontXcounts"))
seurat_obj[["RNA"]]$counts <- decontX_counts

n_before <- ncol(seurat_obj)
seurat_obj <- seurat_obj[, seurat_obj$decontX_contamination < 0.2]
n_after  <- ncol(seurat_obj)

rm(sce_decont, counts_raw, decontX_counts, p_contam); gc()

plan("multisession", workers = 6)

seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = "LogNormalize",
                            scale.factor         = 10000)

seurat_obj <- FindVariableFeatures(seurat_obj,
                                   selection.method = "vst",
                                   nfeatures        = 2000)

top10    <- head(VariableFeatures(seurat_obj), 10)
p_hvg    <- VariableFeaturePlot(seurat_obj)
p_hvg    <- LabelPoints(plot=p_hvg, points=top10, repel=TRUE, xnudge=0, ynudge=0)
pdf(file.path(out_dir, "INT_02_VariableFeatures.pdf"), width=8, height=5)
print(p_hvg)
dev.off()

s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj <- CellCycleScoring(seurat_obj,
                               s.features   = s.genes,
                               g2m.features = g2m.genes,
                               set.ident    = FALSE)
seurat_obj$Phase <- factor(seurat_obj$Phase, levels=c("G1","S","G2M"))

p_cc <- VlnPlot(seurat_obj,
                features = c("S.Score","G2M.Score"),
                group.by = "orig.ident",
                pt.size  = 0, ncol = 2)
pdf(file.path(out_dir, "INT_03_CellCycle_Score.pdf"), width=12, height=5)
print(p_cc)
dev.off()
print(table(seurat_obj$Phase))

plan(sequential)

hvg_features <- VariableFeatures(seurat_obj)

seurat_obj <- ScaleData(seurat_obj,
                        features       = hvg_features,
                        vars.to.regress = c("S.Score","G2M.Score","percent.mt"))

plan("multisession", workers = 6)

seurat_obj <- RunPCA(seurat_obj,
                     features = VariableFeatures(seurat_obj),
                     npcs     = 50,
                     verbose  = FALSE)

p_elbow <- ElbowPlot(seurat_obj, ndims=50)
pdf(file.path(out_dir, "INT_04_PCA_ElbowPlot.pdf"), width=6, height=4)
print(p_elbow)
dev.off()

pdf(file.path(out_dir, "INT_05_PCA_Heatmap.pdf"), width=14, height=10)
DimHeatmap(seurat_obj, dims=1:15, cells=500, balanced=TRUE)
dev.off()

p_pca_pre <- DimPlot(seurat_obj, reduction="pca",
                     group.by="orig.ident", pt.size=0.1) +
  ggtitle("PCA before Harmony (colored by sample)") +
  theme_bw()
pdf(file.path(out_dir, "INT_06_PCA_preBatchCorrection.pdf"), width=9, height=6)
print(p_pca_pre)
dev.off()

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars  = "orig.ident",
  reduction      = "pca",
  reduction.save = "harmony",
  plot_convergence = TRUE,
  max.iter.harmony = 30
)

p_harm_sample <- DimPlot(seurat_obj, reduction="harmony",
                         group.by="orig.ident", pt.size=0.1) +
  ggtitle("Harmony embedding (colored by sample)") +
  theme_bw()

p_harm_group  <- DimPlot(seurat_obj, reduction="harmony",
                         group.by="Group", pt.size=0.1,
                         cols=c("Primary"="#b436659f","Metastatic"="#6fa6cf9f")) +
  ggtitle("Harmony embedding (colored by Group)") +
  theme_bw()

pdf(file.path(out_dir, "INT_07_Harmony_Embedding.pdf"), width=14, height=6)
print(p_harm_sample | p_harm_group)
dev.off()

saveRDS(seurat_obj,
        file = file.path(out_dir, "GSE306201_Step2_Harmony.rds"))

print(table(seurat_obj$Group))

options(stringsAsFactors = FALSE)

library(Seurat)
library(SeuratObject)
library(SCP)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(future)

out_dir <- file.path(PROJECT_DIR, "results/GSE306201")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plan("multisession", workers = 6)
options(future.globals.maxSize = 10000 * 1024^2)

seurat_obj <- readRDS(file.path(out_dir, "GSE306201_Step2_Harmony.rds"))
print(table(seurat_obj$orig.ident))
print(table(seurat_obj$Group))

if (!"harmony" %in% names(seurat_obj@reductions)) {
  stop("[ERROR] harmony reduction ,  02_Integration.R")
}

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims      = 1:30,
  k.param   = 20,
  verbose   = FALSE
)

plan(sequential)

res_seq <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
for (res in res_seq) {
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution   = res,
    algorithm    = 1,
    random.seed  = 42,
    verbose      = FALSE
  )
}

pdf(file.path(out_dir, "CLU_01_Clustree.pdf"), width = 12, height = 10)
print(clustree(seurat_obj@meta.data, prefix = "RNA_snn_res."))
dev.off()

selected_res <- 0.2

seurat_obj <- FindClusters(
  seurat_obj,
  resolution  = selected_res,
  algorithm   = 1,
  random.seed = 42,
  verbose     = FALSE
)

cluster_tab <- table(seurat_obj$seurat_clusters)
print(cluster_tab)

res_cols_all    <- grep("RNA_snn_res\\.", colnames(seurat_obj@meta.data), value = TRUE)
res_col_keep    <- paste0("RNA_snn_res.", selected_res)
res_cols_remove <- setdiff(res_cols_all, res_col_keep)
for (col in res_cols_remove) seurat_obj[[col]] <- NULL

plan("multisession", workers = 6)

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction     = "harmony",
  dims          = 1:10,
  n.neighbors   = 30,
  min.dist      = 0.3,
  spread        = 1,
  seed.use      = 42,
  reduction.name = "umap",
  verbose       = FALSE
)

library(SCP)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)

n_clusters  <- nlevels(seurat_obj$seurat_clusters)
total_cells <- ncol(seurat_obj)

if (n_clusters <= 8) {
  clust_cols <- brewer.pal(n_clusters, "Set2")
} else if (n_clusters <= 12) {
  clust_cols <- brewer.pal(n_clusters, "Set3")
} else {
  clust_cols <- colorRampPalette(
    c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"))
  )(n_clusters)
}
names(clust_cols) <- levels(seurat_obj$seurat_clusters)

n_samples    <- length(unique(seurat_obj$orig.ident))
scp_col_samp <- colorRampPalette(
  c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))
)(n_samples)
names(scp_col_samp) <- sort(unique(seurat_obj$orig.ident))

umap_df         <- as.data.frame(Embeddings(seurat_obj, "umap"))
umap_df$cluster <- seurat_obj$seurat_clusters

centers <- umap_df %>%
  group_by(cluster) %>%
  summarise(
    UMAP_1 = mean(umap_1),
    UMAP_2 = mean(umap_2),
    .groups = "drop"
  )

theme_box <- theme_scp() +
  theme(
    plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title  = element_text(size = 8,  face = "bold"),
    legend.text   = element_text(size = 7.5),
    aspect.ratio  = 1
  )

p1_left <- CellDimPlot(
  seurat_obj,
  group.by  = "seurat_clusters",
  reduction = "umap",
  cols      = clust_cols,
  pt.size   = 0.01,
  alpha     = 0.75,
  label     = FALSE
) +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    ncol = 1
  )) +
  ggtitle("UMAP - Seurat Clusters") +
  theme_box

p1_right <- CellDimPlot(
  seurat_obj,
  group.by  = "seurat_clusters",
  reduction = "umap",
  cols      = clust_cols,
  pt.size   = 0.01,
  alpha     = 0.75,
  label     = FALSE
) +
  geom_text(
    data        = centers,
    aes(x = UMAP_1, y = UMAP_2, label = cluster),
    color       = "black",
    size        = 3.5,
    fontface    = "bold",
    inherit.aes = FALSE
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    ncol = 1
  )) +
  ggtitle("UMAP - Seurat Clusters (Labeled)") +
  theme_box

p1_combined <- p1_left + p1_right +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

pdf(file.path(out_dir, "CLU_02_03_UMAP_Clusters.pdf"), width = 13, height = 6)
print(p1_combined)
dev.off()

theme_arrow <- theme(
  plot.title        = element_text(size = 12, face = "bold", hjust = 0.5),
  legend.title      = element_text(size = 8,  face = "bold"),
  legend.text       = element_text(size = 7.5),
  panel.border      = element_blank(),
  panel.background  = element_blank(),
  legend.background = element_blank(),
  legend.key        = element_blank(),
  axis.line         = element_line(
    arrow     = arrow(length = unit(0.18, "inches"), type = "closed"),
    linewidth = 0.55,
    color     = "black"
  ),
  axis.ticks  = element_blank(),
  axis.text   = element_blank(),
  axis.title  = element_text(face = "bold", size = 9),
  aspect.ratio = 1
)

p2 <- CellDimPlot(
  seurat_obj,
  group.by   = "seurat_clusters",
  reduction  = "umap",
  cols       = clust_cols,
  pt.size    = 0.01,
  alpha      = 0.75,
  label      = FALSE,
  theme_use  = "theme_blank"
) +
  geom_text(
    data        = centers,
    aes(x = UMAP_1, y = UMAP_2, label = cluster),
    color       = "black",
    size        = 3.5,
    fontface    = "bold",
    inherit.aes = FALSE
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    ncol = 1
  )) +
  ggtitle("UMAP - Seurat Clusters") +
  theme_arrow

pdf(file.path(out_dir, "CLU_04_UMAP_Arrow.pdf"), width = 7, height = 6)
print(p2)
dev.off()

p3 <- CellDimPlot(
  seurat_obj,
  group.by  = "orig.ident",
  reduction = "umap",
  cols      = scp_col_samp,
  pt.size   = 0.01,
  alpha     = 0.65,
  label     = FALSE
) +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    ncol = ceiling(n_samples / 15)
  )) +
  ggtitle("UMAP - Sample Origin") +
  theme_scp() +
  theme(
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 8,  face = "bold"),
    legend.text  = element_text(size = 7),
    aspect.ratio = 1
  )

pdf(file.path(out_dir, "CLU_05_UMAP_Sample.pdf"), width = 10, height = 7)
print(p3)
dev.off()

p4 <- CellDimPlot(
  seurat_obj,
  group.by  = "Group",
  reduction = "umap",
  cols      = group_colors,
  pt.size   = 0.01,
  alpha     = 0.65,
  label     = FALSE
) +
  guides(color = guide_legend(
    override.aes = list(size = 5, alpha = 1),
    ncol = 1
  )) +
  ggtitle("UMAP - Primary vs Metastatic") +
  theme_scp() +
  theme(
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 8,  face = "bold"),
    legend.text  = element_text(size = 8),
    aspect.ratio = 1
  )

pdf(file.path(out_dir, "CLU_06_UMAP_Group.pdf"), width = 7, height = 6)
print(p4)
dev.off()

p5 <- CellDimPlot(
  seurat_obj,
  group.by  = "seurat_clusters",
  split.by  = "Group",
  reduction = "umap",
  cols      = clust_cols,
  pt.size   = 0.01,
  alpha     = 0.7,
  label     = FALSE,
  ncol      = 2
) +
  ggtitle("UMAP Split by Group (Primary vs Metastatic)") +
  theme_scp() +
  theme(
    plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
    strip.text   = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 8,  face = "bold"),
    legend.text  = element_text(size = 7.5),
    aspect.ratio = 1
  )

pdf(file.path(out_dir, "CLU_07_UMAP_SplitByGroup.pdf"), width = 13, height = 6)
print(p5)
dev.off()

theme_qc <- theme_scp() +
  theme(
    plot.title   = element_text(size = 10, face = "bold", hjust = 0.5),
    aspect.ratio = 1
  )

pG1 <- FeatureDimPlot(
  seurat_obj, features = "percent.mt",
  reduction = "umap", pt.size = 0.01,
  theme_use = "theme_scp"
) + labs(title = "Mitochondrial %") + theme_qc

pG2 <- FeatureDimPlot(
  seurat_obj, features = "nFeature_RNA",
  reduction = "umap", pt.size = 0.01,
  theme_use = "theme_scp"
) + labs(title = "nFeature_RNA") + theme_qc

pG3 <- FeatureDimPlot(
  seurat_obj, features = "nCount_RNA",
  reduction = "umap", pt.size = 0.01,
  theme_use = "theme_scp"
) + labs(title = "nCount_RNA") + theme_qc

pG4 <- FeatureDimPlot(
  seurat_obj, features = "decontX_contamination",
  reduction = "umap", pt.size = 0.01,
  theme_use = "theme_scp"
) + labs(title = "decontX Contamination") + theme_qc

pdf(file.path(out_dir, "CLU_08_UMAP_QC.pdf"), width = 12, height = 11)
print((pG1 | pG2) / (pG3 | pG4))
dev.off()

pH1 <- DimPlot(
  seurat_obj,
  reduction = "pca",
  group.by  = "orig.ident",
  pt.size   = 0.01,
  cols      = scp_col_samp
) +
  ggtitle("Before Harmony: PCA") +
  theme_scp() +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, size = 11, face = "bold"),
    aspect.ratio    = 1
  )

pH2 <- DimPlot(
  seurat_obj,
  reduction = "harmony",
  group.by  = "orig.ident",
  pt.size   = 0.01,
  cols      = scp_col_samp
) +
  ggtitle("After Harmony") +
  theme_scp() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 11, face = "bold"),
    legend.text  = element_text(size = 7),
    aspect.ratio = 1
  )

pdf(file.path(out_dir, "CLU_09_Harmony_Compare.pdf"), width = 12, height = 5)
print(pH1 | pH2)
dev.off()

cluster_summary <- seurat_obj@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    n_cells      = dplyr::n(),
    n_primary    = sum(Group == "Primary",    na.rm = TRUE),
    n_metastatic = sum(Group == "Metastatic", na.rm = TRUE),
    pct_primary  = round(100 * n_primary    / dplyr::n(), 1),
    pct_meta     = round(100 * n_metastatic / dplyr::n(), 1),
    .groups      = "drop"
  ) %>%
  dplyr::arrange(as.integer(as.character(seurat_clusters)))

print(cluster_summary, n = Inf)
write.csv(
  cluster_summary,
  file      = file.path(out_dir, "CLU_10_Cluster_Summary.csv"),
  row.names = FALSE
)

saveRDS(
  seurat_obj,
  file = file.path(out_dir, "GSE306201_Step3_Clustered_UMAP.rds")
)

library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(SCP)
library(SCpubr)

anno_dir <- "SingleR"
if (!dir.exists(anno_dir)) dir.create(anno_dir, recursive = TRUE)

seurat_obj <- readRDS("GSE306201_Step3_Clustered_UMAP.rds")

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

if (!"data" %in% Layers(seurat_obj[["RNA"]])) {
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor         = 1e4,
    verbose              = FALSE
  )
}

expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
cluster_ids <- seurat_obj$seurat_clusters

ref_hpa <- celldex::HumanPrimaryCellAtlasData()

singler_result <- tryCatch({
  SingleR(
    test     = expr_matrix,
    ref      = ref_hpa,
    labels   = ref_hpa$label.main,
    clusters = cluster_ids
  )
}, error = function(e) {
  NULL
})

if (!is.null(singler_result)) {

  cluster_to_type <- singler_result$labels
  names(cluster_to_type) <- rownames(singler_result)

  print(data.frame(
    Cluster  = names(cluster_to_type),
    CellType = cluster_to_type,
    row.names = NULL
  ))

  cell_types_vec <- setNames(
    cluster_to_type[as.character(cluster_ids)],
    colnames(seurat_obj)
  )

  seurat_obj$SingleR_celltype <- cell_types_vec

  print(sort(table(seurat_obj$SingleR_celltype), decreasing = TRUE))

  write.csv(
    data.frame(
      Cluster  = names(cluster_to_type),
      CellType = cluster_to_type,
      row.names = NULL
    ),
    file      = file.path(anno_dir, "CellCluster_AutoAnno_SingleR.csv"),
    row.names = FALSE
  )

} else {
  seurat_obj$SingleR_celltype <- paste0("Cluster_", as.character(cluster_ids))
}

seurat_obj$SingleR_celltype <- factor(
  seurat_obj$SingleR_celltype,
  levels = sort(unique(seurat_obj$SingleR_celltype))
)

p_square <- SCpubr::do_DimPlot(
  sample          = seurat_obj,
  group.by        = "SingleR_celltype",
  reduction       = "umap",
  label           = TRUE,
  label.box       = TRUE,
  label.size      = 3.5,
  repel           = TRUE,
  legend.position = "right",
  pt.size         = 0.08,
  border.size     = 0.8,
  font.size       = 11
) +
  coord_fixed() +
  theme(
    aspect.ratio = 1,
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.text  = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  ) +
  ggtitle("UMAP - SingleR Cell Type Annotation")

pdf(file.path(anno_dir, "Anno_01_UMAP_Square_Labeled.pdf"),
    width = 10, height = 8)
print(p_square)
dev.off()

p_arrow <- CellDimPlot(
  srt       = seurat_obj,
  group.by  = "SingleR_celltype",
  reduction = "umap",
  pt.size   = 0.05,
  alpha     = 0.8,
  theme_use = "theme_blank"
) +
  guides(color = guide_legend(
    override.aes = list(size = 4, alpha = 1),
    title        = "Cell Type"
  )) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.text = element_text(size = 9)
  ) +
  ggtitle("UMAP - SingleR Cell Type Annotation")

pdf(file.path(anno_dir, "Anno_02_UMAP_ArrowAxes.pdf"),
    width = 9, height = 7)
print(p_arrow)
dev.off()

saveRDS(seurat_obj, file = "GSE306201_Step4_Annotated.rds")

library(Seurat)
library(SCP)
library(ggplot2)
library(dplyr)

panel_dir <- "Marker_Panel"
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

seurat_obj <- readRDS("GSE306201_Step3_Clustered_UMAP.rds")

stopifnot("RNA" %in% names(seurat_obj@assays))
stopifnot("seurat_clusters" %in% colnames(seurat_obj@meta.data))
stopifnot("umap" %in% names(seurat_obj@reductions))

DefaultAssay(seurat_obj) <- "RNA"

if (!"data" %in% Layers(seurat_obj[["RNA"]])) {
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor         = 1e4,
    verbose              = FALSE
  )
} else {
}

marker_panel <- list(

  Epithelial = c("EPCAM", "KRT8", "KRT18", "ESR1", "GATA3"),

  Proliferating_cells = c("MKI67", "TOP2A", "PCNA"),

  Myeloid_cells = c("CD68", "CD14", "LYZ", "C1QA"),

  T_cells = c("CD3D", "CD3E", "TRAC"),

  NK_cells = c("NKG7", "KLRD1", "GNLY"),

  B_cells = c("MS4A1", "CD79A", "CD19"),

  Plasma_cells = c("JCHAIN", "MZB1", "IGKC"),

  Mast_cells = c("TPSAB1", "CPA3", "KIT"),

  Fibroblasts_CAFs = c("DCN", "COL1A1", "FAP", "PDGFRB"),

  Endothelial_cells = c("PECAM1", "VWF", "CDH5"),

  Mural_cells_SMCs = c("ACTA2", "RGS5", "MCAM")
)

marker_panel2 <- list(
  Epithelial_cells = c("EPCAM", "KRT8", "KRT18", "ESR1"),
  Proliferating_cells    = c("MKI67", "TOP2A", "CENPF", "PCNA"),
  T_cells                = c("CD3D", "CD3E", "TRBC1", "CD2"),
  NK_cells               = c("NKG7", "GNLY", "XCL1", "KLRD1"),
  B_cells                = c("MS4A1", "CD79A", "CD74", "CD19"),
  Plasma_cells           = c("JCHAIN", "MZB1", "SDC1", "XBP1"),
  Myeloid_cells          = c("CD68", "LST1", "TYROBP", "LYZ"),
  Mast_cells             = c("TPSAB1", "TPSB2", "KIT", "CPA3"),
  Endothelial_cells      = c("PECAM1", "VWF", "KDR", "CDH5"),
  Mesenchymal_cells      = c("PDGFRB", "COL1A1", "DCN", "LUM")
)

feat_vec_panel <- unique(unlist(marker_panel))

genes_present <- feat_vec_panel[feat_vec_panel %in% rownames(seurat_obj)]
genes_missing <- setdiff(feat_vec_panel, genes_present)

if (length(genes_missing) > 0) {
}

feat_vec_panel <- genes_present

pdf(file.path(panel_dir, "Panel_01_DotPlot_by_cluster.pdf"), width = 16, height = 8)
p1 <- DotPlot(
  seurat_obj,
  features = feat_vec_panel,
  assay    = "RNA",
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  ggtitle("Marker Gene Expression per Cluster (DotPlot)") +
  theme_scp() +
  theme(
    axis.text.x = element_text(size = 9, face = "italic"),
    axis.text.y = element_text(size = 10)
  )
print(p1)
dev.off()

if ("SingleR_celltype" %in% colnames(seurat_obj@meta.data)) {

  pdf(file.path(panel_dir, "Panel_02_DotPlot_by_SingleR.pdf"), width = 6, height = 10)
  p2 <- DotPlot(
    seurat_obj,
    features = feat_vec_panel,
    assay    = "RNA",
    group.by = "SingleR_celltype"
  ) +
    coord_flip() +
    RotatedAxis() +
    ggtitle("Marker Gene Expression per Cell Type (DotPlot)") +
    theme_scp() +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9, face = "italic")
    )
  print(p2)
  dev.off()
} else {
}

n_genes <- length(feat_vec_panel)
n_cols  <- 4
n_rows  <- ceiling(n_genes / n_cols)

pdf(
  file.path(panel_dir, "Panel_03_FeaturePlot_UMAP.pdf"),
  width  = n_cols * 6,
  height = n_rows * 5
)
p3 <- FeaturePlot(
  seurat_obj,
  features  = feat_vec_panel,
  reduction = "umap",
  slot      = "data",
  ncol      = n_cols,
  pt.size   = 0.01
) &
  theme_scp() &
  theme(
    plot.title = element_text(face = "italic", size = 12),
    legend.position = "right"
  )
print(p3)
dev.off()

panel_df <- data.frame(
  Cell_Type = rep(names(marker_panel), sapply(marker_panel, length)),
  Gene      = unlist(marker_panel),
  Present   = unlist(marker_panel) %in% rownames(seurat_obj),
  row.names = NULL
)

write.csv(
  panel_df,
  file      = file.path(panel_dir, "Panel_00_Gene_List.csv"),
  row.names = FALSE
)

if ("SingleR_celltype" %in% colnames(seurat_obj@meta.data)) {
}

library(Seurat)
library(SCP)
library(ggplot2)
library(dplyr)
library(ggrepel)

panel_dir <- "Marker_Panel"
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

if (!exists("seurat_obj")) {
  seurat_obj <- readRDS("GSE306201_Step3_Clustered_UMAP.rds")
}

Idents(seurat_obj) <- seurat_obj$seurat_clusters

cluster_to_type <- c(
  "0"  = "Epithelial cells",
  "1"  = "Mesenchymal cells",
  "2"  = "T cells",
  "3"  = "Myeloid cells",
  "4"  = "Epithelial cells",
  "5"  = "Epithelial cells",
  "6"  = "Endothelial cells",
  "7"  = "Mesenchymal cells",
  "8"  = "Mesenchymal cells",
  "9"  = "Epithelial cells",
  "10" = "NK cells",
  "11" = "B cells",
  "12" = "Epithelial cells",
  "13" = "B cells",
  "14" = "B cells",
  "15" = "Epithelial cells",
  "16" = "Epithelial cells",
  "17" = "Mesenchymal cells",
  "18" = "Epithelial cells"
)

Idents(seurat_obj) <- seurat_obj$seurat_clusters

seurat_obj <- RenameIdents(seurat_obj, cluster_to_type)

seurat_obj$manual_celltype <- as.character(Idents(seurat_obj))

print(table(seurat_obj$manual_celltype))

all_celltypes <- sort(unique(seurat_obj$manual_celltype))

color_palette <- c(
  "Epithelial cells"       = "#b436657f",
  "Mesenchymal cells"      = "#6fa6cf9f",
  "Myeloid cells"          = "#e580279f",
  "T cells"                = "#187d799f",
  "NK cells"               = "#efb4219f",
  "Endothelial cells"      = "#e7097f7f",
  "B cells"                = "#8560af9f"
)

clusterCols <- color_palette[all_celltypes]

missing_types <- all_celltypes[!all_celltypes %in% names(color_palette)]
if (length(missing_types) > 0) {
  extra_cols <- setNames(
    scales::hue_pal()(length(missing_types)),
    missing_types
  )
  clusterCols <- c(clusterCols, extra_cols)
}

pdf(file.path(panel_dir, "Annot_01_UMAP_SCP_clean.pdf"), width = 7, height = 5)
p1 <- CellDimPlot(
  seurat_obj,
  group.by   = "manual_celltype",
  cols       = clusterCols,
  pt.size    = 0.01,
  alpha      = 0.8,
  theme_use  = "theme_blank"
) +
  guides(color = guide_legend(
    override.aes = list(size = 5, alpha = 1),
    ncol = 1
  )) +
  ggtitle("HR+/HER2- Breast Cancer - Cell Type Annotation") +
  theme(
    plot.title       = element_text(size = 13, hjust = 0.5, face = "bold"),
    panel.border     = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    legend.key       = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank()
  )
print(p1)
dev.off()

pdf(file.path(panel_dir, "Annot_02_UMAP_SCP_standard.pdf"), width = 7, height = 5)
p2 <- CellDimPlot(
  seurat_obj,
  group.by  = "manual_celltype",
  cols      = clusterCols,
  pt.size   = 0.01,
  alpha     = 0.8,
  label     = TRUE,
  label_repel = TRUE,
  label_insitu = FALSE
) +
  guides(color = guide_legend(
    override.aes = list(size = 5, alpha = 1),
    ncol = 1
  )) +
  ggtitle("HR+/HER2- Breast Cancer - Cell Type Annotation") +
  theme_scp() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold")
  )
print(p2)
dev.off()

library(ggplot2)
library(ggrepel)
library(dplyr)

plot_df <- cbind(
  as.data.frame(seurat_obj@reductions$umap@cell.embeddings),
  seurat_obj@meta.data
)
colnames(plot_df)[1:2] <- c("umap_1", "umap_2")

plot_df$manual_celltype <- factor(
  plot_df$manual_celltype,
  levels = names(clusterCols)
)

celltype_pos <- plot_df %>%
  group_by(manual_celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

xmin <- min(plot_df$umap_1); xmax <- max(plot_df$umap_1)
ymin <- min(plot_df$umap_2); ymax <- max(plot_df$umap_2)
xrange <- xmax - xmin;       yrange <- ymax - ymin

ax_len_x <- xrange * 0.22
ax_len_y <- yrange * 0.22

ax_x0 <- xmin - xrange * 0.02
ax_y0 <- ymin - yrange * 0.02

pA <- ggplot(plot_df, aes(x = umap_1, y = umap_2)) +

  geom_point(aes(color = manual_celltype), size = 0.01, alpha = 0.6) +

  stat_density_2d(
    aes(color = manual_celltype),
    geom     = "density_2d",
    linewidth = 0.3,
    linetype  = "dashed",
    contour_var = "ndensity",
    breaks   = 0.15
  ) +

  scale_color_manual(
    values = clusterCols,
    guide  = guide_legend(
      override.aes = list(size = 4, alpha = 1, linetype = 0),
      ncol = 1
    )
  ) +

  annotate("segment",
           x = ax_x0, xend = ax_x0 + ax_len_x,
           y = ax_y0, yend = ax_y0,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           linewidth = 0.7, color = "black") +
  annotate("segment",
           x = ax_x0, xend = ax_x0,
           y = ax_y0, yend = ax_y0 + ax_len_y,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           linewidth = 0.7, color = "black") +
  annotate("text",
           x = ax_x0 + ax_len_x / 2,
           y = ax_y0 - yrange * 0.04,
           label = "UMAP_1", fontface = "bold", size = 3.5) +
  annotate("text",
           x = ax_x0 - xrange * 0.04,
           y = ax_y0 + ax_len_y / 2,
           label = "UMAP_2", angle = 90, fontface = "bold", size = 3.5) +

  coord_fixed() +
  theme_void() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin      = margin(15, 15, 15, 20)
  ) +
  labs(color = "Cell Type",
       title = "HR+/HER2- Breast Cancer - Cell Type Annotation")

ggsave(
  file.path(panel_dir, "Annot_03A_UMAP_contour.pdf"),
  plot = pA, width = 6, height = 5.5
)

pB <- ggplot(plot_df, aes(x = umap_1, y = umap_2)) +

  geom_point(aes(color = manual_celltype), size = 0.01, alpha = 0.6) +

  stat_density_2d(
    aes(color = manual_celltype),
    geom        = "density_2d",
    linewidth   = 0.3,
    linetype    = "dashed",
    contour_var = "ndensity",
    breaks      = 0.15
  ) +

  geom_label_repel(
    data          = celltype_pos,
    aes(x = umap_1, y = umap_2,
        label = manual_celltype,
        color = manual_celltype),
    fontface      = "bold",
    size          = 3,
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size  = 0.3,
    fill          = alpha("white", 0.6),
    show.legend   = FALSE,
    max.overlaps  = 30
  ) +

  scale_color_manual(values = clusterCols, guide = "none") +

  annotate("segment",
           x = ax_x0, xend = ax_x0 + ax_len_x,
           y = ax_y0, yend = ax_y0,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           linewidth = 0.7, color = "black") +
  annotate("segment",
           x = ax_x0, xend = ax_x0,
           y = ax_y0, yend = ax_y0 + ax_len_y,
           arrow = arrow(length = unit(0.12, "inches"), type = "closed"),
           linewidth = 0.7, color = "black") +
  annotate("text",
           x = ax_x0 + ax_len_x / 2,
           y = ax_y0 - yrange * 0.04,
           label = "UMAP_1", fontface = "bold", size = 3.5) +
  annotate("text",
           x = ax_x0 - xrange * 0.04,
           y = ax_y0 + ax_len_y / 2,
           label = "UMAP_2", angle = 90, fontface = "bold", size = 3.5) +

  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin     = margin(15, 15, 15, 20)
  ) +
  labs(title = "HR+/HER2- Breast Cancer - Cell Type Annotation")

ggsave(
  file.path(panel_dir, "Annot_03B_UMAP_contour_label.pdf"),
  plot = pB, width = 5.5, height = 5.5
)

saveRDS(
  seurat_obj,
  file = "GSE306201_Step4_the_end_Annotated.rds"
)

annot_table <- data.frame(
  Cluster   = names(cluster_to_type),
  Cell_Type = unname(cluster_to_type),
  row.names = NULL
)
write.csv(
  annot_table,
  file      = file.path(panel_dir, "Annot_00_Cluster_Mapping.csv"),
  row.names = FALSE
)

library(Seurat)
library(tidyverse)
library(patchwork)

setwd(file.path(PROJECT_DIR, "results/GSE306201"))
panel_dir <- "Marker_Panel"
if (!dir.exists(panel_dir)) dir.create(panel_dir, recursive = TRUE)

combined <- readRDS("GSE306201_Step4_the_end_Annotated.rds")

if (!"manual_celltype" %in% colnames(combined@meta.data)) {

  if (!"seurat_clusters" %in% colnames(combined@meta.data)) {
    stop(" manual_celltype,  seurat_clusters, .")
  }

  cluster_to_type <- c(
    "0"  = "Epithelial cells",
    "1"  = "Mesenchymal cells",
    "2"  = "T cells",
    "3"  = "Myeloid cells",
    "4"  = "Epithelial cells",
    "5"  = "Epithelial cells",
    "6"  = "Endothelial cells",
    "7"  = "Mesenchymal cells",
    "8"  = "Mesenchymal cells",
    "9"  = "Epithelial cells",
    "10" = "NK cells",
    "11" = "B cells",
    "12" = "Epithelial cells",
    "13" = "B cells",
    "14" = "B cells",
    "15" = "Epithelial cells",
    "16" = "Epithelial cells",
    "17" = "Mesenchymal cells",
    "18" = "Epithelial cells"
  )

  combined$manual_celltype <- cluster_to_type[as.character(combined$seurat_clusters)]
  combined$manual_celltype <- as.character(combined$manual_celltype)
}

if (all(is.na(combined$manual_celltype))) {
  stop("manual_celltype is entirely NA; check the cluster-to-cell-type mapping.")
}

print(table(combined$manual_celltype, useNA = "ifany"))

cell_order <- c(
  "Epithelial cells",
  "Mesenchymal cells",
  "Myeloid cells",
  "T cells",
  "NK cells",
  "Endothelial cells",
  "B cells"
)

cell_colors <- c(
  "Epithelial cells"  = "#b43665",
  "Mesenchymal cells" = "#6fa6cf",
  "Myeloid cells"     = "#e58027",
  "T cells"           = "#187d79",
  "NK cells"          = "#efb421",
  "Endothelial cells" = "#e70979",
  "B cells"           = "#8560af"
)

present_types <- intersect(cell_order, unique(combined$manual_celltype))
cell_colors <- cell_colors[present_types]
cell_order <- present_types

p2 <- combined@meta.data %>%
  filter(!is.na(manual_celltype)) %>%
  mutate(manual_celltype = factor(manual_celltype, levels = cell_order)) %>%
  ggplot(aes(
    y = forcats::fct_rev(forcats::fct_infreq(manual_celltype)),
    fill = manual_celltype
  )) +
  geom_bar(stat = "count", width = 0.75) +
  labs(x = "Cell count", y = NULL) +
  scale_fill_manual(
    name = "Cell type",
    values = cell_colors,
    breaks = cell_order
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

if (!"orig.ident" %in% colnames(combined@meta.data)) {
  stop(" orig.ident .")
}

cell_counts <- as.data.frame(table(combined$manual_celltype, combined$orig.ident))
colnames(cell_counts) <- c("manual_celltype", "orig.ident", "Freq")
cell_counts <- cell_counts %>%
  filter(!is.na(manual_celltype)) %>%
  mutate(manual_celltype = factor(manual_celltype, levels = cell_order))

p3 <- ggplot(data = cell_counts, aes(
  x = forcats::fct_rev(orig.ident),
  y = Freq,
  fill = manual_celltype
)) +
  geom_bar(position = "fill", stat = "identity", width = 0.75) +
  coord_flip() +
  labs(x = NULL, y = "Cell type frequency") +
  scale_fill_manual(
    name = "Cell type",
    values = cell_colors,
    breaks = cell_order
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

fig5e <- p2 + p3 + plot_layout(guides = "collect") &
  plot_annotation(
    title = "Atlas Composition",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold"
      )
    )
  )

output_file <- file.path(panel_dir, "Fig5E_Atlas_Composition.pdf")

ggsave(
  filename = output_file,
  plot = fig5e,
  width = 12,
  height = 5,
  device = cairo_pdf
)

library(Seurat)
library(irGSEA)
library(ggplot2)
library(paletteer)
library(ggpubr)
library(tidyr)
library(dplyr)

output_dir <- "score"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

scRNA <- readRDS("GSE306201_Step4_the_end_Annotated.rds")
Idents(scRNA) <- scRNA$manual_celltype

rrgene <- read.table(file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt"), header = FALSE, stringsAsFactors = FALSE)[,1]

rrgene <- intersect(rrgene, rownames(scRNA))

rr_list <- list(RR_panel = rrgene)

methods <- c("AUCell", "UCell", "singscore", "AddModuleScore")

res_list <- vector("list", length(methods))
names(res_list) <- methods

rf = file.path(output_dir, "score.Rdata")

if(!file.exists(rf)){
  for (m in methods) {

    res_list[[m]] <- tryCatch(
      irGSEA.score(
        object  = scRNA,
        assay   = "RNA",
        slot    = "data",
        custom  = TRUE,
        geneset = rr_list,
        msigdb  = FALSE,
        method  = m,
        seeds   = 123,
        ncores  = 1,
        minGSSize = 1,
        maxGSSize = 2000
      ),
      error = function(e) {
        NULL
      }
    )

    if(!is.null(res_list[[m]])){
      saveRDS(res_list[[m]], file = file.path(output_dir, paste0("irGSEA_", m, ".rds")))
    }
  }

  rr_all = lapply(1:length(res_list), function(i){

    if(is.null(res_list[[i]])) return(NULL)

    assay_name <- names(res_list[[i]]@assays)[2]
    dat <- GetAssayData(res_list[[i]], assay = assay_name, layer = "data")

    dat_t <- t(as.matrix(dat))

    colnames(dat_t) <- names(res_list)[i]
    return(as.matrix(dat_t))
  })

  rr_all <- rr_all[!sapply(rr_all, is.null)]

  rr_all = do.call(cbind, rr_all)
  save(rr_all, file = rf)
} else {
  load(rf)
}

fai = setdiff(methods,colnames(rr_all))
print(apply(rr_all, 2, range, na.rm = TRUE))

rr_normalized <- apply(rr_all, 2, function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)

  if (is.na(x_min) || x_max == x_min) {
    return(rep(0, length(x)))
  }

  (x - x_min) / (x_max - x_min)
})

Scoring <- rowMeans(rr_normalized, na.rm = TRUE)

drop_cols <- intersect(colnames(scRNA@meta.data), c(colnames(rr_normalized), "Scoring"))
if(length(drop_cols) > 0) scRNA@meta.data <- scRNA@meta.data[, !colnames(scRNA@meta.data) %in% drop_cols]

scRNA@meta.data <- cbind(scRNA@meta.data, rr_normalized, Scoring)

all_scores <- c(colnames(rr_all), "Scoring")

p1 <- DotPlot(scRNA,
              features = all_scores,
              scale = FALSE) +
  RotatedAxis() +
  scale_color_gradientn(colours = paletteer::paletteer_c("grDevices::Zissou 1", 30)) +
  labs(title = "DotPlot (Original Tutorial Style)")

ggsave(file.path(output_dir, "Plot1_DotPlot.pdf"), p1, width = 8, height = 6)

score_cols <- colnames(rr_all)

plot_data <- scRNA@meta.data %>%
  mutate(cell_type = Idents(scRNA)) %>%
  select(all_of(c("cell_type", score_cols))) %>%
  pivot_longer(cols = all_of(score_cols), names_to = "Method", values_to = "Score")

plot_summary <- plot_data %>%
  group_by(cell_type, Method) %>%
  summarise(
    avg_score = mean(Score, na.rm = TRUE),
    pct_expressed = sum(abs(Score) > 0, na.rm = TRUE) / n() * 100,
    .groups = "drop"
  )

plot_summary <- plot_summary %>%
  group_by(Method) %>%
  mutate(
    avg_score_scaled = if (max(avg_score) == min(avg_score)) {
      rep(0, n())
    } else {
      (avg_score - min(avg_score)) / (max(avg_score) - min(avg_score))
    }
  ) %>%
  ungroup()

p2 <- ggplot(plot_summary, aes(x = Method, y = cell_type)) +
  geom_point(aes(size = pct_expressed, color = avg_score_scaled)) +
  scale_size_continuous(range = c(0, 6), name = "Percent") +
  scale_color_paletteer_c("grDevices::Zissou 1", name = "Nom_Score") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  labs(title = "ggplot2 Bubble Plot (High Contrast)")

ggsave(file.path(output_dir, "Plot2_ggplot_Bubble.pdf"), p2, width = 8, height = 6)

gc()

selected_celltypes <- c("Epithelial cells")

work_dir    <- file.path(PROJECT_DIR, "results/GSE306201")
out_dir     <- file.path(work_dir, "RL65")
seurat_file <- file.path(work_dir, "GSE306201_Step4_the_end_Annotated.rds")
gene_file   <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridis)
  library(MASS)
  library(SCP)
})

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

scRNA <- readRDS(seurat_file)

if (!is.null(selected_celltypes)) {
  scRNA <- subset(scRNA, subset = manual_celltype %in% selected_celltypes)
}

rl65_genes     <- readLines(gene_file, warn = FALSE)
rl65_genes     <- unique(trimws(rl65_genes[nchar(trimws(rl65_genes)) > 0]))
rl65_genes_use <- intersect(rl65_genes, rownames(scRNA))
if (length(rl65_genes_use) < 5) stop("gene_set")

scRNA <- AddModuleScore(scRNA,
                        features = list(RL65 = rl65_genes_use),
                        name     = "RL65_Score")

umap_coord           <- as.data.frame(Embeddings(scRNA, "umap"))
colnames(umap_coord) <- c("UMAP_1", "UMAP_2")

keep_cols <- c("Group", "manual_celltype", "seurat_clusters", "RL65_Score1")
meta      <- scRNA@meta.data[, keep_cols, drop = FALSE]

plot_df <- cbind(umap_coord, meta)
plot_df <- plot_df[plot_df$Group %in% c("Primary", "Metastatic"), ]

colnames(plot_df)[colnames(plot_df) == "RL65_Score1"]    <- "RL65_score"
colnames(plot_df)[colnames(plot_df) == "manual_celltype"] <- "Celltype"

plot_df$Group <- factor(plot_df$Group, levels = c("Primary", "Metastatic"))

for (g in c("Primary", "Metastatic")) {
}

x_range   <- range(plot_df$UMAP_1, na.rm = TRUE)
y_range   <- range(plot_df$UMAP_2, na.rm = TRUE)
max_range <- max(diff(x_range), diff(y_range))
x_center  <- mean(x_range)
y_center  <- mean(y_range)
pad       <- max_range * 0.05

xlim_use <- c(x_center - max_range/2 - pad, x_center + max_range/2 + pad)
ylim_use <- c(y_center - max_range/2 - pad, y_center + max_range/2 + pad)

plot_df$RL65_weight <- plot_df$RL65_score - min(plot_df$RL65_score, na.rm = TRUE) + 1e-6

dens_max_vals <- sapply(c("Primary", "Metastatic"), function(g) {
  d  <- plot_df[plot_df$Group == g, ]
  bw <- c(MASS::bandwidth.nrd(d$UMAP_1), MASS::bandwidth.nrd(d$UMAP_2))
  max(MASS::kde2d(d$UMAP_1, d$UMAP_2, h = bw, n = 300)$z, na.rm = TRUE)
})
global_max <- max(dens_max_vals)

galaxyTheme_black <- function(base_size = 26, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(

      axis.line         = element_blank(),
      axis.text.x       = element_text(size = base_size * 0.8, color = "black", lineheight = 0.9),
      axis.text.y       = element_text(size = base_size * 0.8, color = "black", lineheight = 0.9),
      axis.ticks        = element_line(color = "black", linewidth = 0.2),
      axis.title.x      = element_text(size = base_size, color = "black",
                                       margin = margin(10, 0, 0, 0)),
      axis.title.y      = element_text(size = base_size, color = "black", angle = 90,
                                       margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),

      legend.background = element_rect(color = NA, fill = "white"),
      legend.key        = element_rect(color = "grey80", fill = "white"),
      legend.key.size   = unit(1.2, "lines"),
      legend.text       = element_text(size = base_size * 0.8, color = "black"),
      legend.title      = element_text(size = base_size * 0.8, face = "bold",
                                       hjust = 0, color = "black"),
      legend.position   = "right",
      legend.direction  = "vertical",

      panel.background  = element_rect(fill = "black", color = NA),
      panel.border      = element_rect(fill = NA, color = "white", linewidth = 0.8),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),

      strip.background  = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x      = element_text(size = base_size * 0.8, color = "white", face = "bold"),
      strip.text.y      = element_text(size = base_size * 0.8, color = "white", angle = -90),

      plot.background   = element_rect(color = "white", fill = "white"),
      plot.title        = element_text(size = base_size * 1.0, color = "black",
                                       face = "bold", hjust = 0.5),
      plot.margin       = unit(rep(1, 4), "lines")
    )
}

make_galaxy_plot <- function(group_name) {
  d      <- plot_df[plot_df$Group == group_name, ]
  n_cell <- nrow(d)

  ggplot(d, aes(x = UMAP_1, y = UMAP_2)) +
    stat_density_2d(
      aes(fill = after_stat(density), weight = RL65_weight),
      geom    = "raster",
      contour = FALSE,
      n       = 300
    ) +
    geom_point(color = "white", size = 0.01, alpha = 0.5) +
    scale_fill_viridis(
      option = "magma",
      name   = "RL65\ndensity",
      limits = c(0, global_max),
      oob    = scales::squish
    ) +
    scale_x_continuous(limits = xlim_use, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim_use, expand = c(0, 0)) +
    coord_fixed() +
    labs(
      title = paste0(group_name, "\n(n = ", format(n_cell, big.mark = ","), ")"),
      x     = "UMAP_1",
      y     = "UMAP_2"
    ) +
    galaxyTheme_black(base_size = 26)
}

p_g1 <- make_galaxy_plot("Primary")
p_g2 <- make_galaxy_plot("Metastatic")

p_galaxy_final <- p_g1 + p_g2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  file.path(out_dir, "PDF1_RL65_Galaxy_Primary_vs_Metastatic.pdf"),
  plot   = p_galaxy_final,
  width  = 24,
  height = 11,
  bg     = "white"
)

scRNA_p <- subset(scRNA, subset = Group == "Primary")
scRNA_m <- subset(scRNA, subset = Group == "Metastatic")

p_c1 <- CellDimPlot(
  srt              = scRNA_p,
  group.by         = "seurat_clusters",
  reduction        = "umap",
  pt.size          = 0.01,
  label            = FALSE,
  title            = paste0("Primary\n(n = ", format(ncol(scRNA_p), big.mark = ","), ")"),
  legend.position  = "right"
) + xlim(xlim_use) + ylim(ylim_use) + coord_fixed()

p_c2 <- CellDimPlot(
  srt              = scRNA_m,
  group.by         = "seurat_clusters",
  reduction        = "umap",
  pt.size          = 0.01,
  label            = FALSE,
  title            = paste0("Metastatic\n(n = ", format(ncol(scRNA_m), big.mark = ","), ")"),
  legend.position  = "right"
) + xlim(xlim_use) + ylim(ylim_use) + coord_fixed()

p_cluster_final <- p_c1 + p_c2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  file.path(out_dir, "PDF2_Epithelial_Cluster_Primary_vs_Metastatic.pdf"),
  plot   = p_cluster_final,
  width  = 20,
  height = 9,
  bg     = "white"
)

violin_df <- data.frame(
  Group      = plot_df$Group,
  RL65_score = plot_df$RL65_score
)
violin_df$Group <- factor(violin_df$Group, levels = c("Primary", "Metastatic"))

grp_colors <- c("Primary" = "#6fa6cf", "Metastatic" = "#b43665")

wt      <- wilcox.test(RL65_score ~ Group, data = violin_df, exact = FALSE)
pval    <- wt$p.value
p_label <- if (pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))

y_max   <- max(violin_df$RL65_score, na.rm = TRUE)
y_min   <- min(violin_df$RL65_score, na.rm = TRUE)
y_span  <- y_max - y_min
anno_y  <- y_max + y_span * 0.08
tick_h  <- y_span * 0.025

p_violin <- ggplot(violin_df, aes(x = Group, y = RL65_score, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
  geom_jitter(
    aes(color = Group),
    width       = 0.18,
    size        = 0.25,
    alpha       = 0.35,
    show.legend = FALSE
  ) +
  geom_boxplot(
    width         = 0.12,
    outlier.shape = NA,
    fill          = "white",
    alpha         = 0.75,
    color         = "black",
    linewidth     = 0.6
  ) +
  scale_fill_manual(values  = grp_colors) +
  scale_color_manual(values = grp_colors) +

  annotate("segment", x = 1, xend = 2,
           y = anno_y, yend = anno_y,
           color = "black", linewidth = 0.6) +
  annotate("segment", x = 1, xend = 1,
           y = anno_y, yend = anno_y - tick_h,
           color = "black", linewidth = 0.6) +
  annotate("segment", x = 2, xend = 2,
           y = anno_y, yend = anno_y - tick_h,
           color = "black", linewidth = 0.6) +
  annotate("text",
           x = 1.5, y = anno_y + y_span * 0.045,
           label = p_label, size = 5.5, fontface = "bold") +
  labs(
    title = "RL65 Score: Primary vs Metastatic\n(Epithelial cells)",
    x     = NULL,
    y     = "RL65 Score"
  ) +
  coord_cartesian(ylim = c(y_min - y_span * 0.05,
                           anno_y + y_span * 0.15)) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text       = element_text(color = "black", size = 13),
    axis.title.y    = element_text(color = "black", size = 14, face = "bold"),
    axis.line       = element_line(color = "black", linewidth = 0.6),
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "PDF3_RL65_Violin_Primary_vs_Metastatic.pdf"),
  plot   = p_violin,
  width  = 5.5,
  height = 7,
  bg     = "white"
)
