# C5 differential analysis and candidate-gene selection

PROJECT_DIR <- "."

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

BASE_WORK <- file.path(PROJECT_DIR, "results", "GSE306201", "candidate_selection")
dir.create(BASE_WORK, recursive = TRUE, showWarnings = FALSE)
setwd(BASE_WORK)

MAL_RDS <- file.path(PROJECT_DIR, "results", "GSE306201", "Malignant", "Malignant_Epithelial_Reclustered.rds")

OUT_DIR <- "02"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FM_MIN_PCT   <- 0.05
FM_LOGFC_THR <- 0
FM_TEST      <- "wilcox"

MIN_CELLS_PER_CLUSTER <- 10

MAIN_FDR_CUT      <- 0.05
MAIN_LOGFC_CUT    <- 0.25
STABLE_FRAC_UP    <- 0.75
STABLE_MEDIAN_LFC <- 0.25

STABLE_MIN_VALID_DEFAULT <- 6

MAX_SKIPPED_BEFORE_ADJUST <- 2

set.seed(1234)

mal_obj <- readRDS(MAL_RDS)
DefaultAssay(mal_obj) <- "RNA"

if (!"data" %in% Layers(mal_obj[["RNA"]])) {
  mal_obj <- NormalizeData(mal_obj, verbose = FALSE)
}

mal_obj$C5_cluster_label <- paste0("C", mal_obj$seurat_clusters)
mal_obj$C5_state <- ifelse(mal_obj$C5_cluster_label == "C5", "C5", "Other")

print(table(mal_obj$C5_cluster_label))

markers_main <- FindMarkers(
  object          = mal_obj,
  ident.1         = "C5",
  ident.2         = "Other",
  group.by        = "C5_state",
  assay           = "RNA",
  test.use        = FM_TEST,
  min.pct         = FM_MIN_PCT,
  logfc.threshold = FM_LOGFC_THR,
  only.pos        = FALSE,
  verbose         = FALSE
) %>%
  rownames_to_column("gene") %>%
  select(gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
  arrange(p_val_adj)

all_clusters   <- sort(unique(mal_obj$C5_cluster_label))
other_clusters <- setdiff(all_clusters, "C5")
n_candidate_clusters <- length(other_clusters)

pairwise_list <- list()

for (cl in other_clusters) {

  n_ref <- sum(mal_obj$C5_cluster_label == cl)
  if (n_ref < MIN_CELLS_PER_CLUSTER) {
    next
  }

  res <- FindMarkers(
    object          = mal_obj,
    ident.1         = "C5",
    ident.2         = cl,
    group.by        = "C5_cluster_label",
    assay           = "RNA",
    test.use        = FM_TEST,
    min.pct         = FM_MIN_PCT,
    logfc.threshold = FM_LOGFC_THR,
    only.pos        = FALSE,
    verbose         = FALSE
  ) %>%
    rownames_to_column("gene") %>%
    select(gene, avg_log2FC) %>%
    mutate(ref_cluster = cl)

  pairwise_list[[cl]] <- res
}

pairwise_df <- bind_rows(pairwise_list)
n_valid_pairwise_comparisons <- length(pairwise_list)
n_skipped_clusters <- n_candidate_clusters - n_valid_pairwise_comparisons

if (n_skipped_clusters > MAX_SKIPPED_BEFORE_ADJUST) {
  STABLE_MIN_VALID <- ceiling(0.75 * n_valid_pairwise_comparisons)
} else {
  STABLE_MIN_VALID <- STABLE_MIN_VALID_DEFAULT
}

stability_stats <- pairwise_df %>%
  group_by(gene) %>%
  summarise(
    n_valid_comparisons    = n(),
    n_up                   = sum(avg_log2FC > 0, na.rm = TRUE),
    frac_up                = n_up / n_valid_comparisons,
    median_log2FC_pairwise = median(avg_log2FC, na.rm = TRUE),
    min_log2FC_pairwise    = min(avg_log2FC, na.rm = TRUE),
    .groups = "drop"
  )

markers_summary <- markers_main %>%
  left_join(stability_stats, by = "gene") %>%
  mutate(
    n_valid_comparisons = replace_na(n_valid_comparisons, 0),
    is_stable_C5_up = (p_val_adj < MAIN_FDR_CUT) &
      (avg_log2FC >= MAIN_LOGFC_CUT) &
      (n_valid_comparisons >= STABLE_MIN_VALID) &
      (frac_up >= STABLE_FRAC_UP) &
      (median_log2FC_pairwise >= STABLE_MEDIAN_LFC)
  ) %>%
  arrange(desc(is_stable_C5_up), p_val_adj)

n_stable <- sum(markers_summary$is_stable_C5_up, na.rm = TRUE)

write.csv(
  markers_summary,
  file.path(OUT_DIR, "C5_FindMarkerssummaryresults.csv"),
  row.names = FALSE
)

stable_up_genes <- markers_summary %>%
  filter(is_stable_C5_up) %>%
  pull(gene)

writeLines(
  stable_up_genes,
  file.path(OUT_DIR, "C5_stable_up_genes.txt")
)

colnames(mal_obj@meta.data)

if ("package:plyr" %in% search()) {
}

sample_counts <- mal_obj@meta.data %>%
  dplyr::group_by(.data[["orig.ident"]], .data[["C5_state"]]) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = "C5_state",
    values_from = "n_cells",
    values_fill = 0
  ) %>%
  dplyr::rename(sample_id = "orig.ident") %>%
  dplyr::mutate(
    n_total     = C5 + Other,
    C5_fraction = round(C5 / n_total, 4)
  ) %>%
  dplyr::arrange(dplyr::desc(C5))

print(sample_counts)
OTHER_MIN <- 50
thresholds <- c(10, 20, 30)

sensitivity_df <- lapply(thresholds, function(th) {
  passed <- sample_counts %>% filter(C5 >= th, Other >= OTHER_MIN)
  data.frame(C5_threshold = th, Other_threshold = OTHER_MIN,
             n_samples_passed = nrow(passed))
}) %>% bind_rows()

n_at_20 <- sensitivity_df$n_samples_passed[sensitivity_df$C5_threshold == 20]
recommended_threshold <- if (n_at_20 >= 6) 20 else 10

result_bundle <- list(
  markers_summary          = markers_summary,
  stable_up_genes           = stable_up_genes,
  n_valid_pairwise_comparisons = n_valid_pairwise_comparisons,
  STABLE_MIN_VALID_used    = STABLE_MIN_VALID,
  sample_counts            = sample_counts,
  sensitivity_df            = sensitivity_df,
  recommended_C5_threshold  = recommended_threshold
)

saveRDS(result_bundle, file.path(OUT_DIR, "C5_DE_results_bundle.rds"))

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(Matrix)
})

BASE_DIR <- file.path(PROJECT_DIR, "results/GSE306201/candidate_selection")
setwd(BASE_DIR)

MAL_RDS <- file.path(PROJECT_DIR, "results", "GSE306201", "Malignant", "Malignant_Epithelial_Reclustered.rds")
PREV_BUNDLE <- file.path(BASE_DIR, "02/C5_DE_results_bundle.rds")

OUT_DIR <- file.path(BASE_DIR, "03_C5_Pseudobulk")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

if (!file.exists(MAL_RDS))     stop(" MAL_RDS: ", MAL_RDS)
if (!file.exists(PREV_BUNDLE)) stop(" PREV_BUNDLE: ", PREV_BUNDLE)

OTHER_MIN_MAIN <- 50

set.seed(1234)

mal_obj <- readRDS(MAL_RDS)
DefaultAssay(mal_obj) <- "RNA"

mal_obj$C5_cluster_label <- paste0("C", mal_obj$seurat_clusters)
mal_obj$C5_state <- ifelse(mal_obj$C5_cluster_label == "C5", "C5", "Other")

prev_bundle <- readRDS(PREV_BUNDLE)
sample_counts_prev    <- prev_bundle$sample_counts
recommended_threshold <- prev_bundle$recommended_C5_threshold

safe_flatten <- function(x) {
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  as.character(x)
}

orig_ident_clean <- safe_flatten(mal_obj@meta.data$orig.ident)
group_clean       <- safe_flatten(mal_obj@meta.data$Group)

stopifnot(length(orig_ident_clean) == nrow(mal_obj@meta.data))
stopifnot(length(group_clean)      == nrow(mal_obj@meta.data))

raw_map <- unique(data.frame(
  orig.ident      = orig_ident_clean,
  tissue_type_raw = group_clean,
  stringsAsFactors = FALSE
))

id_counts <- table(raw_map$orig.ident)
if (any(id_counts > 1)) {
  stop("sampleGroup, check: ",
       paste(names(id_counts[id_counts > 1]), collapse = ", "))
}

raw_map$tissue_type <- ifelse(
  grepl("prim", tolower(raw_map$tissue_type_raw)), "Primary",
  ifelse(grepl("meta", tolower(raw_map$tissue_type_raw)), "Metastatic",
         raw_map$tissue_type_raw)
)

raw_map$patient_id <- sub("_[^_]+$", "", raw_map$orig.ident)

sample_info <- raw_map[order(raw_map$orig.ident),
                       c("orig.ident", "tissue_type", "patient_id")]

write.csv(
  sample_info,
  file.path(OUT_DIR, "01_samplegroup.csv"),
  row.names = FALSE
)

print(table(sample_info$tissue_type))
print(sample_info)

counts_mat <- LayerData(mal_obj, assay = "RNA", layer = "counts")

group_vec    <- paste(orig_ident_clean, mal_obj$C5_state, sep = "__")
group_factor <- factor(group_vec)
indicator_mat <- Matrix::sparse.model.matrix(~ 0 + group_factor)
colnames(indicator_mat) <- levels(group_factor)

pseudobulk_mat <- as.matrix(counts_mat %*% indicator_mat)
rownames(pseudobulk_mat) <- rownames(counts_mat)

included_samples_main <- sample_counts_prev %>%
  filter(C5 >= recommended_threshold, Other >= OTHER_MIN_MAIN) %>%
  left_join(sample_info, by = c("sample_id" = "orig.ident")) %>%
  filter(!is.na(tissue_type))

print(included_samples_main %>% select(sample_id, C5, Other, tissue_type, patient_id))

write.csv(
  included_samples_main,
  file.path(OUT_DIR, "02_bulksample.csv"),
  row.names = FALSE
)

if (nrow(included_samples_main) < 4) {
  stop("analysissample(<4), analysis, check.")
}

patient_counts_main <- table(included_samples_main$patient_id)
multi_patients <- names(patient_counts_main[patient_counts_main > 1])
n_multi <- length(multi_patients)

print(patient_counts_main)

if (n_multi > 0) {
} else {
}

run_paired_edgeR <- function(pseudobulk_mat, sample_ids) {

  meta_sub <- tidyr::expand_grid(
    sample_id = as.character(sample_ids),
    C5_state = c("Other", "C5")
  ) %>%
    dplyr::mutate(
      colname = paste(sample_id, C5_state, sep = "__"),
      sample_id = factor(sample_id, levels = sample_ids),
      C5_state = factor(C5_state, levels = c("Other", "C5"))
    )

  missing_cols <- setdiff(meta_sub$colname, colnames(pseudobulk_mat))
  if (length(missing_cols) > 0) {
    stop("sample-bulkmissing: ",
         paste(missing_cols, collapse = ", "))
  }

  mat_sub <- pseudobulk_mat[, meta_sub$colname, drop = FALSE]

  y <- edgeR::DGEList(counts = mat_sub)

  keep_genes <- edgeR::filterByExpr(y, group = meta_sub$C5_state)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]

  y <- edgeR::calcNormFactors(y)

  design <- model.matrix(~ sample_id + C5_state, data = meta_sub)

  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)

  if (!"C5_stateC5" %in% colnames(design)) {
    stop("not_foundC5_stateC5")
  }

  qlf <- edgeR::glmQLFTest(fit, coef = "C5_stateC5")

  edgeR::topTags(qlf, n = Inf)$table %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::rename(log2FC = logFC, PValue = PValue, FDR = FDR) %>%
    dplyr::arrange(FDR)
}

main_sample_ids <- included_samples_main$sample_id

de_main <- run_paired_edgeR(pseudobulk_mat, main_sample_ids)

if (!is.data.frame(de_main) || !all(c("gene", "log2FC", "FDR") %in% colnames(de_main))) {
  stop("de_main(gene/log2FC/FDRdata.frame).\n",
       "current: \n",
       paste(capture.output(str(de_main, max.level = 1)), collapse = "\n"))
}

write.csv(
  de_main,
  file.path(OUT_DIR, "03_analysisbulk_edgeRresults.csv"),
  row.names = FALSE
)

c5_up_genes_main <- de_main %>%
  dplyr::filter(FDR < 0.05, log2FC > 0.25) %>%
  dplyr::pull(gene)

writeLines(
  c5_up_genes_main,
  file.path(OUT_DIR, "04_analysisC5gene.txt")
)

pseudobulk_bundle <- list(
  pseudobulk_mat        = pseudobulk_mat,
  sample_info           = sample_info,
  included_samples_main = included_samples_main,
  de_main               = de_main,
  c5_up_genes           = c5_up_genes_main,
  recommended_threshold = recommended_threshold
)

saveRDS(pseudobulk_bundle, file.path(OUT_DIR, "C5_Pseudobulk_results_bundle.rds"))

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(tidyverse)
  library(VennDiagram)
  library(UpSetR)
  library(grid)
})

BASE_DIR <- file.path(PROJECT_DIR, "results/GSE306201/candidate_selection")
OUT_DIR  <- file.path(BASE_DIR, "05_")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

rl65_cor_file <- file.path(
  BASE_DIR,
  "04_Malignant_hdWGCNA_mainline/04_Module_Trait_Analysis/11_hMERL65.csv"
)

if (!file.exists(rl65_cor_file)) {
  stop(": ", rl65_cor_file)
}

rl65_cor <- read.csv(rl65_cor_file, stringsAsFactors = FALSE)

required_cols <- c("module", "cor", "p_value")
if (!all(required_cols %in% colnames(rl65_cor))) {
  stop("resultsmissing: module / cor / p_value")
}

selected_modules <- rl65_cor %>%
  filter(!is.na(cor), !is.na(p_value), cor > 0.30, p_value < 0.05) %>%
  arrange(desc(cor)) %>%
  pull(module) %>%
  unique()

if (length(selected_modules) == 0) {
  stop(" cor > 0.3  P < 0.05 ")
}

hdwgcna_rds <- file.path(
  BASE_DIR,
  "04_Malignant_hdWGCNA_mainline/hdWGCNA_Malignant.rds"
)

if (!file.exists(hdwgcna_rds)) {
  stop(": ", hdwgcna_rds)
}

seurat_obj <- readRDS(hdwgcna_rds)
module_table <- GetModules(seurat_obj)

if (!all(c("module", "gene_name") %in% colnames(module_table))) {
  stop("GetModules()resultsmissingmodulegene_name")
}

module_gene_list <- lapply(selected_modules, function(mod) {
  unique(module_table$gene_name[module_table$module == mod])
})
names(module_gene_list) <- selected_modules

for (mod in selected_modules) {
}

selected_module_union <- unique(unlist(module_gene_list))

stable_fm_file <- file.path(BASE_DIR, "02/C5_stable_up_genes.txt")
if (!file.exists(stable_fm_file)) {
  stop(": ", stable_fm_file)
}

stable_fm_genes <- readLines(stable_fm_file)
stable_fm_genes <- unique(stable_fm_genes[nzchar(stable_fm_genes)])

pb_file <- file.path(BASE_DIR, "03_C5_Pseudobulk/04_analysisC5gene.txt")
if (!file.exists(pb_file)) {
  stop(": ", pb_file)
}

pb_up_genes <- readLines(pb_file)
pb_up_genes <- unique(pb_up_genes[nzchar(pb_up_genes)])

venn_sets <- list(
  RL65_modules_union = selected_module_union,
  C5_FindMarkers     = stable_fm_genes,
  C5_Pseudobulk      = pb_up_genes
)

intersect_genes <- Reduce(intersect, venn_sets)
intersect_genes <- sort(unique(intersect_genes))

if (length(intersect_genes) > 0) {
} else {
}

writeLines(
  intersect_genes,
  file.path(OUT_DIR, "gene.txt")
)

venn_grob <- venn.diagram(
  x = venn_sets,
  category.names = c(
    paste0("RL65-related modules\nunion genes\nn=", length(venn_sets$RL65_modules_union)),
    paste0("C5 FindMarkers\nn=", length(venn_sets$C5_FindMarkers)),
    paste0("C5 Pseudobulk\nn=", length(venn_sets$C5_Pseudobulk))
  ),
  filename = NULL,
  fill = c("#b43665", "#6fa6cf", "#e58027"),
  col = c("#b43665", "#6fa6cf", "#e58027"),
  alpha = 0.45,
  lwd = 2,
  cex = 1.4,
  fontface = "bold",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.col = c("#b43665", "#6fa6cf", "#e58027"),
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.06, 0.06, 0.05),
  margin = 0.08
)

pdf(
  file.path(OUT_DIR, "plot.pdf"),
  width = 8,
  height = 8
)
grid.newpage()
grid.draw(venn_grob)
dev.off()

upset_sets <- c(
  module_gene_list,
  list(
    C5_FindMarkers = stable_fm_genes,
    C5_Pseudobulk  = pb_up_genes
  )
)

upset_input <- fromList(upset_sets)

pdf(
  file.path(OUT_DIR, "UpSetplot.pdf"),
  width = 12,
  height = 7
)

upset(
  upset_input,
  nsets = length(upset_sets),
  nintersects = 30,
  order.by = c("freq", "degree"),
  decreasing = c(TRUE, FALSE),
  sets.bar.color = c(
    rep("#8172B2", length(module_gene_list)),
    "#4C72B0",
    "#55A868"
  ),
  main.bar.color = "#C44E52",
  matrix.color = "#333333",
  mainbar.y.label = "Intersection size",
  sets.x.label = "Set size",
  text.scale = c(1.4, 1.2, 1.2, 1, 1.2, 1.2)
)

dev.off()

suppressPackageStartupMessages({
  library(GSEABase)
})

CANDIDATE_FILE <- file.path(
  OUT_DIR,
  "gene.txt"
)

RLOOP_GMT_FILE <- file.path(PROJECT_DIR, "results/bulk/GSE96058/01_Rloop_regulators.gmt")

RL65_FILE <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")

RESULT_FILE <- file.path(
  OUT_DIR,
  "24geneRloopRL65results.txt"
)

required_files <- c(
  CANDIDATE_FILE,
  RLOOP_GMT_FILE,
  RL65_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    ": \n",
    paste(missing_files, collapse = "\n")
  )
}

clean_genes <- function(x) {

  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)

  unique(
    x[
      !is.na(x) &
        nzchar(x)
    ]
  )
}

candidate_genes <- clean_genes(
  readLines(
    CANDIDATE_FILE,
    warn = FALSE
  )
)

if (length(candidate_genes) == 0) {
  stop("candidategeneloadgene")
}

rloop_gsc <- GSEABase::getGmt(
  RLOOP_GMT_FILE,
  geneIdType = GSEABase::SymbolIdentifier()
)

rloop_gene_list <- lapply(
  seq_along(rloop_gsc),
  function(i) {
    GSEABase::geneIds(
      rloop_gsc[[i]]
    )
  }
)

rloop_genes <- clean_genes(
  unlist(
    rloop_gene_list,
    use.names = FALSE
  )
)

rl65_genes <- clean_genes(
  readLines(
    RL65_FILE,
    warn = FALSE
  )
)

candidate_rloop_intersection <- sort(
  intersect(
    candidate_genes,
    rloop_genes
  )
)

candidate_rl65_intersection <- sort(
  intersect(
    candidate_genes,
    rl65_genes
  )
)

candidate_rloop_rl65_intersection <- sort(
  Reduce(
    intersect,
    list(
      candidate_genes,
      rloop_genes,
      rl65_genes
    )
  )
)

genes_or_none <- function(x) {

  if (length(x) == 0) {
    return("")
  }

  x
}

output_text <- c(
  "24candidategeneR-loopRL65gene_setresults",
  "==================================================",
  "",
  paste0(
    "24candidategene: ",
    length(candidate_genes)
  ),
  paste0(
    "R-loopgene_setgene: ",
    length(rloop_genes)
  ),
  paste0(
    "RL65gene: ",
    length(rl65_genes)
  ),
  "",
  "--------------------------------------------------",
  paste0(
    "[1] 24candidategene  intersection  R-loopgene_set, n=",
    length(candidate_rloop_intersection)
  ),
  "--------------------------------------------------",
  genes_or_none(candidate_rloop_intersection),
  "",
  "--------------------------------------------------",
  paste0(
    "[2] 24candidategene  intersection  RL65, n=",
    length(candidate_rl65_intersection)
  ),
  "--------------------------------------------------",
  genes_or_none(candidate_rl65_intersection),
  "",
  "--------------------------------------------------",
  paste0(
    "[3] 24candidategene  intersection  R-loopgene_set  intersection  RL65, n=",
    length(candidate_rloop_rl65_intersection)
  ),
  "--------------------------------------------------",
  genes_or_none(candidate_rloop_rl65_intersection)
)

output_connection <- file(
  RESULT_FILE,
  open = "w",
  encoding = "UTF-8"
)

writeLines(
  output_text,
  con = output_connection
)

close(output_connection)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

BASE_DIR <- file.path(PROJECT_DIR, "results/GSE306201/candidate_selection")

cor_file <- file.path(
  BASE_DIR,
  "04_Malignant_hdWGCNA_mainline/04_Module_Trait_Analysis/11_hMERL65.csv"
)

hd_file <- file.path(
  BASE_DIR,
  "04_Malignant_hdWGCNA_mainline/hdWGCNA_Malignant.rds"
)

fm_file <- file.path(
  BASE_DIR,
  "02/C5_stable_up_genes.txt"
)

pb_file <- file.path(
  BASE_DIR,
  "03_C5_Pseudobulk/04_analysisC5gene.txt"
)

pb_bundle_file <- file.path(
  BASE_DIR,
  "03_C5_Pseudobulk/C5_Pseudobulk_results_bundle.rds"
)

candidate_file <- file.path(
  BASE_DIR,
  "05_/gene.txt"
)

annotation_file <- file.path(
  BASE_DIR,
  "05_/24geneRloopRL65results.txt"
)

need_files <- c(
  cor_file, hd_file, fm_file,
  pb_file, pb_bundle_file, candidate_file
)

if (any(!file.exists(need_files))) {
  stop(
    ": \n",
    paste(need_files[!file.exists(need_files)], collapse = "\n")
  )
}

cor_df <- read.csv(cor_file)

selected <- cor_df %>%
  filter(
    !is.na(cor),
    !is.na(p_value),
    cor > 0.30,
    p_value < 0.05
  ) %>%
  arrange(desc(cor))

if (nrow(selected) == 0) {
  stop(" r > 0.30  P < 0.05 ")
}

obj <- readRDS(hd_file)
module_table <- GetModules(obj)

module_union <- unique(
  module_table$gene_name[
    module_table$module %in% selected$module
  ]
)

module_text <- paste(
  sprintf("%s  r=%.2f", selected$module, selected$cor),
  collapse = "\n"
)

fm_genes <- unique(readLines(fm_file))
fm_genes <- fm_genes[nzchar(fm_genes)]

pb_genes <- unique(readLines(pb_file))
pb_genes <- pb_genes[nzchar(pb_genes)]

pb_bundle <- readRDS(pb_bundle_file)

c5_threshold <- pb_bundle$recommended_threshold
n_samples <- nrow(pb_bundle$included_samples_main)

candidate_genes <- sort(unique(readLines(candidate_file)))
candidate_genes <- candidate_genes[nzchar(candidate_genes)]

candidate_check <- Reduce(
  intersect,
  list(module_union, fm_genes, pb_genes)
)

if (!setequal(candidate_genes, candidate_check)) {
  warning("gene.txt, check")
}

gene_text <- paste(
  strwrap(
    paste(candidate_genes, collapse = "   -   "),
    width = 88
  ),
  collapse = "\n"
)

rloop_n <- rl65_n <- NA_integer_

if (file.exists(annotation_file)) {
  ann <- readLines(annotation_file, encoding = "UTF-8")

  x1 <- grep("^\\[1\\].*n=", ann)
  x2 <- grep("^\\[2\\].*n=", ann)

  if (length(x1) > 0)
    rloop_n <- as.integer(sub(".*n=([0-9]+).*", "\\1", ann[x1[1]]))

  if (length(x2) > 0)
    rl65_n <- as.integer(sub(".*n=([0-9]+).*", "\\1", ann[x2[1]]))
}

annotation_text <- if (!is.na(rloop_n) && !is.na(rl65_n)) {
  paste0(
    "Post-selection annotation     ",
    "R-loop gene set: ", rloop_n, "/", length(candidate_genes),
    "     Direct RL65 overlap: ", rl65_n, "/", length(candidate_genes)
  )
} else {
  "Post-selection annotation"
}

arrow_style <- arrow(
  length = unit(0.13, "cm"),
  type = "closed"
)

p <- ggplot() +

  annotate(
    "text", x = 8.3, y = 11.65,
    label = "Identification of C5-associated R-loop candidate genes",
    family = "Arial", fontface = "bold", size = 3.35
  ) +

  annotate(
    "rect",
    xmin = 1.0, xmax = 15.6,
    ymin = 10.35, ymax = 11.15,
    fill = "#F3F3F3", color = "#555555", linewidth = 0.45
  ) +
  annotate(
    "text", x = 8.3, y = 10.75,
    label = "Malignant epithelial cells   C0-C8   |   C5-focused analysis",
    family = "Arial", fontface = "bold", size = 2.65
  ) +

  annotate(
    "segment", x = 3.0, xend = 3.0,
    y = 10.35, yend = 9.72,
    linewidth = 0.45, color = "#555555",
    arrow = arrow_style
  ) +
  annotate(
    "segment", x = 8.3, xend = 8.3,
    y = 10.35, yend = 9.72,
    linewidth = 0.45, color = "#555555",
    arrow = arrow_style
  ) +
  annotate(
    "segment", x = 13.6, xend = 13.6,
    y = 10.35, yend = 9.72,
    linewidth = 0.45, color = "#555555",
    arrow = arrow_style
  ) +

  annotate(
    "rect",
    xmin = 0.55, xmax = 5.45,
    ymin = 6.25, ymax = 9.65,
    fill = "#F8E9EE", color = "#B45A78", linewidth = 0.55
  ) +
  annotate(
    "text", x = 3.0, y = 9.25,
    label = "hdWGCNA",
    family = "Arial", fontface = "bold", size = 2.9
  ) +
  annotate(
    "text", x = 3.0, y = 8.72,
    label = "Signed co-expression network\nSample x cluster metacells\nHarmonized module eigengenes",
    family = "Arial", size = 2.12, lineheight = 1.12
  ) +
  annotate(
    "text", x = 3.0, y = 7.83,
    label = "Module-RL65 association\nRepeated-measures correlation\nSelection: r > 0.30, P < 0.05",
    family = "Arial", size = 2.12, lineheight = 1.12
  ) +
  annotate(
    "text", x = 3.0, y = 6.93,
    label = module_text,
    family = "Arial", fontface = "bold",
    size = 2.1, lineheight = 1.08
  ) +
  annotate(
    "text", x = 3.0, y = 6.48,
    label = paste0("Selected-module union   n = ", length(module_union)),
    family = "Arial", fontface = "bold", size = 2.2
  ) +

  annotate(
    "rect",
    xmin = 5.85, xmax = 10.75,
    ymin = 6.25, ymax = 9.65,
    fill = "#E9F1F8", color = "#5A86AD", linewidth = 0.55
  ) +
  annotate(
    "text", x = 8.3, y = 9.25,
    label = "C5 stable FindMarkers",
    family = "Arial", fontface = "bold", size = 2.9
  ) +
  annotate(
    "text", x = 8.3, y = 8.58,
    label = "C5 vs all other malignant cells\nWilcoxon test\nFDR < 0.05   |   avg log2FC >= 0.25",
    family = "Arial", size = 2.12, lineheight = 1.12
  ) +
  annotate(
    "text", x = 8.3, y = 7.55,
    label = "Pairwise direction-stability filter\nUpregulated in >=75% comparisons\nMedian pairwise log2FC >= 0.25",
    family = "Arial", size = 2.12, lineheight = 1.12
  ) +
  annotate(
    "text", x = 8.3, y = 6.55,
    label = paste0("Stable C5-up genes   n = ", length(fm_genes)),
    family = "Arial", fontface = "bold", size = 2.25
  ) +

  annotate(
    "rect",
    xmin = 11.15, xmax = 16.05,
    ymin = 6.25, ymax = 9.65,
    fill = "#F9EFE2", color = "#C78338", linewidth = 0.55
  ) +
  annotate(
    "text", x = 13.6, y = 9.25,
    label = "C5 pseudobulk",
    family = "Arial", fontface = "bold", size = 2.9
  ) +
  annotate(
    "text", x = 13.6, y = 8.57,
    label = "Sample x C5-state aggregation\nPaired edgeR quasi-likelihood model\nC5 vs Other within each sample",
    family = "Arial", size = 2.12, lineheight = 1.12
  ) +
  annotate(
    "text", x = 13.6, y = 7.60,
    label = paste0(
      "C5 >= ", c5_threshold,
      " cells   |   Other >= 50 cells\n",
      "Included samples   n = ", n_samples,
      "\nFDR < 0.05   |   log2FC > 0.25"
    ),
    family = "Arial", size = 2.12, lineheight = 1.12
  ) +
  annotate(
    "text", x = 13.6, y = 6.55,
    label = paste0("Pseudobulk C5-up genes   n = ", length(pb_genes)),
    family = "Arial", fontface = "bold", size = 2.25
  ) +

  annotate(
    "segment", x = 3.0, xend = 7.25,
    y = 6.25, yend = 5.28,
    linewidth = 0.5, color = "#666666",
    arrow = arrow_style
  ) +
  annotate(
    "segment", x = 8.3, xend = 8.3,
    y = 6.25, yend = 5.28,
    linewidth = 0.5, color = "#666666",
    arrow = arrow_style
  ) +
  annotate(
    "segment", x = 13.6, xend = 9.35,
    y = 6.25, yend = 5.28,
    linewidth = 0.5, color = "#666666",
    arrow = arrow_style
  ) +

  annotate(
    "rect",
    xmin = 5.35, xmax = 11.25,
    ymin = 4.35, ymax = 5.25,
    fill = "#EEE9F5", color = "#7865A5", linewidth = 0.65
  ) +
  annotate(
    "text", x = 8.3, y = 4.80,
    label = paste0(
      "Three-way intersection     n = ",
      length(candidate_genes)
    ),
    family = "Arial", fontface = "bold", size = 2.75
  ) +

  annotate(
    "segment", x = 8.3, xend = 8.3,
    y = 4.35, yend = 3.82,
    linewidth = 0.55, color = "#666666",
    arrow = arrow_style
  ) +

  annotate(
    "rect",
    xmin = 1.0, xmax = 15.6,
    ymin = 1.25, ymax = 3.75,
    fill = "#F6F3FA", color = "#66508E", linewidth = 0.7
  ) +
  annotate(
    "text", x = 8.3, y = 3.35,
    label = paste0(
      "Final candidate genes   n = ",
      length(candidate_genes)
    ),
    family = "Arial", fontface = "bold", size = 3.0
  ) +
  annotate(
    "text", x = 8.3, y = 2.43,
    label = gene_text,
    family = "Arial", fontface = "bold",
    size = 2.15, lineheight = 1.15
  ) +

  annotate(
    "text", x = 8.3, y = 0.68,
    label = annotation_text,
    family = "Arial", size = 2.15
  ) +

  coord_fixed(
    xlim = c(0, 16.6),
    ylim = c(0.25, 12),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

ggsave(
  filename = file.path(
    BASE_DIR,
    "05_C5candidategeneplot.pdf"
  ),
  plot = p,
  width = 16.6,
  height = 12,
  units = "cm",
  device = cairo_pdf
)
