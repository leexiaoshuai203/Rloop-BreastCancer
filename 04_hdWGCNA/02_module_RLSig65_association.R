# hdWGCNA module association with RL-Sig65

PROJECT_DIR <- "."

suppressPackageStartupMessages({
  library(Seurat)
  library(hdWGCNA)
  library(tidyverse)
  library(UCell)
  library(pheatmap)
  library(rmcorr)
})

WORK_DIR <- file.path(PROJECT_DIR, "results/GSE306201/hdWGCNA")

OUT_DIR <- file.path(
  WORK_DIR,
  "04_Module_Trait_Analysis"
)

RL65_FILE <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(WORK_DIR)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

seurat_obj <- readRDS(
  file.path(WORK_DIR, "hdWGCNA_Malignant.rds")
)

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj$C5_cluster_label <- paste0(
  "C",
  as.character(seurat_obj$seurat_clusters)
)

wgcna_name <- GetActiveWGCNAName(seurat_obj)
wgcna_params <- seurat_obj@misc[[wgcna_name]]$wgcna_params

if (
  !is.null(wgcna_params) &&
  !is.null(wgcna_params$networkType)
) {

  if (wgcna_params$networkType != "signed") {
    stop("currentsignednetwork, checkinputRDS")
  }

} else {

  warning(
    "loadnetworkType, ConstructNetwork()networkType='signed'"
  )
}

MEs_h <- tryCatch(
  GetMEs(
    seurat_obj,
    harmonized = TRUE
  ),
  error = function(e) NULL
)

if (
  is.null(MEs_h) ||
  ncol(MEs_h) == 0
) {
  stop("harmonized hME")
}

if (any(!is.finite(as.matrix(MEs_h)))) {
  stop("harmonized hMENA, NaNInf")
}

mods_h <- setdiff(
  colnames(MEs_h),
  "grey"
)

RL65_raw <- read.table(
  RL65_FILE,
  header = FALSE,
  stringsAsFactors = FALSE
)[, 1]

RL65_raw <- unique(
  RL65_raw[
    !is.na(RL65_raw) &
      nzchar(RL65_raw)
  ]
)

RL65_use <- intersect(
  RL65_raw,
  rownames(seurat_obj)
)

if (length(RL65_use) == 0) {
  stop("RL65gene_setexpressiongene")
}

rl65_score <- ScoreSignatures_UCell(
  matrix = GetAssayData(
    seurat_obj,
    assay = "RNA",
    slot = "data"
  ),
  features = list(
    RL65 = RL65_use
  ),
  ncores = 1
)

rl_col <- grep(
  "^RL65",
  colnames(rl65_score),
  value = TRUE
)

if (length(rl_col) == 0) {
  stop("RL65 UCell")
}

rl_col <- rl_col[1]

seurat_obj$RL65 <- rl65_score[
  Cells(seurat_obj),
  rl_col
]

MIN_CELLS_PER_GROUP <- 10

cells_use <- intersect(
  Cells(seurat_obj),
  rownames(MEs_h)
)

cell_level_df <- data.frame(
  cell = cells_use,
  orig.ident = as.character(
    seurat_obj@meta.data[cells_use, "orig.ident"]
  ),
  C5_cluster_label = as.character(
    seurat_obj@meta.data[cells_use, "C5_cluster_label"]
  ),
  RL65 = as.numeric(
    seurat_obj@meta.data[cells_use, "RL65"]
  ),
  stringsAsFactors = FALSE
)

for (m in mods_h) {
  cell_level_df[[m]] <- as.numeric(
    MEs_h[cells_use, m]
  )
}

sample_cluster_agg <- cell_level_df %>%
  dplyr::group_by(
    orig.ident,
    C5_cluster_label
  ) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    RL65_mean = mean(
      RL65,
      na.rm = TRUE
    ),
    dplyr::across(
      dplyr::all_of(mods_h),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_cells >= MIN_CELLS_PER_GROUP
  )

sample_group_number <- table(
  sample_cluster_agg$orig.ident
)

valid_samples <- names(
  sample_group_number[
    sample_group_number >= 2
  ]
)

sample_cluster_agg <- sample_cluster_agg %>%
  dplyr::filter(
    orig.ident %in% valid_samples
  )

if (
  length(unique(sample_cluster_agg$orig.ident)) < 3
) {
  stop("analysissample3")
}

write.csv(
  sample_cluster_agg,
  file.path(
    OUT_DIR,
    "10_samplecluster_hMERL65summary.csv"
  ),
  row.names = FALSE
)

mod_rl65_list <- lapply(
  mods_h,
  function(m) {

    dat_m <- data.frame(
      participant = factor(
        sample_cluster_agg$orig.ident
      ),
      RL65 = as.numeric(
        sample_cluster_agg$RL65_mean
      ),
      hME = as.numeric(
        sample_cluster_agg[[m]]
      )
    )

    dat_m <- dat_m[
      complete.cases(dat_m) &
        is.finite(dat_m$RL65) &
        is.finite(dat_m$hME),
      ,
      drop = FALSE
    ]

    participant_number <- table(
      dat_m$participant
    )

    valid_participants <- names(
      participant_number[
        participant_number >= 2
      ]
    )

    dat_m <- dat_m[
      dat_m$participant %in% valid_participants,
      ,
      drop = FALSE
    ]

    dat_m$participant <- droplevels(
      dat_m$participant
    )

    fit <- rmcorr::rmcorr(
      participant = participant,
      measure1 = RL65,
      measure2 = hME,
      dataset = dat_m,
      CI.level = 0.95
    )

    ci_values <- as.numeric(
      fit$CI
    )

    data.frame(
      module = m,
      cor = as.numeric(fit$r),
      CI_low = ci_values[1],
      CI_high = ci_values[2],
      p_value = as.numeric(fit$p),
      n_units = nrow(dat_m),
      n_samples = nlevels(dat_m$participant),
      stringsAsFactors = FALSE
    )
  }
)

mod_rl65_cor <- dplyr::bind_rows(
  mod_rl65_list
) %>%
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(cor)
  )

write.csv(
  mod_rl65_cor,
  file.path(
    OUT_DIR,
    "11_hMERL65.csv"
  ),
  row.names = FALSE
)

plot_cor_df <- mod_rl65_cor %>%
  dplyr::filter(
    is.finite(cor),
    is.finite(FDR)
  )

plot_cor_df$module <- factor(
  plot_cor_df$module,
  levels = plot_cor_df$module
)

plot_cor_df$FDR_plot <- pmax(
  plot_cor_df$FDR,
  .Machine$double.xmin
)

p_bubble <- ggplot(
  plot_cor_df,
  aes(
    x = module,
    y = "RL65",
    size = -log10(FDR_plot),
    color = cor
  )
) +
  geom_point() +
  scale_color_gradient2(
    low = "#6fa6cf",
    mid = "grey90",
    high = "#b43665",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_classic() +
  labs(
    x = "Module",
    y = NULL,
    size = "-log10(FDR)",
    color = "Repeated-measures\ncorrelation",
    title = paste0(
      "Module hME correlation with RL65 UCell score\n",
      "(sample x cluster level)"
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    OUT_DIR,
    "12A_RL65plot.pdf"
  ),
  plot = p_bubble,
  width = 12,
  height = 4,
  units = "in",
  device = "pdf"
)

heatmap_df <- mod_rl65_cor %>%
  dplyr::filter(
    is.finite(cor)
  )

if (nrow(heatmap_df) == 0) {
  stop("plotanalysisresults")
}

cor_mat_for_heatmap <- matrix(
  heatmap_df$cor,
  ncol = 1,
  dimnames = list(
    heatmap_df$module,
    "RL65"
  )
)

format_probability <- function(x) {

  result <- rep(
    "NA",
    length(x)
  )

  result[
    !is.na(x) & x < 0.001
  ] <- "<0.001"

  result[
    !is.na(x) & x >= 0.001
  ] <- sprintf(
    "%.3f",
    x[!is.na(x) & x >= 0.001]
  )

  result
}

cor_display <- heatmap_df$cor

cor_display[
  abs(cor_display) < 0.005
] <- 0

core_text <- matrix(
  paste0(
    "r=",
    sprintf("%.2f", cor_display),
    "\nP=",
    format_probability(heatmap_df$p_value),
    "\nFDR=",
    format_probability(heatmap_df$FDR)
  ),
  ncol = 1,
  dimnames = list(
    heatmap_df$module,
    "RL65"
  )
)

module_names <- rownames(
  cor_mat_for_heatmap
)

is_valid_color <- function(x) {
  tryCatch(
    {
      col2rgb(x)
      TRUE
    },
    error = function(e) FALSE
  )
}

module_colors_vec <- ifelse(
  vapply(
    module_names,
    is_valid_color,
    logical(1)
  ),
  module_names,
  "grey60"
)

color_annotation <- data.frame(
  Module = factor(
    module_names,
    levels = module_names
  )
)

rownames(color_annotation) <- module_names

module_color_map <- setNames(
  module_colors_vec,
  module_names
)

pdf(
  file.path(
    OUT_DIR,
    "12B_Module_Rloop_Enhanced.pdf"
  ),
  width = 8,
  height = max(
    6,
    nrow(cor_mat_for_heatmap) * 0.60 + 3
  )
)

pheatmap(
  cor_mat_for_heatmap,
  color = colorRampPalette(
    c(
      "#053061",
      "#2166AC",
      "#4393C3",
      "#92C5DE",
      "#D1E5F0",
      "#FFFFFF",
      "#FDDBC7",
      "#F4A582",
      "#D6604D",
      "#B2182B",
      "#67001F"
    )
  )(100),
  breaks = seq(
    -1,
    1,
    length.out = 101
  ),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 12,
  main = paste0(
    "Module-RL65 relationship\n",
    "(repeated-measures correlation)"
  ),
  border_color = "white",
  cellwidth = 130,
  cellheight = 42,
  display_numbers = core_text,
  number_color = "black",
  fontsize_number = 7.5,
  annotation_row = color_annotation,
  annotation_colors = list(
    Module = module_color_map
  ),
  annotation_names_row = FALSE,
  legend = TRUE,
  legend_breaks = c(
    -1,
    -0.5,
    0,
    0.5,
    1
  ),
  legend_labels = c(
    "-1",
    "-0.5",
    "0",
    "0.5",
    "1"
  )
)

dev.off()

positive_modules <- mod_rl65_cor %>%
  dplyr::filter(
    is.finite(cor),
    cor > 0
  ) %>%
  dplyr::arrange(
    dplyr::desc(cor)
  )

significant_positive_modules <- mod_rl65_cor %>%
  dplyr::filter(
    is.finite(cor),
    cor > 0,
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  dplyr::arrange(
    dplyr::desc(cor)
  )

if (nrow(positive_modules) == 0) {

  target_module <- NA_character_
  target_module_df <- positive_modules

} else {

  target_module_df <- positive_modules[
    1,
    ,
    drop = FALSE
  ]

  target_module <- target_module_df$module[1]

  if (
    !is.na(target_module_df$FDR[1]) &&
    target_module_df$FDR[1] < 0.05
  ) {

  } else {

  }
}

write.csv(
  significant_positive_modules,
  file.path(
    OUT_DIR,
    "13_RL65.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    mod_rl65_cor = mod_rl65_cor,
    sample_cluster_agg = sample_cluster_agg,
    target_module = target_module,
    target_module_df = target_module_df,
    significant_positive_modules = significant_positive_modules,
    RL65_genes = RL65_use,
    minimum_cells_per_group = MIN_CELLS_PER_GROUP
  ),
  file.path(
    OUT_DIR,
    "RL65_module_correlation_bundle.rds"
  )
)
