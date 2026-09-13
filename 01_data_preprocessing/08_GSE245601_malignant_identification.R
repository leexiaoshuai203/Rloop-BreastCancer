# GSE245601 malignant epithelial-cell identification

PROJECT_DIR <- "."

samples <- unique(as.character(epi_seu[[sample_col, drop = TRUE]]))
samples <- samples[!is.na(samples)]

score_list <- list()
ref_list <- list()
threshold_list <- list()

for (sid in samples) {
  epi_cells <- colnames(epi_seu)[
    as.character(epi_seu[[sample_col, drop = TRUE]]) == sid
  ]
  if (length(epi_cells) < 10) {
    next
  }

  epi_sub <- subset(epi_seu, features = common_genes, cells = epi_cells)
  epi_sub$fastcnv_ref_type <- "Epithelial cells"
  epi_sub <- RenameCells(epi_sub,
                         new.names = paste0("EPI_", colnames(epi_sub)))
  ref_sub <- RenameCells(ref_obj,
                         new.names = paste0("REF_", colnames(ref_obj)))

  cnv_input <- merge(epi_sub, y = ref_sub, merge.data = FALSE)
  if (length(Layers(cnv_input[["RNA"]])) > 1) {
    cnv_input <- JoinLayers(cnv_input, assay = "RNA")
  }

  sample_dir <- file.path(run_dir, sid)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(sample_dir)

  res <- tryCatch(
    fastCNV(
      seuratObj = cnv_input,
      sampleName = sid,
      referenceVar = "fastcnv_ref_type",
      referenceLabel = ref_types,
      assay = "RNA",
      reClusterSeurat = FALSE,
      getCNVPerChromosomeArm = FALSE,
      getCNVClusters = FALSE,
      doPlot = TRUE,
      outputType = "pdf"
    ),
    error = function(e) {
      NULL
    }
  )
  setwd(old_wd)
  if (is.null(res)) next

  if (inherits(res, "Seurat")) {
    cnv_obj <- res
  } else if (is.list(res)) {
    idx_seu <- which(vapply(res, inherits, logical(1), "Seurat"))
    if (!length(idx_seu)) {
      next
    }
    cnv_obj <- res[[idx_seu[1]]]
  } else {
    next
  }

  score_col <- grep("cnv_fraction|cnv.*score|score.*cnv",
                    colnames(cnv_obj[[]]),
                    ignore.case = TRUE, value = TRUE)
  if (!length(score_col)) {
    print(colnames(cnv_obj[[]]))
    next
  }
  score_col <- score_col[1]

  meta <- cnv_obj[[]]
  meta_epi <- meta[grepl("EPI_", rownames(meta), fixed = TRUE), , drop = FALSE]
  meta_ref <- meta[grepl("REF_", rownames(meta), fixed = TRUE), , drop = FALSE]
  meta_epi$Cell <- sub("^.*EPI_", "", rownames(meta_epi))
  meta_ref$Cell <- sub("^.*REF_", "", rownames(meta_ref))
  meta_ref$ReferenceType <- ref_map[meta_ref$Cell]

  ref_score <- as.numeric(meta_ref[[score_col]])
  ref_score <- ref_score[is.finite(ref_score)]
  if (length(ref_score) < 50) {
    next
  }

  threshold <- as.numeric(quantile(ref_score, threshold_q))
  epi_score <- as.numeric(meta_epi[[score_col]])
  epi_label <- ifelse(is.finite(epi_score) & epi_score > threshold,
                      "Tumor", "Normal")

  score_list[[sid]] <- data.frame(
    Cell = meta_epi$Cell, Sample = sid,
    FastCNV_score = epi_score, Threshold = threshold,
    Relative_score = epi_score - threshold,
    FastCNV_label = epi_label
  )

  ref_list[[sid]] <- data.frame(
    Sample = sid, ReferenceType = meta_ref$ReferenceType,
    Relative_score = as.numeric(meta_ref[[score_col]]) - threshold
  ) %>% filter(!is.na(ReferenceType), is.finite(Relative_score))

  threshold_list[[sid]] <- data.frame(
    Sample = sid, ReferenceN = length(ref_score),
    Threshold = threshold, EpithelialN = nrow(meta_epi),
    TumorN = sum(epi_label == "Tumor"),
    NormalN = sum(epi_label == "Normal")
  )

  rm(epi_sub, ref_sub, cnv_input, res, cnv_obj, meta)
  invisible(gc())
}

if (!length(score_list)) stop("samplefastCNV")

score_df <- bind_rows(score_list)
ref_df <- bind_rows(ref_list)
threshold_df <- bind_rows(threshold_list)

write.csv(threshold_df,
          file.path(out_dir, "01_sampleFastCNV.csv"),
          row.names = FALSE)
write.csv(score_df,
          file.path(out_dir, "02_cellFastCNVresults.csv"),
          row.names = FALSE)

epi_seu$fastcnv_score <- NA_real_
epi_seu$fastcnv_threshold <- NA_real_
epi_seu$fastcnv_relative_score <- NA_real_
epi_seu$fastcnv_label <- NA_character_

idx <- match(score_df$Cell, colnames(epi_seu))
ok <- !is.na(idx)
epi_seu$fastcnv_score[idx[ok]] <- score_df$FastCNV_score[ok]
epi_seu$fastcnv_threshold[idx[ok]] <- score_df$Threshold[ok]
epi_seu$fastcnv_relative_score[idx[ok]] <- score_df$Relative_score[ok]
epi_seu$fastcnv_label[idx[ok]] <- score_df$FastCNV_label[ok]

epi_seu$fastcnv_validation <- case_when(
  is.na(epi_seu$fastcnv_label) ~ "No_result",
  epi_seu$tumorselect_label == "Tumor" &
    epi_seu$fastcnv_label == "Tumor" ~ "Concordant_Tumor",
  epi_seu$tumorselect_label == "Normal" &
    epi_seu$fastcnv_label == "Normal" ~ "Concordant_Normal",
  TRUE ~ "Discordant"
)

agreement <- as.data.frame(table(
  PhenoIter = epi_seu$tumorselect_label,
  FastCNV = epi_seu$fastcnv_label,
  useNA = "ifany"
))
write.csv(agreement,
          file.path(out_dir, "03_PhenoIterFastCNVconsistency.csv"),
          row.names = FALSE)

sample_agreement <- epi_seu[[]] %>%
  mutate(Sample = .data[[sample_col]]) %>%
  filter(!is.na(fastcnv_label)) %>%
  group_by(Sample) %>%
  summarise(
    CellN = n(),
    Agreement = mean(tumorselect_label == fastcnv_label),
    PhenoIter_Tumor = mean(tumorselect_label == "Tumor"),
    FastCNV_Tumor = mean(fastcnv_label == "Tumor"),
    .groups = "drop"
  )
write.csv(sample_agreement,
          file.path(out_dir, "03_sampleconsistency.csv"),
          row.names = FALSE)

p_score <- FeaturePlot(
  epi_seu, features = "fastcnv_score",
  reduction = "umap", pt.size = 0.01, order = TRUE
) +
  scale_color_gradientn(colors = c("grey90", "#efb421", "#b43665")) +
  ggtitle("FastCNV score")

p_relative <- FeaturePlot(
  epi_seu, features = "fastcnv_relative_score",
  reduction = "umap", pt.size = 0.01, order = TRUE
) +
  scale_color_gradient2(low = "#6fa6cf", mid = "grey90",
                        high = "#b43665", midpoint = 0) +
  ggtitle("FastCNV score relative to sample threshold")

ggsave(file.path(out_dir, "04_FastCNV_Score_UMAP.pdf"),
       p_score | p_relative, width = 13, height = 6)

p_label <- DimPlot(
  epi_seu, reduction = "umap", group.by = "fastcnv_label",
  cols = c("Tumor" = "#b43665", "Normal" = "#6fa6cf"),
  pt.size = 0.02
) + ggtitle("FastCNV classification")

p_validation <- DimPlot(
  epi_seu, reduction = "umap", group.by = "fastcnv_validation",
  cols = c("Concordant_Tumor" = "#b43665",
           "Concordant_Normal" = "#6fa6cf",
           "Discordant" = "#efb421",
           "No_result" = "grey85"),
  pt.size = 0.02
) + ggtitle("PhenoIter and FastCNV agreement")

ggsave(file.path(out_dir, "05_FastCNV_Label_Validation_UMAP.pdf"),
       p_label | p_validation, width = 13, height = 6)

epi_vln <- epi_seu[[]] %>%
  transmute(Group = tumorselect_label,
            RelativeScore = fastcnv_relative_score) %>%
  filter(!is.na(RelativeScore))
ref_vln <- ref_df %>%
  transmute(Group = ReferenceType,
            RelativeScore = Relative_score)
vln_df <- bind_rows(epi_vln, ref_vln)

group_cols <- c(
  "Tumor" = "#b43665", "Normal" = "#6fa6cf",
  "T cells" = "#187d79", "B cells" = "#8560af",
  "Macrophages" = "#efb421", "DCs" = "#c3347e",
  "Mast cells" = "#e58027", "Endothelial cells" = "#69ae31"
)

p_vln <- ggplot(vln_df, aes(Group, RelativeScore, fill = Group)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = group_cols) +
  labs(title = "FastCNV validation of PhenoIter labels",
       x = NULL,
       y = "CNV score - sample-specific 99% reference threshold") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(out_dir, "06_PhenoIter_FastCNV_Violin.pdf"),
       p_vln, width = 9, height = 6)

epi_seu@misc$FastCNV_validation <- list(
  reference_types = ref_types,
  threshold_method = "Per-sample reference 99th percentile",
  threshold_table = threshold_df,
  sample_agreement = sample_agreement
)

saveRDS(
  epi_seu,
  file.path(out_dir, "GSE245601_Epithelial_TumorSelect_FastCNV.rds")
)

high_conf_tumor <- subset(
  epi_seu,
  subset = fastcnv_validation == "Concordant_Tumor"
)
saveRDS(
  high_conf_tumor,
  file.path(out_dir,
            "GSE245601_HighConfidence_Malignant_Epithelial.rds")
)
