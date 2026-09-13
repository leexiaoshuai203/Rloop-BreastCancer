# RL-Sig65 functional and grade analyses

PROJECT_DIR <- "."

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggprism)
library(ReactomePA)
library(reactome.db)

input_file <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")
output_dir <- file.path(PROJECT_DIR, "results/bulk/discovery/05_gene/RL_Sig65enrichmentanalysis")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)

gene_list <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)[, 1]

gene_ids <- bitr(gene_list,
                 fromType = "SYMBOL",
                 toType = c("ENTREZID", "ENSEMBL"),
                 OrgDb = org.Hs.eg.db)

go_bp <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

go_cc <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

go_mf <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

tryCatch({
  reactome_enrich <- enrichPathway(gene = gene_ids$ENTREZID,
                                   organism = "human",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.2)

  reactome_enrich <- setReadable(reactome_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}, error = function(e) {
  reactome_enrich <<- NULL
})

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  go_bp_df <- go_bp@result %>% mutate(ONTOLOGY = "BP")
} else {
  go_bp_df <- data.frame()
}

if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
  go_cc_df <- go_cc@result %>% mutate(ONTOLOGY = "CC")
} else {
  go_cc_df <- data.frame()
}

if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
  go_mf_df <- go_mf@result %>% mutate(ONTOLOGY = "MF")
} else {
  go_mf_df <- data.frame()
}

go_all <- bind_rows(go_bp_df, go_cc_df, go_mf_df)
if (nrow(go_all) > 0) {
  write.csv(go_all, "GO_enrich.csv", row.names = FALSE)
}

if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  kegg_df <- kegg_enrich@result
  write.csv(kegg_df, "KEGG_enrich.csv", row.names = FALSE)
} else {
  kegg_df <- data.frame()
}

if (!is.null(reactome_enrich) && nrow(reactome_enrich@result) > 0) {
  reactome_df <- reactome_enrich@result %>% mutate(ONTOLOGY = "Reactome")
  write.csv(reactome_df, "Reactome_enrich.csv", row.names = FALSE)
} else {
  reactome_df <- data.frame()
}

go_df   <- read.csv("GO_enrich.csv", header = TRUE, stringsAsFactors = FALSE)
kegg_df <- read.csv("KEGG_enrich.csv", header = TRUE, stringsAsFactors = FALSE)

if (file.exists("Reactome_enrich.csv")) {
  reactome_df <- read.csv("Reactome_enrich.csv", header = TRUE, stringsAsFactors = FALSE)
  has_reactome <- TRUE
} else {
  reactome_df <- data.frame()
  has_reactome <- FALSE
}

if (has_reactome && nrow(reactome_df) > 0) {
  ontology_level <- c("BP", "CC", "MF", "KEGG", "Reactome")
  pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39C12")
} else {
  ontology_level <- c("BP", "CC", "MF", "KEGG")
  pal <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
}

go_selected <- go_df %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup()

kegg_selected <- kegg_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  mutate(ONTOLOGY = "KEGG")

if (has_reactome && nrow(reactome_df) > 0) {
  reactome_selected <- reactome_df %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    mutate(ONTOLOGY = "Reactome")
  pathway_list <- list(go_selected, kegg_selected, reactome_selected)
} else {
  pathway_list <- list(go_selected, kegg_selected)
}

use_pathway <- bind_rows(pathway_list) %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(ontology_level))) %>%
  arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = unique(Description)),
         geneID = map_chr(str_split(geneID, "/"), ~ paste(head(.x, 10), collapse = "/"))) %>%
  tibble::rowid_to_column("index")

width     <- 1
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 5

rect.data <- use_pathway %>%
  count(ONTOLOGY, name = "n") %>%
  mutate(
    xmin = -4 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

p <- ggplot(use_pathway,
            aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_col(aes(y = Description), width = 0.6, alpha = 0.8) +
  geom_text(aes(x = 0.05, label = Description),
            hjust = 0, size = 5) +
  geom_text(aes(x = 0.1, label = geneID, colour = ONTOLOGY),
            hjust = 0, vjust = 2.4, size = 3.5,
            fontface = "italic", show.legend = FALSE) +
  geom_point(aes(x = -width, size = Count), shape = 21) +
  geom_text(aes(x = -width, label = Count), size = 3) +
  scale_size_continuous(name = "Count", range = c(5, 12)) +
  geom_rect(aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax, fill = ONTOLOGY),
            data = rect.data, inherit.aes = FALSE) +
  geom_text(aes(x = (xmin + xmax) / 2,
                y = (ymin + ymax) / 2, label = ONTOLOGY),
            data = rect.data, angle = 90, inherit.aes = FALSE) +
  annotate("segment",
           x = 0, xend = xaxis_max, y = 0, yend = 0,
           linewidth = 1.5) +
  labs(y = NULL, x = "-log10(p.adjust)") +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, xaxis_max, 10),
                     expand = expansion(mult = c(0.05, 0))) +
  guides(fill  = guide_legend(reverse = TRUE),
         colour = guide_legend(reverse = TRUE)) +
  theme_prism() +
  theme(axis.text.y = element_blank(),
        axis.line   = element_blank(),
        axis.ticks.y= element_blank())

ggsave("RL_Sig65_enrichment_barplot.pdf", p, width = 8, height = 12)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(GSVA)
  library(GSEABase)
  library(ggplot2)
  library(ggpubr)
})

select <- dplyr::select
filter <- dplyr::filter

base_dir <- file.path(PROJECT_DIR, "results")
input_gene_file <- file.path(base_dir, "bulk/discovery/05_gene/Final_validated_Rloop_signature_genes.txt")
output_base_dir <- file.path(base_dir, "bulk/discovery/05_gene/RL-Sig65MIK67")

datasets <- list(
  GSE96058 = list(
    expr_file = file.path(base_dir, "bulk/GSE96058/GSE96058_exp.csv"),
    clin_file = file.path(base_dir, "bulk/GSE96058/GSE96058_HRp_HERn_clinical.csv"),
    output_dir = file.path(output_base_dir, "GSE96058"),
    grade_pattern = "G",
    remove_sample_suffix = FALSE,
    phenotypes = c("grade", "PAM50")
  ),
  GSE81538 = list(
    expr_file = file.path(base_dir, "bulk/GSE81538/GSE81538_exp.csv"),
    clin_file = file.path(base_dir, "bulk/GSE81538/GSE81538_HRp_HERn_clinical.csv"),
    output_dir = file.path(output_base_dir, "GSE81538"),
    grade_pattern = "",
    remove_sample_suffix = FALSE,
    phenotypes = c("grade", "PAM50")
  ),
  GSE25066 = list(
    expr_file = file.path(base_dir, "bulk/GSE25066/GSE25066_exp.csv"),
    clin_file = file.path(base_dir, "bulk/GSE25066/GSE25066_HRp_HERn_clinical.csv"),
    output_dir = file.path(output_base_dir, "GSE25066"),
    grade_pattern = "",
    remove_sample_suffix = FALSE,
    phenotypes = c("grade", "PAM50")
  ),
  METABRIC = list(
    expr_file = file.path(base_dir, "bulk/METBRIC/METBRIC_expr_HighLow.csv"),
    clin_file = file.path(base_dir, "bulk/METBRIC/05_METBRIC_HRp_HERn_clinical_2.csv"),
    output_dir = file.path(output_base_dir, "METABRIC"),
    grade_pattern = "",
    remove_sample_suffix = TRUE,
    phenotypes = c("grade")
  )
)

sig_genes <- read.table(input_gene_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
sig_genes <- unique(na.omit(sig_genes))

plot_roc_compare <- function(df, outcome_col, outcome_name, output_dir) {

  dd <- df[, c("Rloop_ssGSEA_Score", "RL_Sig65_score", "MKI67_expr", outcome_col)]
  dd <- dd[complete.cases(dd), ]

  if (nrow(dd) < 10 || length(unique(dd[[outcome_col]])) < 2) {
    return(NULL)
  }

  roc_rloop <- roc(response = dd[[outcome_col]],
                   predictor = dd$Rloop_ssGSEA_Score,
                   direction = "<", quiet = TRUE)

  roc_sig65 <- roc(response = dd[[outcome_col]],
                   predictor = dd$RL_Sig65_score,
                   direction = "<", quiet = TRUE)

  roc_mki67 <- roc(response = dd[[outcome_col]],
                   predictor = dd$MKI67_expr,
                   direction = "<", quiet = TRUE)

  auc_rloop <- as.numeric(auc(roc_rloop))
  ci_rloop  <- ci.auc(roc_rloop)

  auc_sig65 <- as.numeric(auc(roc_sig65))
  ci_sig65  <- ci.auc(roc_sig65)

  auc_mki67 <- as.numeric(auc(roc_mki67))
  ci_mki67  <- ci.auc(roc_mki67)

  p <- ggroc(list("R-loop" = roc_rloop, "RL-Sig65" = roc_sig65, "MKI67" = roc_mki67),
             aes = c("color")) +
    geom_line(size = 1.2, alpha = 0.8) +
    annotate("text", x = 0.3, y = 0.38,
             label = sprintf("R-loop AUC = %.3f\n95%%CI: %.3f-%.3f",
                             auc_rloop, ci_rloop[1], ci_rloop[3]),
             hjust = 0, color = "#218380", size = 4) +
    annotate("text", x = 0.3, y = 0.22,
             label = sprintf("RL-Sig65 AUC = %.3f\n95%%CI: %.3f-%.3f",
                             auc_sig65, ci_sig65[1], ci_sig65[3]),
             hjust = 0, color = "#B5179E", size = 4) +
    annotate("text", x = 0.3, y = 0.06,
             label = sprintf("MKI67 AUC = %.3f\n95%%CI: %.3f-%.3f",
                             auc_mki67, ci_mki67[1], ci_mki67[3]),
             hjust = 0, color = "#F77F00", size = 4) +
    scale_color_manual(values = c("R-loop" = "#218380",
                                  "RL-Sig65" = "#B5179E",
                                  "MKI67" = "#F77F00")) +
    labs(title = outcome_name, x = "Specificity", y = "Sensitivity") +
    theme_bw(base_size = 12) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 13, color = "black"),
      panel.border = element_rect(linewidth = 1.5, color = "black"),
      legend.position = "right"
    )

  ggsave(file.path(output_dir, paste0(outcome_name, "_ROC_compare.pdf")),
         p, width = 6, height = 5)

  return(p)
}

for (dataset_name in names(datasets)) {

  config <- datasets[[dataset_name]]

  if (!dir.exists(config$output_dir)) {
    dir.create(config$output_dir, recursive = TRUE)
  }

  expr_data <- read.csv(config$expr_file, header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE, row.names = 1)
  expr_mat <- as.matrix(expr_data)
  mode(expr_mat) <- "numeric"

  clin <- read.csv(config$clin_file, header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE)

  gene_sets_list <- list(RL_Sig65 = sig_genes)

  ssgsea_param <- ssgseaParam(
    exprData  = expr_mat,
    geneSets  = gene_sets_list,
    alpha     = 0.25,
    normalize = FALSE
  )

  score_mat <- gsva(ssgsea_param, verbose = FALSE)
  score_mat <- t(apply(score_mat, 1, function(x) (x - min(x)) / (max(x) - min(x))))

  score_df <- data.frame(
    sample = colnames(score_mat),
    RL_Sig65_score = as.numeric(score_mat[1, ]),
    stringsAsFactors = FALSE
  )

  if (config$remove_sample_suffix) {
    score_df$sample <- sub("_.*$", "", score_df$sample)
  }

  if (!"MKI67" %in% rownames(expr_mat)) {
    next
  }

  mki67_df <- data.frame(
    sample = colnames(expr_mat),
    MKI67_expr = as.numeric(expr_mat["MKI67", ]),
    stringsAsFactors = FALSE
  )

  if (config$remove_sample_suffix) {
    mki67_df$sample <- sub("_.*$", "", mki67_df$sample)
  }

  merged <- clin %>%
    left_join(score_df, by = "sample") %>%
    left_join(mki67_df, by = "sample") %>%
    filter(!is.na(RL_Sig65_score), !is.na(Rloop_ssGSEA_Score), !is.na(MKI67_expr))

  phenotype_configs <- list()

  if ("grade" %in% config$phenotypes && "grade" %in% colnames(merged)) {
    if (config$grade_pattern == "G") {
      merged <- merged %>%
        mutate(
          Grade_binary = case_when(
            grepl("G3", grade, ignore.case = TRUE) ~ 1,
            grepl("G1|G2", grade, ignore.case = TRUE) ~ 0,
            TRUE ~ NA_real_
          )
        )
    } else {
      merged <- merged %>%
        mutate(
          Grade_binary = case_when(
            grepl("3", grade, ignore.case = TRUE) ~ 1,
            grepl("1|2", grade, ignore.case = TRUE) ~ 0,
            TRUE ~ NA_real_
          )
        )
    }
    phenotype_configs[["Grade_binary"]] <- list(name = "Grade3_vs_Grade12")
  }

  if ("PAM50" %in% config$phenotypes && "PAM50" %in% colnames(merged)) {
    merged <- merged %>%
      mutate(
        PAM50_binary = case_when(
          PAM50 == "LumB" ~ 1,
          PAM50 == "LumA" ~ 0,
          TRUE ~ NA_real_
        )
      )
    phenotype_configs[["PAM50_binary"]] <- list(name = "LumB_vs_LumA")
  }

  for (pheno_col in names(phenotype_configs)) {
    pheno_info <- phenotype_configs[[pheno_col]]
    plot_roc_compare(merged, pheno_col, pheno_info$name, config$output_dir)
  }

}

WORK_DIR <- file.path(PROJECT_DIR, "results/bulk/discovery/05_gene")

SIGNATURE_FILE <- file.path(WORK_DIR, "Final_validated_Rloop_signature_genes.txt")
EXP_FILE       <- file.path(PROJECT_DIR, "results/bulk/GSE96058/GSE96058_exp.csv")
CANCER_SEA_GMT <- file.path(PROJECT_DIR, "results/bulk/GSE96058/CancerSEA.gmt")

OUTPUT_PREFIX1 <- "RL65_CancerSEA_scatter_7x2"
OUTPUT_PREFIX2 <- "RL65_CancerSEA_lollipop"

POINT_SIZE  <- 1.5
POINT_ALPHA <- 0.5
LINE_WIDTH  <- 0.8
SE_ALPHA    <- 0.25

CORR_TEXT_SIZE <- 4
PVAL_TEXT_SIZE  <- 3.5
FDR_TEXT_SIZE   <- 3.3

AXIS_TEXT_SIZE  <- 9
AXIS_TITLE_SIZE <- 11
STRIP_TEXT_SIZE <- 10

FIGURE_WIDTH1  <- 10
FIGURE_HEIGHT1 <- 22

FIGURE_WIDTH2  <- 10
FIGURE_HEIGHT2 <- 5.5

DPI <- 300

PHENOTYPE_COLORS <- c(
  "Angiogenesis"    = "#c58cbd",
  "Apoptosis"       = "#1f9953",
  "Cell cycle"      = "#efb421",
  "Differentiation" = "#e58027",
  "DNA damage"      = "#7daada",
  "DNA repair"      = "#4779bd",
  "EMT"             = "#f38185",
  "Hypoxia"         = "#f9ba4e",
  "Inflammation"    = "#1faa9f",
  "Invasion"        = "#c49874",
  "Metastasis"      = "#ee194b",
  "Proliferation"   = "#b43665",
  "Quiescence"      = "#6fa6cf",
  "Stemness"        = "#187d79"
)

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggforce)
  library(tidyverse)
  library(cowplot)
})

setwd(WORK_DIR)

sig_raw <- fread(SIGNATURE_FILE, data.table = FALSE, header = FALSE)

signature_genes <- unique(na.omit(trimws(as.character(sig_raw[, 1]))))
signature_genes <- signature_genes[signature_genes != ""]

if (length(signature_genes) < 2) {
  stop("signature genes number, check Final_validated_Rloop_signature_genes.txt .")
}

expr_data <- fread(EXP_FILE, data.table = FALSE)

gene_names <- as.character(expr_data[, 1])
expr_data   <- expr_data[, -1, drop = FALSE]

expr_matrix <- as.matrix(expr_data)
mode(expr_matrix) <- "numeric"
rownames(expr_matrix) <- gene_names

expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)) & rownames(expr_matrix) != "", , drop = FALSE]

if (any(duplicated(rownames(expr_matrix)))) {
  expr_matrix <- as.matrix(aggregate(expr_matrix,
                                     by = list(Gene = rownames(expr_matrix)),
                                     FUN = mean, na.rm = TRUE))
  rownames(expr_matrix) <- expr_matrix[, 1]
  expr_matrix <- expr_matrix[, -1, drop = FALSE]
  mode(expr_matrix) <- "numeric"
}

signature_gene_sets <- list(RL_65 = signature_genes)

zscore_param <- zscoreParam(exprData = expr_matrix, geneSets = signature_gene_sets)
sig_score_mat <- gsva(zscore_param, verbose = FALSE)

signature_scores <- as.numeric(sig_score_mat[1, ])
names(signature_scores) <- colnames(expr_matrix)
signature_scores <- signature_scores[!is.na(signature_scores)]

signature_scores <- as.numeric(scale(signature_scores))
names(signature_scores) <- colnames(expr_matrix)

gene_sets_gmt <- getGmt(CANCER_SEA_GMT, geneIdType = SymbolIdentifier())

gene_sets_list <- list()
for (i in seq_along(gene_sets_gmt)) {
  gs <- gene_sets_gmt[[i]]
  gene_sets_list[[setName(gs)]] <- geneIds(gs)
}

zscore_param_cancersea <- zscoreParam(exprData = expr_matrix, geneSets = gene_sets_list)
gsva_result <- gsva(zscore_param_cancersea, verbose = FALSE)

phenotype_scores <- t(scale(t(gsva_result)))

common_samples <- intersect(names(signature_scores), colnames(phenotype_scores))

if (length(common_samples) == 0) {
  names(signature_scores) <- trimws(gsub('^"|"$', "", names(signature_scores)))
  colnames(phenotype_scores) <- trimws(gsub('^"|"$', "", colnames(phenotype_scores)))
  common_samples <- intersect(names(signature_scores), colnames(phenotype_scores))
}

if (length(common_samples) < 10) {
  stop("sample < 10, checksample.")
}

sig_vec   <- signature_scores[common_samples]
pheno_mat <- phenotype_scores[, common_samples, drop = FALSE]

phenotype_names <- rownames(pheno_mat)

cor_results <- data.frame(
  Phenotype   = phenotype_names,
  Correlation = numeric(length(phenotype_names)),
  Pvalue      = numeric(length(phenotype_names)),
  FDR         = numeric(length(phenotype_names)),
  stringsAsFactors = FALSE
)

for (i in seq_along(phenotype_names)) {
  pheno_vec <- as.numeric(pheno_mat[i, ])
  cor_test <- cor.test(sig_vec, pheno_vec, method = "spearman")
  cor_results$Correlation[i] <- as.numeric(cor_test$estimate)
  cor_results$Pvalue[i]      <- cor_test$p.value
}

cor_results$FDR <- p.adjust(cor_results$Pvalue, method = "BH")

write.csv(cor_results,
          file = file.path(WORK_DIR, "RL65_CancerSEA_correlation.csv"),
          row.names = FALSE, quote = FALSE)

sig_range   <- range(sig_vec, na.rm = TRUE)
pheno_range <- range(as.vector(pheno_mat), na.rm = TRUE)

y_limits <- c(floor(sig_range[1] - 0.2), ceiling(sig_range[2] + 0.2))
x_limits <- c(floor(pheno_range[1] - 0.2), ceiling(pheno_range[2] + 0.2))

y_breaks <- seq(y_limits[1], y_limits[2], length.out = 5)
x_breaks <- seq(x_limits[1], x_limits[2], length.out = 5)

make_scatter <- function(pheno_name, pheno_vec, sig_vec,
                         cor_val, p_val, fdr_val, color,
                         is_left_col, is_bottom_row,
                         x_limits, y_limits, x_breaks, y_breaks) {

  df <- data.frame(
    PhenotypeScore = pheno_vec,
    SignatureScore = sig_vec
  )

  format_p <- function(p) {
    if (p < 1e-4) {
      sprintf("p=%.1e", p)
    } else if (p < 0.001) {
      sprintf("p=%.4f", p)
    } else if (p < 0.01) {
      sprintf("p=%.3f", p)
    } else {
      sprintf("p=%.2f", p)
    }
  }

  format_fdr <- function(f) {
    if (f < 1e-4) {
      sprintf("FDR=%.1e", f)
    } else if (f < 0.001) {
      sprintf("FDR=%.4f", f)
    } else if (f < 0.01) {
      sprintf("FDR=%.3f", f)
    } else {
      sprintf("FDR=%.2f", f)
    }
  }

  p_text   <- format_p(p_val)
  fdr_text  <- format_fdr(fdr_val)
  cor_text  <- sprintf("R=%.2f", cor_val)

  p <- ggplot(df, aes(x = PhenotypeScore, y = SignatureScore)) +
    geom_point(color = color, alpha = POINT_ALPHA, size = POINT_SIZE) +
    geom_smooth(method = "lm", se = TRUE,
                color = color, fill = color,
                alpha = SE_ALPHA, linewidth = LINE_WIDTH) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    labs(title = pheno_name, x = NULL, y = NULL) +
    annotate("text",
             x = x_limits[1] + 0.05 * diff(x_limits),
             y = y_limits[2] - 0.06 * diff(y_limits),
             label = cor_text,
             hjust = 0, vjust = 1,
             size = CORR_TEXT_SIZE,
             fontface = "italic") +
    annotate("text",
             x = x_limits[1] + 0.05 * diff(x_limits),
             y = y_limits[2] - 0.18 * diff(y_limits),
             label = p_text,
             hjust = 0, vjust = 1,
             size = PVAL_TEXT_SIZE,
             fontface = "italic") +
    annotate("text",
             x = x_limits[1] + 0.05 * diff(x_limits),
             y = y_limits[2] - 0.30 * diff(y_limits),
             label = fdr_text,
             hjust = 0, vjust = 1,
             size = FDR_TEXT_SIZE,
             fontface = "italic") +
    coord_fixed(ratio = 1.4, clip = "on") +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = STRIP_TEXT_SIZE,
        margin = margin(3, 0, 5, 0)
      ),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line        = element_blank(),
      axis.ticks       = element_line(linewidth = 0.4, color = "black"),
      axis.text        = element_text(size = AXIS_TEXT_SIZE, color = "black"),
      axis.text.y      = if (is_left_col) element_text(size = AXIS_TEXT_SIZE) else element_blank(),
      axis.ticks.y     = if (is_left_col) element_line() else element_blank(),
      axis.text.x      = if (is_bottom_row) element_text(size = AXIS_TEXT_SIZE) else element_blank(),
      axis.ticks.x     = if (is_bottom_row) element_line() else element_blank(),
      plot.margin      = margin(2, 2, 2, 2)
    )

  return(p)
}

plot_list <- list()

for (i in seq_along(phenotype_names)) {
  pname  <- phenotype_names[i]
  pvec   <- as.numeric(pheno_mat[i, ])
  corval <- cor_results$Correlation[i]
  pval   <- cor_results$Pvalue[i]
  fdrval <- cor_results$FDR[i]

  col <- PHENOTYPE_COLORS[pname]
  if (is.na(col)) col <- "#6495ED"

  is_left_col   <- (i - 1) %% 2 == 0
  is_bottom_row <- i > 12

  plot_list[[i]] <- make_scatter(
    pname, pvec, sig_vec,
    corval, pval, fdrval, col,
    is_left_col, is_bottom_row,
    x_limits, y_limits, x_breaks, y_breaks
  )
}

scatter_grid <- plot_grid(
  plotlist = plot_list,
  ncol = 2,
  align = "hv"
)

scatter_final <- ggdraw() +
  draw_plot(scatter_grid, x = 0.08, y = 0.08, width = 0.90, height = 0.88) +
  draw_label("RL-65 signature score (z-score)",
             x = 0.02, y = 0.52, angle = 90,
             size = 12, fontface = "bold") +
  draw_label("CancerSEA phenotype score (z-score)",
             x = 0.53, y = 0.02,
             size = 12, fontface = "bold")

pdf_file1 <- file.path(WORK_DIR, paste0(OUTPUT_PREFIX1, ".pdf"))

ggsave(pdf_file1, scatter_final,
       width = FIGURE_WIDTH1,
       height = FIGURE_HEIGHT1,
       dpi = DPI)

cor_results$Significance <- "NS"
cor_results$Significance[cor_results$FDR < 0.05]  <- "*"
cor_results$Significance[cor_results$FDR < 0.01]  <- "**"
cor_results$Significance[cor_results$FDR < 0.001] <- "***"

lolli_df <- cor_results %>%
  dplyr::arrange(Correlation) %>%
  dplyr::mutate(
    Phenotype = factor(Phenotype, levels = Phenotype)
  )

lolli_colors <- PHENOTYPE_COLORS[as.character(lolli_df$Phenotype)]
names(lolli_colors) <- as.character(lolli_df$Phenotype)

p1 <- ggplot(lolli_df, aes(x = Phenotype, y = Correlation, color = Phenotype)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "grey70") +
  geom_segment(aes(x = Phenotype, xend = Phenotype, y = 0, yend = Correlation),
               linewidth = 1.4) +
  geom_point(size = 4) +
  geom_text(aes(label = Significance),
            vjust = ifelse(lolli_df$Correlation >= 0, -0.8, 1.4),
            size = 5, color = "black") +
  scale_color_manual(values = lolli_colors) +
  scale_y_continuous(
    limits = c(min(-0.1, min(lolli_df$Correlation) - 0.1),
               max(0.1,  max(lolli_df$Correlation) + 0.12)),
    breaks = seq(-1, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(y = "Spearman correlation coefficient", x = NULL) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 16, color = "grey10"),
    axis.text.y  = element_text(size = 13, color = "#808181"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y  = element_line(linewidth = 1.2, color = "#808181"),
    axis.ticks.y = element_line(linewidth = 1.2, color = "#808181"),
    axis.ticks.length.y = unit(0.2, "cm"),
    plot.margin = margin(5, 5, 0, 5)
  )

p2 <- ggplot(lolli_df, aes(x = Phenotype, y = 0.1, fill = Phenotype)) +
  geom_tile(width = 0.95, height = 0.2) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = lolli_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_rect(linewidth = 1.2),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      size = 12,
      color = "#808181",
      family = "sans",
      hjust = 1,
      vjust = 1,
      angle = 45
    ),
    plot.margin = margin(0, 5, 5, 5)
  )

lollipops <- p1 / p2 + plot_layout(heights = c(1, 0.08))

pdf_file2 <- file.path(WORK_DIR, paste0(OUTPUT_PREFIX2, ".pdf"))

ggsave(file = pdf_file2, plot = lollipops,
       width = FIGURE_WIDTH2,
       height = FIGURE_HEIGHT2,
       dpi = DPI)

gc()
set.seed(1234)

required_pkgs <- c(
  "dplyr", "tidyr", "readr", "tibble", "purrr",
  "ggplot2", "pROC", "GSVA", "GSEABase", "patchwork"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(pROC)
  library(GSVA)
  library(GSEABase)
  library(patchwork)
})

select <- dplyr::select
filter <- dplyr::filter

WORK_DIR <- file.path(PROJECT_DIR, "results/bulk/discovery")
PROJECT_DIR <- dirname(WORK_DIR)

SIGNATURE_FILE <- file.path(
  WORK_DIR,
  "05_gene",
  "Final_validated_Rloop_signature_genes.txt"
)

OUT_DIR <- file.path(WORK_DIR, "07_FIG2D_ROC_DeLong")
FIG_DIR <- file.path(OUT_DIR, "Figures")
TAB_DIR <- file.path(OUT_DIR, "Tables")
AUDIT_DIR <- file.path(OUT_DIR, "Audit")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

ROC_DIRECTION <- "<"

NI_MARGIN <- 0.05

datasets <- list(
  GSE96058 = list(
    expr_file = file.path(
      PROJECT_DIR,
      "bulk/GSE96058/GSE96058_exp.csv"
    ),
    clin_file = file.path(
      PROJECT_DIR,
      "bulk/GSE96058/GSE96058_HRp_HERn_clinical.csv"
    ),
    grade_pattern = "G",
    remove_sample_suffix = FALSE
  ),
  GSE81538 = list(
    expr_file = file.path(
      PROJECT_DIR,
      "bulk/GSE81538/GSE81538_exp.csv"
    ),
    clin_file = file.path(
      PROJECT_DIR,
      "bulk/GSE81538/GSE81538_HRp_HERn_clinical.csv"
    ),
    grade_pattern = "",
    remove_sample_suffix = FALSE
  ),
  GSE25066 = list(
    expr_file = file.path(
      PROJECT_DIR,
      "bulk/GSE25066/GSE25066_exp.csv"
    ),
    clin_file = file.path(
      PROJECT_DIR,
      "bulk/GSE25066/GSE25066_HRp_HERn_clinical.csv"
    ),
    grade_pattern = "",
    remove_sample_suffix = FALSE
  ),
  METABRIC = list(
    expr_file = file.path(
      PROJECT_DIR,
      "bulk/METBRIC/METBRIC_expr_HighLow.csv"
    ),
    clin_file = file.path(
      PROJECT_DIR,
      "bulk/METBRIC/05_METBRIC_HRp_HERn_clinical_2.csv"
    ),
    grade_pattern = "",
    remove_sample_suffix = TRUE
  )
)

check_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found:\n", path)
  }
}

clean_sample_id <- function(x, remove_suffix = FALSE) {
  x <- trimws(as.character(x))
  x <- gsub('^"|"$', "", x)
  if (remove_suffix) {
    x <- sub("_.*$", "", x)
  }
  x
}

safe_minmax <- function(x) {
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) {
    return(rep(0.5, length(x)))
  }
  (x - r[1]) / (r[2] - r[1])
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
  sprintf("%.4f", p)
}

extract_grade_binary <- function(x, pattern = "") {
  x <- trimws(as.character(x))
  dplyr::case_when(
    pattern == "G" & grepl("G3", x, ignore.case = TRUE) ~ 1,
    pattern == "G" & grepl("G1|G2", x, ignore.case = TRUE) ~ 0,
    pattern != "G" & grepl("(^|[^0-9])3([^0-9]|$)", x, ignore.case = TRUE) ~ 1,
    pattern != "G" & grepl("(^|[^0-9])(1|2)([^0-9]|$)", x, ignore.case = TRUE) ~ 0,
    TRUE ~ NA_real_
  )
}

read_expression_matrix <- function(path) {
  check_file(path, "Expression file")
  x <- read.csv(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  mat <- as.matrix(x)
  mode(mat) <- "numeric"

  if (anyDuplicated(rownames(mat)) > 0) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE) /
      as.numeric(table(factor(rownames(mat), levels = unique(rownames(mat)))))
  }

  mat
}

collapse_duplicate_samples <- function(df, value_cols) {
  if (!anyDuplicated(df$sample)) return(df)

  df %>%
    group_by(sample) %>%
    summarise(
      across(all_of(value_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
}

calculate_rlsig65 <- function(expr_mat, signature_genes) {
  matched <- intersect(signature_genes, rownames(expr_mat))

  if (length(matched) < 10) {
    stop(
      "Too few RL-Sig65 genes matched the expression matrix: ",
      length(matched), "/", length(signature_genes)
    )
  }

  param <- ssgseaParam(
    exprData = expr_mat,
    geneSets = list(RL_Sig65 = matched),
    alpha = 0.25,
    normalize = FALSE
  )

  score_mat <- gsva(param, verbose = FALSE)
  score <- safe_minmax(as.numeric(score_mat[1, ]))

  list(
    score = data.frame(
      sample = colnames(expr_mat),
      RL_Sig65_score = score,
      stringsAsFactors = FALSE
    ),
    n_matched = length(matched),
    matched_genes = matched
  )
}

make_roc <- function(outcome, predictor) {
  pROC::roc(
    response = outcome,
    predictor = predictor,
    levels = c(0, 1),
    direction = ROC_DIRECTION,
    quiet = TRUE,
    auc = TRUE
  )
}

roc_summary_row <- function(dataset, predictor_name, roc_obj, n_total, n_control, n_case) {
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95, method = "delong"))

  tibble(
    Dataset = dataset,
    Outcome = "Grade 3 vs Grade 1-2",
    Predictor = predictor_name,
    N = n_total,
    N_Grade1_2 = n_control,
    N_Grade3 = n_case,
    AUC = auc_val,
    AUC_CI_low = ci_val[1],
    AUC_CI_high = ci_val[3]
  )
}

delong_compare <- function(dataset, name1, roc1, name2, roc2) {
  test <- pROC::roc.test(
    roc1,
    roc2,
    method = "delong",
    paired = TRUE,
    alternative = "two.sided",
    conf.level = 0.95
  )

  auc1 <- as.numeric(pROC::auc(roc1))
  auc2 <- as.numeric(pROC::auc(roc2))

  ci_diff <- if (!is.null(test$conf.int)) {
    as.numeric(test$conf.int)
  } else {
    c(NA_real_, NA_real_)
  }

  tibble(
    Dataset = dataset,
    Predictor_1 = name1,
    Predictor_2 = name2,
    AUC_1 = auc1,
    AUC_2 = auc2,
    Delta_AUC_1_minus_2 = auc1 - auc2,
    Delta_AUC_CI_low = ci_diff[1],
    Delta_AUC_CI_high = ci_diff[2],
    DeLong_Z = unname(as.numeric(test$statistic)),
    P_value = test$p.value
  )
}

fit_incremental_logistic <- function(dataset, dd) {
  model_df <- dd %>%
    transmute(
      Grade_binary,
      z_RL_Sig65 = as.numeric(scale(RL_Sig65_score)),
      z_MKI67 = as.numeric(scale(MKI67_expr))
    )

  fit_mki67 <- glm(
    Grade_binary ~ z_MKI67,
    data = model_df,
    family = binomial()
  )

  fit_combined <- glm(
    Grade_binary ~ z_MKI67 + z_RL_Sig65,
    data = model_df,
    family = binomial()
  )

  lrt <- anova(fit_mki67, fit_combined, test = "LRT")
  coef_tab <- summary(fit_combined)$coefficients

  beta <- coef_tab["z_RL_Sig65", "Estimate"]
  se <- coef_tab["z_RL_Sig65", "Std. Error"]
  p <- coef_tab["z_RL_Sig65", "Pr(>|z|)"]

  tibble(
    Dataset = dataset,
    N = nrow(model_df),
    RL_Sig65_adjusted_OR_per_SD = exp(beta),
    OR_CI_low = exp(beta - 1.96 * se),
    OR_CI_high = exp(beta + 1.96 * se),
    RL_Sig65_adjusted_P = p,
    LRT_P_MKI67_vs_MKI67_plus_RL_Sig65 = lrt$`Pr(>Chi)`[2],
    AIC_MKI67 = AIC(fit_mki67),
    AIC_MKI67_plus_RL_Sig65 = AIC(fit_combined)
  )
}

plot_roc_panel <- function(dataset, roc_list, delong_df) {
  key_tests <- delong_df %>%
    filter(
      Predictor_1 == "RL-Sig65",
      Predictor_2 %in% c("Parent R-loop-related score", "MKI67 expression")
    )

  p_parent <- key_tests %>%
    filter(Predictor_2 == "Parent R-loop-related score") %>%
    pull(P_value)

  p_mki67 <- key_tests %>%
    filter(Predictor_2 == "MKI67 expression") %>%
    pull(P_value)

  subtitle_text <- paste0(
    "DeLong: RL-Sig65 vs parent P=", format_p(p_parent),
    "; vs MKI67 P=", format_p(p_mki67)
  )

  auc_labels <- map_chr(names(roc_list), function(nm) {
    sprintf("%s, AUC %.3f", nm, as.numeric(pROC::auc(roc_list[[nm]])))
  })
  names(auc_labels) <- names(roc_list)

  pROC::ggroc(roc_list, aes = "color", legacy.axes = FALSE, linewidth = 1.0) +
    scale_color_manual(
      values = c(
        "RL-Sig65" = "#B5179E",
        "Parent R-loop-related score" = "#218380",
        "MKI67 expression" = "#F77F00"
      ),
      labels = auc_labels
    ) +
    geom_abline(
      intercept = 1,
      slope = 1,
      linetype = 2,
      linewidth = 0.5,
      color = "grey55"
    ) +
    labs(
      title = dataset,
      subtitle = subtitle_text,
      x = "Specificity",
      y = "Sensitivity",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.8, color = "black")
    ) +
    coord_equal()
}

check_file(SIGNATURE_FILE, "RL-Sig65 gene file")

signature_genes <- read.table(
  SIGNATURE_FILE,
  header = FALSE,
  stringsAsFactors = FALSE
)[, 1]

signature_genes <- unique(trimws(as.character(signature_genes)))
signature_genes <- signature_genes[
  !is.na(signature_genes) & signature_genes != ""
]

all_roc_summary <- list()
all_delong <- list()
all_incremental <- list()
all_audit <- list()
all_plots <- list()

for (dataset_name in names(datasets)) {

  cfg <- datasets[[dataset_name]]
  check_file(cfg$expr_file, paste0(dataset_name, " expression file"))
  check_file(cfg$clin_file, paste0(dataset_name, " clinical file"))

  expr_mat <- read_expression_matrix(cfg$expr_file)

  colnames(expr_mat) <- clean_sample_id(
    colnames(expr_mat),
    remove_suffix = cfg$remove_sample_suffix
  )

  clin <- read.csv(
    cfg$clin_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_clin_cols <- c("sample", "grade", "Rloop_ssGSEA_Score")
  missing_clin_cols <- setdiff(required_clin_cols, colnames(clin))
  if (length(missing_clin_cols) > 0) {
    stop(
      dataset_name, " clinical file is missing columns: ",
      paste(missing_clin_cols, collapse = ", ")
    )
  }

  clin$sample <- clean_sample_id(
    clin$sample,
    remove_suffix = cfg$remove_sample_suffix
  )

  if (anyDuplicated(clin$sample) > 0) {
    warning(dataset_name, ": duplicated clinical sample IDs detected; retaining first row.")
    clin <- clin[!duplicated(clin$sample), , drop = FALSE]
  }

  sig_res <- calculate_rlsig65(expr_mat, signature_genes)
  score_df <- sig_res$score
  score_df$sample <- clean_sample_id(
    score_df$sample,
    remove_suffix = cfg$remove_sample_suffix
  )
  score_df <- collapse_duplicate_samples(score_df, "RL_Sig65_score")

  if (!"MKI67" %in% rownames(expr_mat)) {
    warning(dataset_name, ": MKI67 was not found; dataset skipped.")
    next
  }

  mki67_df <- tibble(
    sample = clean_sample_id(
      colnames(expr_mat),
      remove_suffix = cfg$remove_sample_suffix
    ),
    MKI67_expr = as.numeric(expr_mat["MKI67", ])
  )
  mki67_df <- collapse_duplicate_samples(mki67_df, "MKI67_expr")

  merged <- clin %>%
    select(sample, grade, Rloop_ssGSEA_Score, everything()) %>%
    left_join(score_df, by = "sample") %>%
    left_join(mki67_df, by = "sample") %>%
    mutate(
      Grade_binary = extract_grade_binary(grade, cfg$grade_pattern)
    )

  dd <- merged %>%
    select(
      sample,
      grade,
      Grade_binary,
      Rloop_ssGSEA_Score,
      RL_Sig65_score,
      MKI67_expr
    ) %>%
    filter(complete.cases(.)) %>%
    distinct(sample, .keep_all = TRUE)

  n_total <- nrow(dd)
  n_control <- sum(dd$Grade_binary == 0)
  n_case <- sum(dd$Grade_binary == 1)

  if (n_total < 20 || n_control < 5 || n_case < 5) {
    warning(dataset_name, ": insufficient Grade 1-2 or Grade 3 samples; skipped.")
    next
  }

  write.csv(
    dd,
    file.path(AUDIT_DIR, paste0(dataset_name, "_paired_ROC_analysis_population.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  roc_rl65 <- make_roc(dd$Grade_binary, dd$RL_Sig65_score)
  roc_parent <- make_roc(dd$Grade_binary, dd$Rloop_ssGSEA_Score)
  roc_mki67 <- make_roc(dd$Grade_binary, dd$MKI67_expr)

  roc_list <- list(
    "RL-Sig65" = roc_rl65,
    "Parent R-loop-related score" = roc_parent,
    "MKI67 expression" = roc_mki67
  )

  cohort_roc_summary <- bind_rows(
    roc_summary_row(
      dataset_name, "RL-Sig65", roc_rl65,
      n_total, n_control, n_case
    ),
    roc_summary_row(
      dataset_name, "Parent R-loop-related score", roc_parent,
      n_total, n_control, n_case
    ),
    roc_summary_row(
      dataset_name, "MKI67 expression", roc_mki67,
      n_total, n_control, n_case
    )
  )

  cohort_delong <- bind_rows(
    delong_compare(
      dataset_name,
      "RL-Sig65", roc_rl65,
      "Parent R-loop-related score", roc_parent
    ),
    delong_compare(
      dataset_name,
      "RL-Sig65", roc_rl65,
      "MKI67 expression", roc_mki67
    ),
    delong_compare(
      dataset_name,
      "Parent R-loop-related score", roc_parent,
      "MKI67 expression", roc_mki67
    )
  )

  cohort_incremental <- fit_incremental_logistic(dataset_name, dd)

  cohort_audit <- tibble(
    Dataset = dataset_name,
    Expression_genes = nrow(expr_mat),
    Expression_samples = ncol(expr_mat),
    Clinical_rows = nrow(clin),
    RL_Sig65_genes_total = length(signature_genes),
    RL_Sig65_genes_matched = sig_res$n_matched,
    Paired_ROC_N = n_total,
    N_Grade1_2 = n_control,
    N_Grade3 = n_case
  )

  p <- plot_roc_panel(dataset_name, roc_list, cohort_delong)
  ggsave(
    file.path(FIG_DIR, paste0(dataset_name, "_Grade3_ROC_DeLong.pdf")),
    p,
    width = 6.4,
    height = 5.5
  )
  ggsave(
    file.path(FIG_DIR, paste0(dataset_name, "_Grade3_ROC_DeLong.png")),
    p,
    width = 6.4,
    height = 5.5,
    dpi = 600
  )

  all_roc_summary[[dataset_name]] <- cohort_roc_summary
  all_delong[[dataset_name]] <- cohort_delong
  all_incremental[[dataset_name]] <- cohort_incremental
  all_audit[[dataset_name]] <- cohort_audit
  all_plots[[dataset_name]] <- p
}

roc_summary_df <- bind_rows(all_roc_summary)
delong_df <- bind_rows(all_delong)
incremental_df <- bind_rows(all_incremental)
audit_df <- bind_rows(all_audit)

if (nrow(delong_df) == 0) {
  stop("No cohort completed the ROC analysis.")
}

delong_df <- delong_df %>%
  group_by(Dataset) %>%
  mutate(P_BH_within_cohort = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    P_BH_global = p.adjust(P_value, method = "BH"),
    NI_margin_parent = if_else(
      Predictor_1 == "RL-Sig65" &
        Predictor_2 == "Parent R-loop-related score",
      NI_MARGIN,
      NA_real_
    ),
    Exploratory_noninferiority_vs_parent = case_when(
      Predictor_1 == "RL-Sig65" &
        Predictor_2 == "Parent R-loop-related score" &
        !is.na(Delta_AUC_CI_low) &
        Delta_AUC_CI_low > -NI_MARGIN ~ "Yes",
      Predictor_1 == "RL-Sig65" &
        Predictor_2 == "Parent R-loop-related score" &
        !is.na(Delta_AUC_CI_low) ~ "No",
      TRUE ~ NA_character_
    )
  )

write.csv(
  roc_summary_df,
  file.path(TAB_DIR, "Supplementary_Table_ROC_AUC_95CI.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  delong_df,
  file.path(TAB_DIR, "Supplementary_Table_Paired_DeLong_comparisons.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  incremental_df,
  file.path(TAB_DIR, "Supplementary_Table_Incremental_value_over_MKI67.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  audit_df,
  file.path(TAB_DIR, "Supplementary_Table_ROC_analysis_audit.csv"),
  row.names = FALSE,
  quote = FALSE
)

if (length(all_plots) > 0) {
  plot_order <- intersect(
    c("GSE96058", "GSE81538", "METABRIC", "GSE25066"),
    names(all_plots)
  )

  combined_plot <- wrap_plots(
    all_plots[plot_order],
    ncol = 2,
    guides = "collect"
  ) &
    theme(legend.position = "bottom")

  ggsave(
    file.path(FIG_DIR, "Figure_2D_Grade3_ROC_four_cohorts.pdf"),
    combined_plot,
    width = 12.5,
    height = 10.5
  )
  ggsave(
    file.path(FIG_DIR, "Figure_2D_Grade3_ROC_four_cohorts.png"),
    combined_plot,
    width = 12.5,
    height = 10.5,
    dpi = 600
  )
}

forest_df <- delong_df %>%
  filter(
    Predictor_1 == "RL-Sig65",
    Predictor_2 %in% c("Parent R-loop-related score", "MKI67 expression")
  ) %>%
  mutate(
    Comparison = recode(
      Predictor_2,
      "Parent R-loop-related score" = "RL-Sig65 vs parent score",
      "MKI67 expression" = "RL-Sig65 vs MKI67"
    ),
    Dataset = factor(
      Dataset,
      levels = rev(c("GSE96058", "GSE81538", "METABRIC", "GSE25066"))
    )
  )

p_forest <- ggplot(
  forest_df,
  aes(
    x = Delta_AUC_1_minus_2,
    y = Dataset,
    xmin = Delta_AUC_CI_low,
    xmax = Delta_AUC_CI_high
  )
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_errorbarh(height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ Comparison, ncol = 1) +
  labs(
    x = expression(Delta * "AUC (RL-Sig65 minus comparator)"),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(FIG_DIR, "Supplementary_Figure_AUC_difference_forest.pdf"),
  p_forest,
  width = 7.2,
  height = 6.8
)
ggsave(
  file.path(FIG_DIR, "Supplementary_Figure_AUC_difference_forest.png"),
  p_forest,
  width = 7.2,
  height = 6.8,
  dpi = 600
)

print(
  delong_df %>%
    filter(Predictor_1 == "RL-Sig65") %>%
    select(
      Dataset,
      Predictor_2,
      Delta_AUC_1_minus_2,
      Delta_AUC_CI_low,
      Delta_AUC_CI_high,
      P_value,
      P_BH_global,
      Exploratory_noninferiority_vs_parent
    )
)

pkgs <- c("ggplot2", "dplyr", "patchwork", "scales", "Cairo")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

infile  <- file.path(PROJECT_DIR, "results/bulk/discovery/07_FIG2D_ROC_DeLong/Tables/Supplementary_Table_Incremental_value_over_MKI67.csv")
outdir  <- dirname(infile)
outfile <- file.path(outdir, "Forest_Incremental_RL_Sig65_over_MKI67.pdf")

df <- read.csv(infile, check.names = FALSE)

df <- df %>%
  mutate(
    Dataset   = factor(Dataset, levels = rev(Dataset)),
    ypos      = as.numeric(Dataset),
    OR_label  = sprintf("%.2f (%.2f\u2013%.2f)",
                        RL_Sig65_adjusted_OR_per_SD,
                        OR_CI_low, OR_CI_high),
    P_label   = ifelse(RL_Sig65_adjusted_P < 0.0001,
                       formatC(RL_Sig65_adjusted_P, format = "e", digits = 2),
                       sprintf("%.4f", RL_Sig65_adjusted_P)),
    LRT_label = ifelse(LRT_P_MKI67_vs_MKI67_plus_RL_Sig65 < 0.0001,
                       formatC(LRT_P_MKI67_vs_MKI67_plus_RL_Sig65, format = "e", digits = 2),
                       sprintf("%.4f", LRT_P_MKI67_vs_MKI67_plus_RL_Sig65)),
    dAIC      = AIC_MKI67 - AIC_MKI67_plus_RL_Sig65,
    dAIC_label = sprintf("+%.1f", dAIC),
    N_label   = as.character(N)
  )

nrow_data <- nrow(df)
yrange    <- c(0.3, nrow_data + 1.2)

col_point  <- "#1A6DAF"
col_ci     <- "#2980B9"
col_ref    <- "#BDC3C7"
col_head   <- "#2C3E50"
col_or     <- "#1A5276"
col_pval   <- "#922B21"
col_aic    <- "#1A7A4A"
col_stripe <- c("#F8FAFB", "#FFFFFF")

p_left <- ggplot(df, aes(y = ypos)) +

  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = seq(0.5, nrow_data - 0.5, 1),
           ymax = seq(1.5, nrow_data + 0.5, 1),
           fill = col_stripe[1], alpha = 0.7) +

  annotate("text", x = 0.05, y = nrow_data + 0.85,
           label = "Dataset", hjust = 0, vjust = 0,
           size = 4.6, fontface = "bold", color = col_head) +
  annotate("text", x = 0.60, y = nrow_data + 0.85,
           label = "N", hjust = 0.5, vjust = 0,
           size = 4.6, fontface = "bold", color = col_head) +

  geom_text(aes(x = 0.05, label = as.character(Dataset)),
            hjust = 0, size = 4.3, fontface = "bold", color = "#1f2b3e") +

  geom_text(aes(x = 0.60, label = N_label),
            hjust = 0.5, size = 4.1, color = "#444444") +

  scale_x_continuous(limits = c(0, 0.85), expand = c(0, 0)) +
  scale_y_continuous(limits = yrange, expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 8, r = 0, b = 8, l = 10))

p_forest <- ggplot(df, aes(y = ypos, x = RL_Sig65_adjusted_OR_per_SD)) +

  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = seq(0.5, nrow_data - 0.5, 1),
           ymax = seq(1.5, nrow_data + 0.5, 1),
           fill = col_stripe[1], alpha = 0.7) +

  annotate("rect",
           xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#EBF5FB", alpha = 0.4) +

  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.8, color = "#95A5A6") +

  geom_errorbarh(aes(xmin = OR_CI_low, xmax = OR_CI_high),
                 height = 0.22, linewidth = 1.3, color = col_ci) +

  geom_point(aes(size = sqrt(N)),
             shape = 21, stroke = 1.3,
             fill = col_point, color = "#0D3D6B") +

  scale_size_continuous(range = c(3.5, 6.5), guide = "none") +

  scale_x_log10(
    limits = c(1.1, 15),
    breaks = c(1, 2, 3, 5, 10),
    labels = c("1", "2", "3", "5", "10"),
    expand = c(0.02, 0.02)
  ) +
  scale_y_continuous(limits = yrange, expand = c(0, 0)) +

  annotate("text", x = sqrt(1.1 * 15), y = nrow_data + 0.85,
           label = "Odds Ratio (95% CI)", hjust = 0.5, vjust = 0,
           size = 4.6, fontface = "bold", color = col_head) +

  labs(
    x = "Adjusted OR per SD increase (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.5),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(size = 11.5, color = "#333333"),
    axis.title.x = element_text(size = 12, face = "bold",
                                margin = margin(t = 8)),
    plot.margin  = margin(t = 8, r = 6, b = 8, l = 6)
  )

x_or   <- 0.02
x_p    <- 0.42
x_lrt  <- 0.62
x_aic  <- 0.84

p_right <- ggplot(df, aes(y = ypos)) +

  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = seq(0.5, nrow_data - 0.5, 1),
           ymax = seq(1.5, nrow_data + 0.5, 1),
           fill = col_stripe[1], alpha = 0.7) +

  annotate("text", x = x_or,  y = nrow_data + 0.85,
           label = "OR (95% CI)", hjust = 0, vjust = 0,
           size = 4.2, fontface = "bold", color = col_head) +
  annotate("text", x = x_p,   y = nrow_data + 0.85,
           label = "P value", hjust = 0, vjust = 0,
           size = 4.2, fontface = "bold", color = col_head) +
  annotate("text", x = x_lrt, y = nrow_data + 0.85,
           label = "LRT-P", hjust = 0, vjust = 0,
           size = 4.2, fontface = "bold", color = col_head) +
  annotate("text", x = x_aic, y = nrow_data + 0.85,
           label = "\u0394AIC", hjust = 0.5, vjust = 0,
           size = 4.2, fontface = "bold", color = col_head) +

  geom_text(aes(x = x_or, label = OR_label),
            hjust = 0, size = 3.9, color = col_or, fontface = "bold") +

  geom_text(aes(x = x_p, label = P_label),
            hjust = 0, size = 3.8, color = col_pval) +

  geom_text(aes(x = x_lrt, label = LRT_label),
            hjust = 0, size = 3.8, color = col_pval) +

  geom_tile(aes(x = x_aic, width = 0.14, height = 0.62,
                fill = dAIC), color = NA) +
  geom_text(aes(x = x_aic, label = dAIC_label),
            hjust = 0.5, size = 3.9, fontface = "bold", color = "white") +

  scale_fill_gradient(low = "#AED6F1", high = "#1A5276",
                      name = "\u0394AIC", guide = "none") +

  scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  scale_y_continuous(limits = yrange, expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(t = 8, r = 10, b = 8, l = 4))

combined <- p_left + p_forest + p_right +
  plot_layout(widths = c(1.8, 2.8, 3.2)) +
  plot_annotation(
    title    = "Incremental Prognostic Value of RL_Sig65 over MKI67",
    subtitle = "Adjusted odds ratios per SD increase across independent breast cancer cohorts",
    caption  = paste0(
      "Point size reflects cohort sample size.  ",
      "\u0394AIC = AIC(MKI67\u2009model) \u2212 AIC(MKI67\u2009+\u2009RL_Sig65\u2009model); ",
      "larger \u0394AIC indicates greater incremental model improvement.\n",
      "LRT-P: likelihood-ratio test comparing nested models."
    ),
    theme = theme(
      plot.title    = element_text(size = 18, face = "bold",
                                   hjust = 0.5, color = "#1C2833",
                                   margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12.5, hjust = 0.5,
                                   color = "#566573",
                                   margin = margin(b = 12)),
      plot.caption  = element_text(size = 9.5, color = "#717D7E",
                                   hjust = 0.5,
                                   margin = margin(t = 10)),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 20, r = 20, b = 15, l = 20)
    )
  )

ggsave(
  filename = outfile,
  plot     = combined,
  width    = 17,
  height   = 6,
  device   = cairo_pdf,
  dpi      = 300
)

gc()

required_pkgs <- c("openxlsx", "dplyr", "tidyr", "data.table")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

select <- dplyr::select
filter <- dplyr::filter

WORK_DIR <- file.path(PROJECT_DIR, "results/bulk/discovery")
PROJECT_DIR <- dirname(WORK_DIR)

OUT_DIR <- file.path(WORK_DIR, "07_FIG2D_ROC_DeLong")
TABLE_DIR <- file.path(OUT_DIR, "Tables")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_XLSX <- file.path(
  OUT_DIR,
  "Supplementary_Tables_Figure2_RL-Sig65.xlsx"
)

SIGNATURE_FILE <- file.path(
  WORK_DIR,
  "05_gene",
  "Final_validated_Rloop_signature_genes.txt"
)

DEG_FILE <- file.path(
  WORK_DIR,
  "04_C1C2_DEG_Analysis",
  "GSE96058_C2_vs_C1_limma_results.csv"
)

VALIDATION_FILE <- file.path(
  WORK_DIR,
  "05_gene",
  "Candidate_genes_validation_results.csv"
)

RLOOP_GMT_FILE <- file.path(
  PROJECT_DIR,
  "01.R-loopsample",
  "GSE96058",
  "01_Rloop_regulators.gmt"
)

GSE25066_EXPR_FILE <- file.path(
  PROJECT_DIR,
  "01.R-loopsample",
  "GSE25066",
  "GSE25066_exp.csv"
)

ROC_FILE <- file.path(
  TABLE_DIR,
  "Supplementary_Table_ROC_AUC_95CI.csv"
)

DELONG_FILE <- file.path(
  TABLE_DIR,
  "Supplementary_Table_Paired_DeLong_comparisons.csv"
)

AUDIT_FILE <- file.path(
  TABLE_DIR,
  "Supplementary_Table_ROC_analysis_audit.csv"
)

INCREMENTAL_FILE <- file.path(
  TABLE_DIR,
  "Supplementary_Table_Incremental_value_over_MKI67.csv"
)

assert_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "\nMissing ", label, ":\n", path,
      "\n\nPlease confirm that the upstream Figure 2 analysis has been completed."
    )
  }
}

required_paths <- c(
  SIGNATURE_FILE,
  DEG_FILE,
  VALIDATION_FILE,
  RLOOP_GMT_FILE,
  GSE25066_EXPR_FILE,
  ROC_FILE,
  DELONG_FILE,
  AUDIT_FILE,
  INCREMENTAL_FILE
)

required_labels <- c(
  "RL-Sig65 gene list",
  "GSE96058 differential-expression result",
  "candidate-gene validation result",
  "R-loopBase GMT file",
  "GSE25066 expression matrix",
  "ROC AUC result",
  "paired DeLong result",
  "ROC audit result",
  "incremental-model result"
)

invisible(Map(assert_file, required_paths, required_labels))

require_columns <- function(df, required, file_label) {
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "\n", file_label, " is missing the following columns:\n",
      paste(missing_cols, collapse = ", ")
    )
  }
}

clean_gene_vector <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & x != ""])
}

read_gmt_genes <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]

  genes <- unlist(
    lapply(lines, function(one_line) {
      parts <- strsplit(one_line, "\t", fixed = TRUE)[[1]]
      if (length(parts) <= 2) {
        return(character(0))
      }
      parts[-c(1, 2)]
    }),
    use.names = FALSE
  )

  clean_gene_vector(genes)
}

cohort_order <- c("GSE96058", "GSE81538", "METABRIC", "GSE25066")

signature_genes <- read.table(
  SIGNATURE_FILE,
  header = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)[, 1]

signature_genes <- clean_gene_vector(signature_genes)

if (length(signature_genes) != 65) {
  warning(
    "The signature file contains ", length(signature_genes),
    " unique genes rather than 65."
  )
}

deg_df <- read.csv(
  DEG_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

require_columns(
  deg_df,
  c("GeneSymbol", "logFC", "adj.P.Val"),
  basename(DEG_FILE)
)

deg_sub <- deg_df %>%
  transmute(
    Gene = as.character(GeneSymbol),
    GSE96058_log2FC = as.numeric(logFC),
    GSE96058_FDR = as.numeric(adj.P.Val)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

validation_df <- read.csv(
  VALIDATION_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

require_columns(
  validation_df,
  c(
    "Gene",
    "MeanDiff_GSE81538",
    "Correlation_GSE81538",
    "MeanDiff_METABRIC",
    "Correlation_METABRIC"
  ),
  basename(VALIDATION_FILE)
)

validation_sub <- validation_df %>%
  transmute(
    Gene = as.character(Gene),
    GSE81538_mean_difference = as.numeric(MeanDiff_GSE81538),
    GSE81538_Spearman_rho = as.numeric(Correlation_GSE81538),
    METABRIC_mean_difference = as.numeric(MeanDiff_METABRIC),
    METABRIC_Spearman_rho = as.numeric(Correlation_METABRIC)
  ) %>%
  distinct(Gene, .keep_all = TRUE)

rloopbase_genes <- read_gmt_genes(RLOOP_GMT_FILE)

gse25066_gene_column <- data.table::fread(
  GSE25066_EXPR_FILE,
  select = 1,
  data.table = FALSE,
  check.names = FALSE,
  showProgress = FALSE
)[[1]]

gse25066_genes <- clean_gene_vector(gse25066_gene_column)

sheet_genes <- data.frame(
  Gene = signature_genes,
  stringsAsFactors = FALSE
) %>%
  left_join(deg_sub, by = "Gene") %>%
  left_join(validation_sub, by = "Gene") %>%
  mutate(
    RloopBase_overlap = ifelse(Gene %in% rloopbase_genes, "Yes", "No"),
    Available_in_GSE25066 = ifelse(Gene %in% gse25066_genes, "Yes", "No")
  ) %>%
  select(
    Gene,
    GSE96058_log2FC,
    GSE96058_FDR,
    GSE81538_mean_difference,
    GSE81538_Spearman_rho,
    METABRIC_mean_difference,
    METABRIC_Spearman_rho,
    RloopBase_overlap,
    Available_in_GSE25066
  )

roc_df <- read.csv(
  ROC_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

delong_df <- read.csv(
  DELONG_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

audit_df <- read.csv(
  AUDIT_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

require_columns(
  roc_df,
  c(
    "Dataset",
    "Predictor",
    "AUC",
    "AUC_CI_low",
    "AUC_CI_high"
  ),
  basename(ROC_FILE)
)

require_columns(
  delong_df,
  c(
    "Dataset",
    "Predictor_1",
    "Predictor_2",
    "Delta_AUC_1_minus_2",
    "Delta_AUC_CI_low",
    "Delta_AUC_CI_high",
    "P_value"
  ),
  basename(DELONG_FILE)
)

require_columns(
  audit_df,
  c(
    "Dataset",
    "RL_Sig65_genes_total",
    "RL_Sig65_genes_matched",
    "Paired_ROC_N",
    "N_Grade1_2",
    "N_Grade3"
  ),
  basename(AUDIT_FILE)
)

roc_wide <- roc_df %>%
  mutate(
    Predictor_key = dplyr::recode(
      Predictor,
      "RL-Sig65" = "RL_Sig65",
      "Parent R-loop-related score" = "Original_Rloop_related_score",
      "MKI67 expression" = "MKI67_expression"
    )
  ) %>%
  filter(
    Predictor_key %in% c(
      "RL_Sig65",
      "Original_Rloop_related_score",
      "MKI67_expression"
    )
  ) %>%
  select(
    Dataset,
    Predictor_key,
    AUC,
    AUC_CI_low,
    AUC_CI_high
  ) %>%
  pivot_wider(
    names_from = Predictor_key,
    values_from = c(AUC, AUC_CI_low, AUC_CI_high),
    names_glue = "{Predictor_key}_{.value}"
  )

delong_core <- delong_df %>%
  filter(
    Predictor_1 == "RL-Sig65",
    Predictor_2 %in% c(
      "Parent R-loop-related score",
      "MKI67 expression"
    )
  ) %>%
  mutate(
    Comparator_key = dplyr::recode(
      Predictor_2,
      "Parent R-loop-related score" = "Original_Rloop_related_score",
      "MKI67 expression" = "MKI67_expression"
    ),
    BH_adjusted_P = p.adjust(P_value, method = "BH")
  )

if (nrow(delong_core) != 8) {
  warning(
    "Expected 8 core DeLong comparisons, but found ",
    nrow(delong_core), "."
  )
}

delong_wide <- delong_core %>%
  select(
    Dataset,
    Comparator_key,
    Delta_AUC_1_minus_2,
    Delta_AUC_CI_low,
    Delta_AUC_CI_high,
    P_value,
    BH_adjusted_P
  ) %>%
  pivot_wider(
    names_from = Comparator_key,
    values_from = c(
      Delta_AUC_1_minus_2,
      Delta_AUC_CI_low,
      Delta_AUC_CI_high,
      P_value,
      BH_adjusted_P
    ),
    names_glue = "RL_Sig65_vs_{Comparator_key}_{.value}"
  )

audit_sub <- audit_df %>%
  transmute(
    Dataset = as.character(Dataset),
    N = as.integer(Paired_ROC_N),
    N_Grade1_2 = as.integer(N_Grade1_2),
    N_Grade3 = as.integer(N_Grade3),
    RL_Sig65_genes_available = as.integer(RL_Sig65_genes_matched),
    RL_Sig65_genes_total = as.integer(RL_Sig65_genes_total)
  )

sheet_roc <- audit_sub %>%
  left_join(roc_wide, by = "Dataset") %>%
  left_join(delong_wide, by = "Dataset") %>%
  mutate(
    Dataset = factor(Dataset, levels = cohort_order)
  ) %>%
  arrange(Dataset) %>%
  mutate(
    Dataset = as.character(Dataset)
  ) %>%
  select(
    Dataset,
    N,
    N_Grade1_2,
    N_Grade3,
    RL_Sig65_genes_available,
    RL_Sig65_genes_total,

    RL_Sig65_AUC,
    RL_Sig65_AUC_CI_low,
    RL_Sig65_AUC_CI_high,

    Original_Rloop_related_score_AUC,
    Original_Rloop_related_score_AUC_CI_low,
    Original_Rloop_related_score_AUC_CI_high,

    MKI67_expression_AUC,
    MKI67_expression_AUC_CI_low,
    MKI67_expression_AUC_CI_high,

    RL_Sig65_vs_Original_Rloop_related_score_Delta_AUC_1_minus_2,
    RL_Sig65_vs_Original_Rloop_related_score_Delta_AUC_CI_low,
    RL_Sig65_vs_Original_Rloop_related_score_Delta_AUC_CI_high,
    RL_Sig65_vs_Original_Rloop_related_score_P_value,
    RL_Sig65_vs_Original_Rloop_related_score_BH_adjusted_P,

    RL_Sig65_vs_MKI67_expression_Delta_AUC_1_minus_2,
    RL_Sig65_vs_MKI67_expression_Delta_AUC_CI_low,
    RL_Sig65_vs_MKI67_expression_Delta_AUC_CI_high,
    RL_Sig65_vs_MKI67_expression_P_value,
    RL_Sig65_vs_MKI67_expression_BH_adjusted_P
  )

incremental_df <- read.csv(
  INCREMENTAL_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

require_columns(
  incremental_df,
  c(
    "Dataset",
    "N",
    "RL_Sig65_adjusted_OR_per_SD",
    "OR_CI_low",
    "OR_CI_high",
    "RL_Sig65_adjusted_P",
    "LRT_P_MKI67_vs_MKI67_plus_RL_Sig65",
    "AIC_MKI67",
    "AIC_MKI67_plus_RL_Sig65"
  ),
  basename(INCREMENTAL_FILE)
)

sheet_incremental <- incremental_df %>%
  transmute(
    Dataset = factor(Dataset, levels = cohort_order),
    N = as.integer(N),
    RL_Sig65_adjusted_OR_per_SD = as.numeric(RL_Sig65_adjusted_OR_per_SD),
    OR_CI_low = as.numeric(OR_CI_low),
    OR_CI_high = as.numeric(OR_CI_high),
    RL_Sig65_adjusted_P = as.numeric(RL_Sig65_adjusted_P),
    LRT_P_MKI67_vs_MKI67_plus_RL_Sig65 = as.numeric(
      LRT_P_MKI67_vs_MKI67_plus_RL_Sig65
    ),
    AIC_MKI67 = as.numeric(AIC_MKI67),
    AIC_MKI67_plus_RL_Sig65 = as.numeric(AIC_MKI67_plus_RL_Sig65),
    Delta_AIC = as.numeric(AIC_MKI67) -
      as.numeric(AIC_MKI67_plus_RL_Sig65)
  ) %>%
  arrange(Dataset) %>%
  mutate(
    Dataset = as.character(Dataset)
  )

wb <- createWorkbook()

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "Bottom",
  borderColour = "#808080"
)

body_style <- createStyle(
  valign = "center"
)

decimal_style <- createStyle(
  numFmt = "0.000"
)

pvalue_style <- createStyle(
  numFmt = "0.00E+00"
)

integer_style <- createStyle(
  numFmt = "0"
)

write_plain_sheet <- function(
    wb,
    sheet_name,
    data,
    widths
) {
  addWorksheet(wb, sheet_name, gridLines = FALSE)

  writeData(
    wb,
    sheet = sheet_name,
    x = data,
    startRow = 1,
    startCol = 1,
    colNames = TRUE,
    rowNames = FALSE,
    keepNA = FALSE,
    headerStyle = header_style,
    borders = "rows"
  )

  addFilter(
    wb,
    sheet = sheet_name,
    row = 1,
    cols = seq_len(ncol(data))
  )

  freezePane(
    wb,
    sheet = sheet_name,
    firstActiveRow = 2,
    firstActiveCol = 2
  )

  setRowHeights(
    wb,
    sheet = sheet_name,
    rows = 1,
    heights = 42
  )

  setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_len(ncol(data)),
    widths = widths
  )

  if (nrow(data) > 0) {
    addStyle(
      wb,
      sheet = sheet_name,
      style = body_style,
      rows = 2:(nrow(data) + 1),
      cols = seq_len(ncol(data)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }
}

write_plain_sheet(
  wb = wb,
  sheet_name = "RL-Sig65_genes",
  data = sheet_genes,
  widths = c(15, 18, 15, 22, 20, 22, 20, 18, 21)
)

if (nrow(sheet_genes) > 0) {
  addStyle(
    wb,
    "RL-Sig65_genes",
    decimal_style,
    rows = 2:(nrow(sheet_genes) + 1),
    cols = c(2, 4, 5, 6, 7),
    gridExpand = TRUE,
    stack = TRUE
  )
  addStyle(
    wb,
    "RL-Sig65_genes",
    pvalue_style,
    rows = 2:(nrow(sheet_genes) + 1),
    cols = 3,
    gridExpand = TRUE,
    stack = TRUE
  )
}

write_plain_sheet(
  wb = wb,
  sheet_name = "ROC_DeLong",
  data = sheet_roc,
  widths = c(
    14, 10, 13, 11, 19, 17,
    rep(15, 9),
    rep(20, 10)
  )
)

if (nrow(sheet_roc) > 0) {
  addStyle(
    wb,
    "ROC_DeLong",
    integer_style,
    rows = 2:(nrow(sheet_roc) + 1),
    cols = 2:6,
    gridExpand = TRUE,
    stack = TRUE
  )

  addStyle(
    wb,
    "ROC_DeLong",
    decimal_style,
    rows = 2:(nrow(sheet_roc) + 1),
    cols = c(7:18, 21:23),
    gridExpand = TRUE,
    stack = TRUE
  )

  addStyle(
    wb,
    "ROC_DeLong",
    pvalue_style,
    rows = 2:(nrow(sheet_roc) + 1),
    cols = c(19, 20, 24, 25),
    gridExpand = TRUE,
    stack = TRUE
  )
}

write_plain_sheet(
  wb = wb,
  sheet_name = "Incremental_models",
  data = sheet_incremental,
  widths = c(14, 10, 23, 14, 14, 18, 28, 15, 24, 14)
)

if (nrow(sheet_incremental) > 0) {
  addStyle(
    wb,
    "Incremental_models",
    integer_style,
    rows = 2:(nrow(sheet_incremental) + 1),
    cols = 2,
    gridExpand = TRUE,
    stack = TRUE
  )

  addStyle(
    wb,
    "Incremental_models",
    decimal_style,
    rows = 2:(nrow(sheet_incremental) + 1),
    cols = c(3:5, 8:10),
    gridExpand = TRUE,
    stack = TRUE
  )

  addStyle(
    wb,
    "Incremental_models",
    pvalue_style,
    rows = 2:(nrow(sheet_incremental) + 1),
    cols = c(6, 7),
    gridExpand = TRUE,
    stack = TRUE
  )
}

saveWorkbook(
  wb,
  file = OUTPUT_XLSX,
  overwrite = TRUE
)

if (!file.exists(OUTPUT_XLSX)) {
  stop("The XLSX workbook was not created:\n", OUTPUT_XLSX)
}
