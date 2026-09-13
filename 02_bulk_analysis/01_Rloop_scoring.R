# R-loop scoring and phenotype associations

PROJECT_DIR <- "."

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(GSEABase)
  library(GSVA)
  library(ggpubr)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(data.table)
})
select <- dplyr::select
filter <- dplyr::filter

setwd(file.path(PROJECT_DIR, "results/bulk/GSE96058"))

expr_data <- read.csv("GSE96058_exp.csv",
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      check.names = FALSE,
                      row.names = 1)
expr_mat <- as.matrix(expr_data)
mode(expr_mat) <- "numeric"

geneSets <- getGmt("01_Rloop_regulators.gmt", geneIdType = SymbolIdentifier())

gene_sets_list <- list()
for (i in seq_along(geneSets)) {
  gs <- geneSets[[i]]
  gene_sets_list[[setName(gs)]] <- geneIds(gs)}

for (pathway_name in names(gene_sets_list)) {
  genes   <- gene_sets_list[[pathway_name]]
  matched <- intersect(genes, rownames(expr_mat))
}

ssgsea_param <- ssgseaParam(
  exprData  = expr_mat,
  geneSets  = gene_sets_list,
  alpha     = 0.25,
  normalize = FALSE
)
gsvaResult <- gsva(ssgsea_param, verbose = FALSE)

normalize <- function(x) { (x - min(x)) / (max(x) - min(x)) }
gsvaResult <- t(apply(gsvaResult, 1, normalize))

gsvaOut <- as.data.frame(gsvaResult)
colname_row <- as.data.frame(t(colnames(gsvaOut)), stringsAsFactors = FALSE)
colnames(colname_row) <- colnames(gsvaOut)
gsvaOut <- rbind(colname_row, gsvaOut)
rownames(gsvaOut)[1] <- "ID"

write.csv(gsvaOut, file = "GSE96058_ssGSEA.csv", quote = FALSE, row.names = TRUE)

ssgsea_result <- fread("GSE96058_ssGSEA.csv", data.table = FALSE, skip = 1)
rownames(ssgsea_result) <- ssgsea_result[, 1]
ssgsea_result <- ssgsea_result[, -1]

ssgsea_t <- as.data.frame(t(ssgsea_result))
ssgsea_t$Sample <- rownames(ssgsea_t)

score_col    <- "Rloop_regulators"
median_score <- median(ssgsea_t[[score_col]], na.rm = TRUE)
group_data <- data.frame(
  Sample            = rownames(ssgsea_t),
  Rloop_ssGSEA_Score = ssgsea_t[[score_col]],
  Rloop_Group       = ifelse(ssgsea_t[[score_col]] >= median_score, "High", "Low"),
  stringsAsFactors  = FALSE)
group_data$Rloop_Group <- factor(group_data$Rloop_Group, levels = c("Low", "High"))

clini_data <- read.csv("GSE96058_clinisurv.csv",
                       header = T,
                       stringsAsFactors = FALSE,
                       check.names = FALSE,
                       row.names = 1)
clini_data$sample<-rownames(clini_data)
rownames(clini_data)<-NULL

colnames(clini_data)[13] <- "sample"

clinical_updated <- clini_data %>%
  left_join(group_data, by = c("sample" = "Sample")) %>%
  na.omit()
write.csv(clinical_updated, "GSE96058_HRp_HERn_clinical.csv", row.names = FALSE)

valid_samples <- clinical_updated$sample

common_samples <- intersect(colnames(expr_mat), valid_samples)

common_samples <- clinical_updated$sample[clinical_updated$sample %in% common_samples]

expr_mat_common <- expr_mat[, common_samples, drop = FALSE]

group_map <- clinical_updated %>%
  dplyr::select(sample, Rloop_Group)
group_vec <- group_map$Rloop_Group
names(group_vec) <- group_map$sample

new_sample_names <- paste0(common_samples, "_", group_vec[common_samples])
colnames(expr_mat_common) <- new_sample_names

expr_out <- as.data.frame(expr_mat_common, check.names = FALSE)

write.csv(expr_out,
          file = "GSE96058_expr_HighLow.csv",
          quote = FALSE)

write.table(expr_out,
            file = "GSE96058_expr_HighLow.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

clinical <- read.csv("GSE96058_HRp_HERn_clinical.csv", stringsAsFactors = FALSE)

clinical <- clinical %>%
  mutate(
    grade     = trimws(grade),
    grade     = as.character(grade),
    PAM50     = as.character(PAM50),
    ki67      = as.character(ki67),
    chemo_tre = as.character(chemo_tre),
    Rloop_Group = as.character(Rloop_Group),
    age       = as.numeric(age),
    Age_Group = case_when(
      is.na(age) ~ NA_character_,
      age < 50   ~ "<50",
      age >= 50  ~ ">=50"
    ),
    Age_Group = factor(Age_Group, levels = c("<50", ">=50"))
  )

clinical_sorted <- clinical %>%
  arrange(Rloop_ssGSEA_Score) %>%
  mutate(sample_order = row_number())

color_maps <- list(
  grade = c(
    "G1" = "#FEE391",
    "G2" = "#FE9929",
    "G3" = "#CC4C02"
  ),
  PAM50 = c(
    "LumA" = "#187d79",
    "LumB" = "#66C2A4"
  ),
  ki67 = c(
    "0" = "#DEEBF7",
    "1" = "#2171B5"
  ),
  chemo_tre = c(
    "0" = "#EFEDF5",
    "1" = "#6A51A3"
  ),
  Rloop_Group = c(
    "High" = "#b43665",
    "Low"  = "#6fa6cf"
  ),
  Age_Group = c(
    "<50"  = "#8DA0CB",
    ">=50" = "#FC8D62"
  )
)

x_max    <- nrow(clinical_sorted)
x_limits <- c(1, x_max)

get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) == 0) return(NULL)
  tmp$grobs[[leg]]
}

p_main <- ggplot(clinical_sorted, aes(x = sample_order, y = Rloop_ssGSEA_Score)) +
  geom_col(fill = "#000033", alpha = 0.9, width = 1) +
  scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(y = "R-loop\nScore") +
  theme_bw() +
  theme(
    axis.title.x       = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.title.y       = element_text(angle = 0, vjust = 0.5, size = 11, face = "bold", color = "black"),
    axis.text.y        = element_text(size = 9, face = "bold", color = "black"),
    axis.ticks.y       = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.margin        = margin(b = 0, l = 10, r = 10, t = 5)
  )

make_bar <- function(df, fill_var, color_map, y_label, show_leg = FALSE) {
  ggplot(df, aes(x = sample_order, y = 1, fill = .data[[fill_var]])) +
    geom_col(width = 1) +
    scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = color_map, name = y_label, na.value = "#CCCCCC") +
    labs(y = y_label) +
    theme_minimal() +
    theme(
      axis.text.x   = element_blank(),
      axis.text.y   = element_blank(),
      axis.title.x  = element_blank(),
      axis.title.y  = element_text(angle = 0, vjust = 0.5, hjust = 1,
                                   size = 10, face = "bold", color = "black",
                                   margin = margin(r = 5)),
      axis.ticks    = element_blank(),
      panel.grid    = element_blank(),
      plot.margin   = margin(t = 0, l = 10, r = 10, b = 0),
      legend.position   = if (show_leg) "bottom" else "none",
      legend.title      = element_text(size = 9, face = "bold", color = "black"),
      legend.text       = element_text(size = 8, color = "black"),
      legend.key.size   = unit(0.5, "lines"),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width  = unit(0.4, "cm")
    )
}

p_age   <- make_bar(clinical_sorted, "Age_Group",   color_maps$Age_Group,   "Age")
p_grade <- make_bar(clinical_sorted, "grade",        color_maps$grade,       "Grade")
p_pam50 <- make_bar(clinical_sorted, "PAM50",        color_maps$PAM50,       "PAM50")
p_ki67  <- make_bar(clinical_sorted, "ki67",         color_maps$ki67,        "Ki67")
p_chemo <- make_bar(clinical_sorted, "chemo_tre",    color_maps$chemo_tre,   "Chemo")
p_group <- make_bar(clinical_sorted, "Rloop_Group",  color_maps$Rloop_Group, "Group")

legend_age   <- get_legend(make_bar(clinical_sorted, "Age_Group",  color_maps$Age_Group,   "Age",   show_leg = TRUE))
legend_grade <- get_legend(make_bar(clinical_sorted, "grade",       color_maps$grade,       "Grade", show_leg = TRUE))
legend_pam50 <- get_legend(make_bar(clinical_sorted, "PAM50",       color_maps$PAM50,       "PAM50", show_leg = TRUE))
legend_ki67  <- get_legend(make_bar(clinical_sorted, "ki67",        color_maps$ki67,        "Ki67",  show_leg = TRUE))
legend_chemo <- get_legend(make_bar(clinical_sorted, "chemo_tre",   color_maps$chemo_tre,   "Chemo", show_leg = TRUE))
legend_group <- get_legend(make_bar(clinical_sorted, "Rloop_Group", color_maps$Rloop_Group, "Group", show_leg = TRUE))

aligned_plots <- align_plots(
  p_main, p_age, p_grade, p_pam50, p_ki67, p_chemo, p_group,
  align = "v",
  axis  = "tblr"
)

p_anno_combined <- plot_grid(
  aligned_plots[[2]],
  aligned_plots[[3]],
  aligned_plots[[4]],
  aligned_plots[[5]],
  aligned_plots[[6]],
  aligned_plots[[7]],
  ncol        = 1,
  rel_heights = rep(0.7, 6)
)

main_with_anno <- plot_grid(
  aligned_plots[[1]],
  p_anno_combined,
  ncol        = 1,
  rel_heights = c(3.5, 2.8)
)

all_legends <- plot_grid(
  legend_age, legend_grade, legend_pam50, legend_ki67, legend_chemo, legend_group,
  nrow       = 1,
  rel_widths = c(1.1, 1.0, 0.9, 0.8, 0.9, 1.1),
  scale      = 0.95
)

final_plot <- plot_grid(
  main_with_anno,
  all_legends,
  ncol        = 1,
  rel_heights = c(6.5, 1)
)

ggsave("GSE96058_Rloop_waterfall.pdf", final_plot, width = 12, height = 5, dpi = 600)

library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

clinical_box <- clinical %>%
  mutate(
    grade = trimws(grade),
    grade = factor(grade, levels = c("G1", "G2", "G3")),
    PAM50 = factor(PAM50, levels = c("LumA", "LumB")),
    ki67  = factor(ki67, levels = c("0", "1"), labels = c("Low", "High")),
    chemo_tre = factor(chemo_tre, levels = c("0", "1"), labels = c("No", "Yes")),
    age   = as.numeric(age),
    Age_Group = factor(
      case_when(is.na(age) ~ NA_character_, age < 50 ~ "<50", TRUE ~ ">=50"),
      levels = c("<50", ">=50")
    )
  )

plot_violin <- function(data, x_var, x_lab, fill_colors) {

  n_label <- data %>%
    filter(!is.na(.data[[x_var]])) %>%
    group_by(.data[[x_var]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0(as.character(.data[[x_var]]), "\n(n=", n, ")"))
  label_vec <- setNames(n_label$label, n_label[[x_var]])

  data_for_stat <- data %>%
    filter(!is.na(.data[[x_var]])) %>%
    group_by(.data[[x_var]]) %>%
    filter(n() >= 3) %>%
    ungroup() %>%
    droplevels()

  group_levels <- levels(data_for_stat[[x_var]])
  if (length(group_levels) < 2) return(NULL)

  pw     <- pairwise.wilcox.test(
    data_for_stat$Rloop_ssGSEA_Score,
    data_for_stat[[x_var]],
    p.adjust.method = "BH"
  )
  pw_mat <- pw$p.value
  idx    <- which(!is.na(pw_mat), arr.ind = TRUE)

  pw_df  <- data.frame(
    group1 = colnames(pw_mat)[idx[, 2]],
    group2 = rownames(pw_mat)[idx[, 1]],
    p.adj  = pw_mat[idx],
    stringsAsFactors = FALSE
  ) %>%
    mutate(p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ))

  y_max  <- max(data_for_stat$Rloop_ssGSEA_Score, na.rm = TRUE)
  y_step <- diff(range(data_for_stat$Rloop_ssGSEA_Score, na.rm = TRUE)) * 0.12
  pw_df$y.position <- y_max + y_step * seq_len(nrow(pw_df))

  ggplot(data %>% filter(!is.na(.data[[x_var]])),
         aes(x = .data[[x_var]], y = Rloop_ssGSEA_Score, fill = .data[[x_var]])) +
    geom_violin(trim = FALSE, alpha = 0.75, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.15, fill = "white", color = "black",
                 outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.1, size = 0.5, alpha = 0.35, color = "grey30") +
    scale_fill_manual(values = fill_colors) +
    scale_x_discrete(labels = label_vec) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = x_lab, y = "R-loop ssGSEA Score") +
    theme_bw() +
    theme(
      aspect.ratio       = 1,
      axis.title         = element_text(size = 10, face = "bold", color = "black"),
      axis.text.x        = element_text(size = 8.5, color = "black"),
      axis.text.y        = element_text(size = 9,   color = "black"),
      axis.ticks         = element_line(color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position    = "none"
    ) +
    stat_pvalue_manual(pw_df, label = "p.signif", tip.length = 0.01, size = 4)
}

p_vio_grade <- plot_violin(clinical_box, "grade", "Grade",
                           c("G1" = "#efb421", "G2" = "#e58027", "G3" = "#b43665"))

p_vio_pam50 <- plot_violin(clinical_box, "PAM50", "PAM50",
                           c("LumA" = "#D95F02", "LumB" = "#1B9E77"))

p_vio_ki67 <- plot_violin(clinical_box, "ki67", "Ki67",
                          c("Low" = "#DEEBF7", "High" = "#2171B5"))

p_vio_chemo <- plot_violin(clinical_box, "chemo_tre", "Chemotherapy",
                           c("No" = "#EFEDF5", "Yes" = "#6A51A3"))

p_vio_age <- plot_violin(clinical_box, "Age_Group", "Age Group",
                         c("<50" = "#8DA0CB", ">=50" = "#FC8D62"))

vio_combined <- plot_grid(
  p_vio_grade, p_vio_pam50, p_vio_ki67, p_vio_chemo, p_vio_age,
  ncol = 5, nrow = 1,
  labels = c("A", "B", "C", "D", "E"), label_size = 12
)

ggsave("GSE96058_Rloop_violin.pdf",
       vio_combined, width = 22, height = 5, dpi = 600)

expr_data <- read.table("GSE96058_expr_HighLow.txt", header = TRUE, sep = "\t", check.names = FALSE)

colnames(expr_data) <- sub("_.*", "", colnames(expr_data))
rownames(expr_data)<-expr_data[,1]
expr_data<-expr_data[,-1]

write.table(expr_data, file = "GSE96058_expr_NO_HighLow.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

WORK_DIR <- file.path(PROJECT_DIR, "results/bulk/GSE96058")

EXP_FILE <- "GSE96058_expr_NO_HighLow.txt"
RLOOP_SCORE_FILE <- "GSE96058_ssGSEA.csv"
CANCER_SEA_GMT <- "CancerSEA.gmt"

OUTPUT_PREFIX <- "Rloop_CancerSEA"

POINT_SIZE <- 1.5
POINT_ALPHA <- 0.5
LINE_WIDTH <- 0.8
SE_ALPHA <- 0.25

CORR_TEXT_SIZE <- 4
PVAL_TEXT_SIZE <- 3.5
FDR_TEXT_SIZE  <- 3.3

AXIS_TEXT_SIZE  <- 9
AXIS_TITLE_SIZE <- 11
STRIP_TEXT_SIZE <- 10

FIGURE_WIDTH  <- 18
FIGURE_HEIGHT <- 6
DPI <- 300

PHENOTYPE_COLORS <- c(
  "Angiogenesis" =    "#c58cbd",
  "Apoptosis" =       "#1f9953",
  "Cell cycle" =      "#efb421",
  "Differentiation" = "#e58027",
  "DNA damage" =      "#7daada",
  "DNA repair" =      "#4779bd",
  "EMT" =             "#f38185",
  "Hypoxia" =         "#f9ba4e",
  "Inflammation" =    "#1faa9f",
  "Invasion" =        "#c49874",
  "Metastasis" =      "#ee194b",
  "Proliferation" =   "#b43665",
  "Quiescence" =      "#6fa6cf",
  "Stemness" =        "#187d79"
)

suppressPackageStartupMessages({
  library(GSVA)
  library(GSEABase)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

setwd(WORK_DIR)

rloop_raw <- fread(RLOOP_SCORE_FILE, data.table = FALSE)
has_rowname_col <- is.na(suppressWarnings(as.numeric(rloop_raw[3, 1])))

if (has_rowname_col) {
  sample_names <- as.character(rloop_raw[2, -1])
  rloop_values <- as.numeric(rloop_raw[3, -1])
} else {
  sample_names <- as.character(rloop_raw[2, ])
  rloop_values <- as.numeric(rloop_raw[3, ])
}

rloop_scores <- rloop_values
names(rloop_scores) <- sample_names
rloop_scores <- rloop_scores[!is.na(rloop_scores) & !is.na(names(rloop_scores))]

expr_data <- fread(EXP_FILE, data.table = FALSE)
gene_names <- expr_data[, 1]
expr_data <- expr_data[, -1]
expr_matrix <- as.matrix(expr_data)
rownames(expr_matrix) <- gene_names
mode(expr_matrix) <- "numeric"

gene_sets_gmt <- getGmt(CANCER_SEA_GMT, geneIdType = SymbolIdentifier())
gene_sets_list <- list()
for (i in seq_along(gene_sets_gmt)) {
  gs <- gene_sets_gmt[[i]]
  gene_sets_list[[setName(gs)]] <- geneIds(gs)
}

zscore_param <- zscoreParam(exprData = expr_matrix, geneSets = gene_sets_list)
gsva_result <- gsva(zscore_param, verbose = FALSE)
phenotype_scores <- t(scale(t(gsva_result)))

common_samples <- intersect(names(rloop_scores), colnames(phenotype_scores))

if (length(common_samples) == 0) {
  rloop_names_clean <- trimws(gsub("^\"|\"$", "", names(rloop_scores)))
  names(rloop_scores) <- rloop_names_clean
  pheno_names_clean <- trimws(colnames(phenotype_scores))
  colnames(phenotype_scores) <- pheno_names_clean
  common_samples <- intersect(names(rloop_scores), colnames(phenotype_scores))
}

if (length(common_samples) < 10) {
  stop("sample < 10, checksample.")
}

rloop_vec <- as.numeric(scale(rloop_scores[common_samples]))
names(rloop_vec) <- common_samples

pheno_mat <- phenotype_scores[, common_samples]

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
  cor_test <- cor.test(rloop_vec, pheno_vec, method = "spearman")
  cor_results$Correlation[i] <- cor_test$estimate
  cor_results$Pvalue[i]      <- cor_test$p.value
}

cor_results$FDR <- p.adjust(cor_results$Pvalue, method = "BH")

write.csv(cor_results, paste0(OUTPUT_PREFIX, "_correlation.csv"),
          row.names = FALSE, quote = FALSE)

rloop_range <- range(rloop_vec, na.rm = TRUE)
pheno_range <- range(as.vector(pheno_mat), na.rm = TRUE)

y_limits <- c(floor(rloop_range[1] - 0.2), ceiling(rloop_range[2] + 0.2))
x_limits <- c(floor(pheno_range[1] - 0.2), ceiling(pheno_range[2] + 0.2))

y_breaks <- seq(y_limits[1], y_limits[2], length.out = 5)
x_breaks <- seq(x_limits[1], x_limits[2], length.out = 5)

make_scatter <- function(pheno_name, pheno_vec, rloop_vec,
                         cor_val, p_val, fdr_val, color,
                         is_left_col, is_bottom_row,
                         x_limits, y_limits, x_breaks, y_breaks) {

  df <- data.frame(Phenotype = pheno_vec, Rloop = rloop_vec)

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
  fdr_text <- format_fdr(fdr_val)
  cor_text <- sprintf("R=%.2f", cor_val)

  p <- ggplot(df, aes(x = Phenotype, y = Rloop)) +
    geom_point(color = color, alpha = POINT_ALPHA, size = POINT_SIZE) +
    geom_smooth(method = "lm", se = TRUE,
                color = color,
                fill  = color,
                alpha = SE_ALPHA,
                linewidth = LINE_WIDTH) +

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
        hjust = 0.5,
        face  = "bold",
        size  = STRIP_TEXT_SIZE,
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

  is_left_col   <- (i - 1) %% 7 == 0
  is_bottom_row <- i > 7

  plot_list[[i]] <- make_scatter(
    pname, pvec, rloop_vec,
    corval, pval, fdrval, col,
    is_left_col, is_bottom_row,
    x_limits, y_limits, x_breaks, y_breaks
  )
}

combined <- wrap_plots(plot_list, nrow = 2, ncol = 7)

final_plot <- combined +
  plot_annotation(
    title = NULL,
    caption = "Phenotype Score (z-score)",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = AXIS_TITLE_SIZE,
                                  face = "bold",
                                  margin = margin(t = 10)),
      plot.margin  = margin(5, 5, 5, 20)
    )
  )

ggsave(paste0(OUTPUT_PREFIX, "_scatter.pdf"), final_plot,
       width = FIGURE_WIDTH, height = FIGURE_HEIGHT, dpi = DPI)
