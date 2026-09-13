# Neoadjuvant response analysis

PROJECT_DIR <- "."

gc()
set.seed(1234)

INPUT_DIR <- file.path(PROJECT_DIR, "data/processed/NAC")
OUTPUT_DIR <- file.path(PROJECT_DIR, "results/NAC")

RL_SIG65_FILE <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")
RLOOP_GMT_FILE <- file.path(PROJECT_DIR, "results/bulk/GSE96058/01_Rloop_regulators.gmt")
HALLMARK_GMT_FILE <- Sys.getenv(
  "HALLMARK_GMT",
  unset = file.path(PROJECT_DIR, "data/gene_sets/h.all.v2026.1.Hs.symbols.gmt")
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cran_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "cowplot",
  "pROC", "metafor", "brglm2", "openxlsx", "scales", "stringr"
)
bioc_packages <- c("GSVA", "GSEABase", "limma", "clusterProfiler")

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_cran) > 0) install.packages(missing_cran, dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(pROC)
  library(metafor)
  library(brglm2)
  library(openxlsx)
  library(scales)
  library(stringr)
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(clusterProfiler)
})

select <- dplyr::select
filter <- dplyr::filter

assert_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ":\n", path)
}
assert_file(RL_SIG65_FILE, "RL-Sig65 file")
assert_file(RLOOP_GMT_FILE, "R-loop GMT")
assert_file(HALLMARK_GMT_FILE, "Hallmark GMT")

base_theme <- theme_classic(base_family = "Arial")

safe_zscore <- function(x) {
  x <- as.numeric(x)
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

format_p <- function(p) {
  p <- as.numeric(p)
  out <- rep("P = NA", length(p))
  ok <- !is.na(p) & is.finite(p)
  out[ok & p < 0.001] <- "P < 0.001"
  idx <- ok & p >= 0.001
  out[idx] <- paste0("P = ", formatC(p[idx], format = "f", digits = 3))
  out
}

format_fdr <- function(p) {
  p <- as.numeric(p)
  out <- rep("FDR = NA", length(p))
  ok <- !is.na(p) & is.finite(p)
  out[ok & p < 0.001] <- "FDR < 0.001"
  idx <- ok & p >= 0.001
  out[idx] <- paste0("FDR = ", formatC(p[idx], format = "f", digits = 3))
  out
}

fdr_to_star <- function(fdr) {
  sapply(fdr, function(x) {
    if (is.na(x)) return("ns")
    if (x < 1e-4) return("****")
    if (x < 1e-3) return("***")
    if (x < 1e-2) return("**")
    if (x < 5e-2) return("*")
    "ns"
  })
}

extract_gse <- function(x) {
  hit <- regmatches(x, regexpr("GSE[0-9]+", x, ignore.case = TRUE))
  hit <- toupper(hit)
  ifelse(nchar(hit) > 0, hit, NA_character_)
}

extract_response <- function(sample_names, file_name = "") {
  response <- rep(NA_character_, length(sample_names))
  response[grepl("([_.-])pCR$", sample_names, ignore.case = TRUE)] <- "pCR"
  response[grepl("([_.-])RD$", sample_names, ignore.case = TRUE)] <- "RD"
  if (all(is.na(response))) {
    if (grepl("pCR", file_name, ignore.case = TRUE)) response[] <- "pCR"
    if (grepl("([_.-])RD", file_name, ignore.case = TRUE)) response[] <- "RD"
  }
  response
}

clean_gene_vector <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & x != ""])
}

read_gmt_as_list <- function(path) {
  gsc <- GSEABase::getGmt(path, geneIdType = GSEABase::SymbolIdentifier())
  out <- lapply(gsc, GSEABase::geneIds)
  names(out) <- vapply(gsc, GSEABase::setName, character(1))
  out
}

read_original_rloop_genes <- function(path) {
  sets <- read_gmt_as_list(path)
  clean_gene_vector(unlist(sets, use.names = FALSE))
}

clean_hallmark_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("Dna", "DNA", stringr::str_to_title(tolower(x)))
  x <- gsub("E2F", "E2F", x)
  x <- gsub("G2M", "G2M", x)
  x
}

collapse_duplicate_genes <- function(mat) {
  genes <- rownames(mat)
  if (!anyDuplicated(genes)) return(mat)
  summed <- rowsum(mat, group = genes, reorder = FALSE, na.rm = TRUE)
  counts <- as.numeric(table(factor(genes, levels = rownames(summed))))
  summed / counts
}

read_expression_file <- function(path) {
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  genes <- trimws(as.character(dt[[1]]))
  sample_names <- colnames(dt)[-1]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  rownames(mat) <- genes
  keep_gene <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep_gene, , drop = FALSE]
  na_fraction <- rowMeans(is.na(mat))
  mat <- mat[na_fraction <= 0.5, , drop = FALSE]
  for (i in seq_len(nrow(mat))) {
    miss <- is.na(mat[i, ])
    if (any(miss)) {
      repl <- median(mat[i, ], na.rm = TRUE)
      if (!is.finite(repl)) repl <- 0
      mat[i, miss] <- repl
    }
  }
  mat <- collapse_duplicate_genes(mat)
  cohort <- extract_gse(basename(path))
  if (is.na(cohort)) {
    cohort_candidates <- unique(na.omit(extract_gse(sample_names)))
    if (length(cohort_candidates) == 1) cohort <- cohort_candidates else stop("Cannot determine cohort for: ", path)
  }
  has_prefix <- grepl(cohort, sample_names, ignore.case = TRUE)
  sample_names[!has_prefix] <- paste0(cohort, "__", sample_names[!has_prefix])
  sample_names <- make.unique(sample_names)
  colnames(mat) <- sample_names
  response <- extract_response(sample_names, basename(path))
  if (any(is.na(response))) stop("Cannot determine response labels in ", basename(path))
  list(path = path, cohort = cohort, expression = mat, response = setNames(response, sample_names))
}

combine_same_cohort_files <- function(file_objects) {
  cohort_names <- unique(vapply(file_objects, `[[`, character(1), "cohort"))
  result <- list()
  for (cohort in cohort_names) {
    cohort_objects <- file_objects[vapply(file_objects, function(x) x$cohort == cohort, logical(1))]
    if (length(cohort_objects) == 1) {
      result[[cohort]] <- cohort_objects[[1]]
      next
    }
    common_genes <- Reduce(intersect, lapply(cohort_objects, function(x) rownames(x$expression)))
    mats <- lapply(cohort_objects, function(x) x$expression[common_genes, , drop = FALSE])
    combined_mat <- do.call(cbind, mats)
    combined_response <- unlist(lapply(cohort_objects, `[[`, "response"), use.names = TRUE)
    if (anyDuplicated(colnames(combined_mat))) {
      new_names <- make.unique(colnames(combined_mat))
      names(combined_response) <- new_names
      colnames(combined_mat) <- new_names
    }
    result[[cohort]] <- list(
      path = paste(vapply(cohort_objects, `[[`, character(1), "path"), collapse = " | "),
      cohort = cohort,
      expression = combined_mat,
      response = combined_response[colnames(combined_mat)]
    )
  }
  result
}

calculate_ssgsea <- function(expr_mat, gene_sets) {
  gene_sets <- lapply(gene_sets, function(x) intersect(x, rownames(expr_mat)))
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) >= 5]
  if (exists("ssgseaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param <- GSVA::ssgseaParam(exprData = expr_mat, geneSets = gene_sets, alpha = 0.25, normalize = FALSE)
    result <- GSVA::gsva(param, verbose = FALSE)
  } else {
    result <- GSVA::gsva(expr = expr_mat, gset.idx.list = gene_sets, method = "ssgsea", ssgsea.norm = FALSE, verbose = FALSE)
  }
  as.matrix(result)
}

fit_bias_reduced_logistic <- function(data, formula_object) {
  glm(formula_object, data = data, family = binomial(), method = brglm2::brglmFit)
}

extract_model_term <- function(fit, term_name) {
  coef_tab <- summary(fit)$coefficients
  if (!term_name %in% rownames(coef_tab)) {
    return(data.frame(beta = NA_real_, SE = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P_value = NA_real_))
  }
  beta <- unname(coef_tab[term_name, "Estimate"])
  se <- unname(coef_tab[term_name, "Std. Error"])
  p <- unname(coef_tab[term_name, "Pr(>|z|)"])
  data.frame(beta = beta, SE = se, OR = exp(beta), CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 1.96 * se), P_value = p)
}

fit_random_effects_meta <- function(rows, label) {
  usable <- rows %>% filter(is.finite(beta), is.finite(SE), SE > 0)
  if (nrow(usable) < 2) {
    return(data.frame(Cohort = "Random-effects pooled", Effect = label, N = sum(usable$N), N_pCR = sum(usable$N_pCR), N_RD = sum(usable$N_RD), beta = NA_real_, SE = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P_value = NA_real_, Q = NA_real_, Q_P = NA_real_, I2 = NA_real_))
  }
  meta_fit <- metafor::rma.uni(yi = usable$beta, sei = usable$SE, method = "REML", slab = usable$Cohort)
  data.frame(Cohort = "Random-effects pooled", Effect = label, N = sum(usable$N), N_pCR = sum(usable$N_pCR), N_RD = sum(usable$N_RD), beta = as.numeric(meta_fit$b), SE = meta_fit$se, OR = exp(as.numeric(meta_fit$b)), CI_low = exp(meta_fit$ci.lb), CI_high = exp(meta_fit$ci.ub), P_value = meta_fit$pval, Q = meta_fit$QE, Q_P = meta_fit$QEp, I2 = meta_fit$I2)
}

make_roc_object <- function(y, predictor) {
  pROC::roc(response = y, predictor = predictor, levels = c(0, 1), direction = "<", quiet = TRUE, auc = TRUE)
}

summarize_roc <- function(roc_object, cohort, model_name, n, n_pcr, n_rd) {
  auc_value <- as.numeric(pROC::auc(roc_object))
  auc_ci <- as.numeric(pROC::ci.auc(roc_object, conf.level = 0.95, method = "delong"))
  data.frame(Cohort = cohort, Model = model_name, N = n, N_pCR = n_pcr, N_RD = n_rd, AUC = auc_value, AUC_CI_low = auc_ci[1], AUC_CI_high = auc_ci[3])
}

roc_to_curve_df <- function(roc_object, cohort, model_name) {
  xy <- pROC::coords(roc_object, x = "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  data.frame(Cohort = cohort, Model = model_name, Specificity = as.numeric(xy$specificity), Sensitivity = as.numeric(xy$sensitivity))
}

calculate_pr_curve <- function(y, predictor, cohort, model_name) {
  y <- as.integer(y)
  predictor <- as.numeric(predictor)
  ord <- order(predictor, decreasing = TRUE)
  y_ord <- y[ord]
  total_pos <- sum(y_ord == 1)
  total_neg <- sum(y_ord == 0)
  if (total_pos == 0 || total_neg == 0) {
    return(list(curve = data.frame(), summary = data.frame(Cohort = cohort, Model = model_name, N = length(y), N_pCR = total_pos, N_RD = total_neg, pCR_prevalence = mean(y), AUPRC = NA_real_)))
  }
  tp <- cumsum(y_ord == 1)
  fp <- cumsum(y_ord == 0)
  recall <- tp / total_pos
  precision <- tp / pmax(tp + fp, 1)
  recall <- c(0, recall)
  precision <- c(1, precision)
  auprc <- sum(diff(recall) * (head(precision, -1) + tail(precision, -1)) / 2, na.rm = TRUE)
  list(curve = data.frame(Cohort = cohort, Model = model_name, Recall = recall, Precision = precision), summary = data.frame(Cohort = cohort, Model = model_name, N = length(y), N_pCR = total_pos, N_RD = total_neg, pCR_prevalence = mean(y), AUPRC = auprc))
}

save_pdf <- function(filename, plot, width, height) {
  output_path <- file.path(OUTPUT_DIR, filename)
  if (capabilities("cairo")) {
    ggsave(output_path, plot = plot, width = width, height = height, device = cairo_pdf, units = "in")
  } else {
    ggsave(output_path, plot = plot, width = width, height = height, device = pdf, units = "in")
  }
}

rlsig65_genes <- clean_gene_vector(readLines(RL_SIG65_FILE, warn = FALSE, encoding = "UTF-8"))
rloopbase_genes <- read_original_rloop_genes(RLOOP_GMT_FILE)
hallmark_sets_full <- read_gmt_as_list(HALLMARK_GMT_FILE)
names(hallmark_sets_full) <- gsub("^HALLMARK_", "", names(hallmark_sets_full))
hallmark_term2gene <- stack(hallmark_sets_full)
colnames(hallmark_term2gene) <- c("gene_symbol", "gs_name")
hallmark_term2gene <- hallmark_term2gene[, c("gs_name", "gene_symbol")]

input_files <- list.files(INPUT_DIR, pattern = "\\.txt$", full.names = TRUE)
if (length(input_files) == 0) stop("No txt files found in INPUT_DIR")
file_objects <- lapply(input_files, read_expression_file)
cohort_objects <- combine_same_cohort_files(file_objects)
cohort_order <- names(cohort_objects)
cohort_order <- cohort_order[order(as.numeric(gsub("GSE", "", cohort_order)))]

sample_score_list <- list()
cohort_summary_list <- list()
hallmark_score_list <- list()
for (cohort in cohort_order) {
  obj <- cohort_objects[[cohort]]
  expr <- obj$expression
  if (!"MKI67" %in% rownames(expr)) stop("MKI67 absent in ", cohort)
  main_sets <- list(RL_Sig65 = intersect(rlsig65_genes, rownames(expr)), Original_Rloop_related_score = intersect(rloopbase_genes, rownames(expr)))
  main_scores <- calculate_ssgsea(expr, main_sets)
  hallmark_sets_cohort <- lapply(hallmark_sets_full, function(x) intersect(x, rownames(expr)))
  hallmark_sets_cohort <- hallmark_sets_cohort[vapply(hallmark_sets_cohort, length, integer(1)) >= 10]
  hallmark_scores <- calculate_ssgsea(expr, hallmark_sets_cohort)
  sample_names <- intersect(names(obj$response), colnames(main_scores))
  cohort_df <- data.frame(Sample = sample_names, Cohort = cohort, Response = factor(obj$response[sample_names], levels = c("RD", "pCR")), Response_Binary = ifelse(obj$response[sample_names] == "pCR", 1L, 0L), RD_Binary = ifelse(obj$response[sample_names] == "RD", 1L, 0L), RL_Sig65_raw = as.numeric(main_scores["RL_Sig65", sample_names]), Original_Rloop_raw = as.numeric(main_scores["Original_Rloop_related_score", sample_names]), MKI67_raw = as.numeric(expr["MKI67", sample_names]))
  cohort_df$RL_Sig65_z <- safe_zscore(cohort_df$RL_Sig65_raw)
  cohort_df$Original_Rloop_z <- safe_zscore(cohort_df$Original_Rloop_raw)
  cohort_df$MKI67_z <- safe_zscore(cohort_df$MKI67_raw)
  cohort_df <- cohort_df[complete.cases(cohort_df), , drop = FALSE]
  sample_score_list[[cohort]] <- cohort_df
  hallmark_sample_df <- as.data.frame(t(hallmark_scores[, cohort_df$Sample, drop = FALSE]), check.names = FALSE)
  hallmark_sample_df$Sample <- rownames(hallmark_sample_df)
  hallmark_sample_df$Cohort <- cohort
  hallmark_score_list[[cohort]] <- hallmark_sample_df
  cohort_summary_list[[cohort]] <- data.frame(Cohort = cohort, Input_file = obj$path, N = nrow(cohort_df), N_pCR = sum(cohort_df$Response == "pCR"), N_RD = sum(cohort_df$Response == "RD"), pCR_rate = mean(cohort_df$Response == "pCR"), Expression_genes = nrow(expr), RL_Sig65_genes_matched = length(main_sets$RL_Sig65), RL_Sig65_genes_total = length(rlsig65_genes), Original_Rloop_genes_matched = length(main_sets$Original_Rloop_related_score), Original_Rloop_genes_total = length(rloopbase_genes), Hallmark_sets_scored = nrow(hallmark_scores))
}
sample_scores <- bind_rows(sample_score_list) %>% mutate(Cohort = factor(Cohort, levels = cohort_order), Response = factor(Response, levels = c("RD", "pCR")))
cohort_summary <- bind_rows(cohort_summary_list)
hallmark_sample_scores <- bind_rows(hallmark_score_list)

univariable_definitions <- data.frame(Effect = c("RL-Sig65", "Original R-loop-related score", "MKI67 expression"), Column = c("RL_Sig65_z", "Original_Rloop_z", "MKI67_z"))
association_rows <- list(); adjusted_rows <- list()
for (cohort in cohort_order) {
  cohort_df <- sample_scores %>% filter(Cohort == cohort) %>% mutate(y = Response_Binary)
  for (j in seq_len(nrow(univariable_definitions))) {
    effect_label <- univariable_definitions$Effect[j]
    score_col <- univariable_definitions$Column[j]
    model_data <- data.frame(y = cohort_df$y, score = cohort_df[[score_col]])
    fit <- fit_bias_reduced_logistic(model_data, y ~ score)
    row <- extract_model_term(fit, "score")
    row$Cohort <- cohort; row$Effect <- effect_label; row$N <- nrow(cohort_df); row$N_pCR <- sum(cohort_df$Response == "pCR"); row$N_RD <- sum(cohort_df$Response == "RD")
    association_rows[[paste(cohort, effect_label)]] <- row %>% select(Cohort, Effect, N, N_pCR, N_RD, beta, SE, OR, CI_low, CI_high, P_value)
  }
  adjusted_data <- data.frame(y = cohort_df$y, RL_Sig65_z = cohort_df$RL_Sig65_z, MKI67_z = cohort_df$MKI67_z)
  adjusted_fit <- fit_bias_reduced_logistic(adjusted_data, y ~ MKI67_z + RL_Sig65_z)
  adjusted_row <- extract_model_term(adjusted_fit, "RL_Sig65_z")
  adjusted_row$Cohort <- cohort; adjusted_row$Effect <- "RL-Sig65 adjusted for MKI67"; adjusted_row$N <- nrow(cohort_df); adjusted_row$N_pCR <- sum(cohort_df$Response == "pCR"); adjusted_row$N_RD <- sum(cohort_df$Response == "RD")
  adjusted_rows[[cohort]] <- adjusted_row %>% select(Cohort, Effect, N, N_pCR, N_RD, beta, SE, OR, CI_low, CI_high, P_value)
}
association_results <- bind_rows(association_rows)
adjusted_results <- bind_rows(adjusted_rows)
association_meta <- bind_rows(
  fit_random_effects_meta(association_results %>% filter(Effect == "RL-Sig65"), "RL-Sig65"),
  fit_random_effects_meta(association_results %>% filter(Effect == "Original R-loop-related score"), "Original R-loop-related score"),
  fit_random_effects_meta(association_results %>% filter(Effect == "MKI67 expression"), "MKI67 expression")
)
adjusted_meta <- fit_random_effects_meta(adjusted_results, "RL-Sig65 adjusted for MKI67")
association_all <- bind_rows(association_results, association_meta)
adjusted_all <- bind_rows(adjusted_results, adjusted_meta)

score_model_definitions <- data.frame(Model = c("RL-Sig65", "Original R-loop-related score", "MKI67 expression"), Column = c("RL_Sig65_z", "Original_Rloop_z", "MKI67_z"))
cohort_roc_summary_list <- list(); cohort_roc_curve_list <- list(); cohort_pr_summary_list <- list(); cohort_pr_curve_list <- list()
for (cohort in cohort_order) {
  cohort_df <- sample_scores %>% filter(Cohort == cohort) %>% arrange(Sample)
  y <- cohort_df$Response_Binary
  for (j in seq_len(nrow(score_model_definitions))) {
    model_label <- score_model_definitions$Model[j]
    score_col <- score_model_definitions$Column[j]
    roc_obj <- make_roc_object(y, cohort_df[[score_col]])
    cohort_roc_summary_list[[paste(cohort, model_label)]] <- summarize_roc(roc_obj, cohort, model_label, nrow(cohort_df), sum(y == 1), sum(y == 0))
    cohort_roc_curve_list[[paste(cohort, model_label)]] <- roc_to_curve_df(roc_obj, cohort, model_label)
    pr_result <- calculate_pr_curve(y, cohort_df[[score_col]], cohort, model_label)
    cohort_pr_summary_list[[paste(cohort, model_label)]] <- pr_result$summary
    cohort_pr_curve_list[[paste(cohort, model_label)]] <- pr_result$curve
  }
}
cohort_roc_summary <- bind_rows(cohort_roc_summary_list)
cohort_roc_curves <- bind_rows(cohort_roc_curve_list)
cohort_pr_summary <- bind_rows(cohort_pr_summary_list)
cohort_pr_curves <- bind_rows(cohort_pr_curve_list)

loco_model_definitions <- list("RL-Sig65" = y ~ RL_Sig65_z, "Original R-loop-related score" = y ~ Original_Rloop_z, "MKI67 expression" = y ~ MKI67_z, "RL-Sig65 + MKI67" = y ~ RL_Sig65_z + MKI67_z)
loco_prediction_list <- list()
for (held_out in cohort_order) {
  train_df <- sample_scores %>% filter(Cohort != held_out) %>% mutate(y = Response_Binary)
  test_df <- sample_scores %>% filter(Cohort == held_out) %>% mutate(y = Response_Binary)
  for (model_name in names(loco_model_definitions)) {
    fit <- fit_bias_reduced_logistic(train_df, loco_model_definitions[[model_name]])
    predicted <- as.numeric(predict(fit, newdata = test_df, type = "response"))
    loco_prediction_list[[paste(held_out, model_name)]] <- data.frame(Sample = test_df$Sample, Cohort = held_out, Response = as.character(test_df$Response), Response_Binary = test_df$Response_Binary, Model = model_name, Predicted_probability = predicted)
  }
}
loco_predictions <- bind_rows(loco_prediction_list) %>% mutate(Cohort = factor(Cohort, levels = cohort_order), Model = factor(Model, levels = c("RL-Sig65", "Original R-loop-related score", "MKI67 expression", "RL-Sig65 + MKI67")))

loco_roc_summary_list <- list(); loco_roc_curve_list <- list(); loco_pr_summary_list <- list(); loco_pr_curve_list <- list(); loco_roc_objects <- list()
for (model_name in levels(loco_predictions$Model)) {
  model_df <- loco_predictions %>% filter(Model == model_name) %>% arrange(Cohort, Sample)
  roc_obj <- make_roc_object(model_df$Response_Binary, model_df$Predicted_probability)
  loco_roc_objects[[model_name]] <- roc_obj
  loco_roc_summary_list[[model_name]] <- summarize_roc(roc_obj, "LOCO combined", model_name, nrow(model_df), sum(model_df$Response_Binary == 1), sum(model_df$Response_Binary == 0))
  loco_roc_curve_list[[model_name]] <- roc_to_curve_df(roc_obj, "LOCO combined", model_name)
  pr_result <- calculate_pr_curve(model_df$Response_Binary, model_df$Predicted_probability, "LOCO combined", model_name)
  loco_pr_summary_list[[model_name]] <- pr_result$summary
  loco_pr_curve_list[[model_name]] <- pr_result$curve
}
loco_roc_summary <- bind_rows(loco_roc_summary_list)
loco_roc_curves <- bind_rows(loco_roc_curve_list)
loco_pr_summary <- bind_rows(loco_pr_summary_list)
loco_pr_curves <- bind_rows(loco_pr_curve_list)

key_delong_comparisons <- list(c("RL-Sig65", "Original R-loop-related score"), c("RL-Sig65", "MKI67 expression"), c("RL-Sig65 + MKI67", "MKI67 expression"), c("RL-Sig65 + MKI67", "RL-Sig65"))
loco_delong_list <- list()
for (comparison in key_delong_comparisons) {
  model_1 <- comparison[1]; model_2 <- comparison[2]
  test <- pROC::roc.test(loco_roc_objects[[model_1]], loco_roc_objects[[model_2]], method = "delong", paired = TRUE)
  ci_diff <- if (!is.null(test$conf.int)) as.numeric(test$conf.int) else c(NA_real_, NA_real_)
  loco_delong_list[[paste(model_1, model_2)]] <- data.frame(Model_1 = model_1, Model_2 = model_2, AUC_1 = as.numeric(pROC::auc(loco_roc_objects[[model_1]])), AUC_2 = as.numeric(pROC::auc(loco_roc_objects[[model_2]])), Delta_AUC_1_minus_2 = as.numeric(pROC::auc(loco_roc_objects[[model_1]]) - pROC::auc(loco_roc_objects[[model_2]])), Delta_AUC_CI_low = ci_diff[1], Delta_AUC_CI_high = ci_diff[2], DeLong_P = test$p.value)
}
loco_delong <- bind_rows(loco_delong_list) %>% mutate(DeLong_BH = p.adjust(DeLong_P, method = "BH"))

residual_de_list <- list(); cohort_gsea_list <- list()
for (cohort in cohort_order) {
  obj <- cohort_objects[[cohort]]
  cohort_df <- sample_scores %>% filter(Cohort == cohort) %>% arrange(Sample)
  expr <- obj$expression[, cohort_df$Sample, drop = FALSE]
  design_data <- data.frame(RL_Sig65_z = cohort_df$RL_Sig65_z, RD_Binary = cohort_df$RD_Binary)
  design <- model.matrix(~ RL_Sig65_z + RD_Binary, data = design_data)
  fit <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit, trend = TRUE)
  de <- limma::topTable(fit, coef = "RD_Binary", number = Inf, sort.by = "none")
  de$Gene <- rownames(de); de$Cohort <- cohort; de$N <- nrow(cohort_df); de$N_pCR <- sum(cohort_df$Response == "pCR"); de$N_RD <- sum(cohort_df$Response == "RD")
  residual_de_list[[cohort]] <- de %>% transmute(Cohort, Gene, N, N_pCR, N_RD, logFC_RD_vs_pCR_adjusted_RL_Sig65 = logFC, AveExpr, t, P_value = P.Value, FDR = adj.P.Val)
  gene_list <- de$t; names(gene_list) <- de$Gene
  gene_list <- gene_list[is.finite(gene_list) & !is.na(names(gene_list)) & names(gene_list) != ""]
  gene_list <- sort(gene_list, decreasing = TRUE); gene_list <- gene_list[!duplicated(names(gene_list))]
  gsea_fit <- suppressWarnings(clusterProfiler::GSEA(geneList = gene_list, TERM2GENE = hallmark_term2gene, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, eps = 0, verbose = FALSE, seed = TRUE))
  cohort_gsea <- as.data.frame(gsea_fit)
  if (nrow(cohort_gsea) > 0) {
    cohort_gsea$Cohort <- cohort
    cohort_gsea_list[[cohort]] <- cohort_gsea %>% transmute(Cohort, Pathway_ID = ID, Pathway = clean_hallmark_name(ID), setSize, enrichmentScore, NES, P_value = pvalue, FDR = p.adjust, qvalue, core_enrichment)
  }
}
residual_de_all <- bind_rows(residual_de_list)
residual_meta <- residual_de_all %>% mutate(P_value_bounded = pmin(pmax(P_value, 1e-300), 1), Signed_Z = sign(logFC_RD_vs_pCR_adjusted_RL_Sig65) * qnorm(P_value_bounded / 2, lower.tail = FALSE), Weight = sqrt(N)) %>% group_by(Gene) %>% summarise(Cohorts_available = n(), Mean_logFC = mean(logFC_RD_vs_pCR_adjusted_RL_Sig65, na.rm = TRUE), Direction_consistency = abs(mean(sign(logFC_RD_vs_pCR_adjusted_RL_Sig65), na.rm = TRUE)), Meta_Z = sum(Weight * Signed_Z, na.rm = TRUE) / sqrt(sum(Weight^2, na.rm = TRUE)), .groups = "drop") %>% filter(Cohorts_available >= 3) %>% mutate(Meta_P = 2 * pnorm(-abs(Meta_Z)), Meta_FDR = p.adjust(Meta_P, method = "BH"), Direction = ifelse(Meta_Z > 0, "Higher in RD", "Higher in pCR"))
meta_gene_list <- residual_meta$Meta_Z; names(meta_gene_list) <- residual_meta$Gene; meta_gene_list <- meta_gene_list[is.finite(meta_gene_list)]; meta_gene_list <- sort(meta_gene_list, decreasing = TRUE); meta_gene_list <- meta_gene_list[!duplicated(names(meta_gene_list))]
meta_gsea_fit <- suppressWarnings(clusterProfiler::GSEA(geneList = meta_gene_list, TERM2GENE = hallmark_term2gene, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, eps = 0, verbose = FALSE, seed = TRUE))
meta_gsea <- as.data.frame(meta_gsea_fit) %>% transmute(Pathway_ID = ID, Pathway = clean_hallmark_name(ID), setSize, enrichmentScore, NES, P_value = pvalue, FDR = p.adjust, qvalue, core_enrichment, Direction = ifelse(NES > 0, "Enriched in RD", "Enriched in pCR")) %>% arrange(FDR, desc(abs(NES)))
cohort_gsea_all <- bind_rows(cohort_gsea_list)
positive_candidates <- meta_gsea %>% filter(NES > 0) %>% arrange(FDR, desc(NES)) %>% slice_head(n = 4)
negative_candidates <- meta_gsea %>% filter(NES < 0) %>% arrange(FDR, NES) %>% slice_head(n = 4)
selected_pathways <- bind_rows(positive_candidates, negative_candidates) %>% distinct(Pathway_ID, .keep_all = TRUE)
selected_pathway_ids <- selected_pathways$Pathway_ID

pathway_association_rows <- list()
for (cohort in cohort_order) {
  cohort_samples <- sample_scores %>% filter(Cohort == cohort) %>% arrange(Sample)
  cohort_hallmark <- hallmark_sample_scores %>% filter(Cohort == cohort)
  rownames(cohort_hallmark) <- cohort_hallmark$Sample
  for (pathway_id in selected_pathway_ids) {
    pathway_col <- gsub("^HALLMARK_", "", pathway_id)
    if (!pathway_col %in% colnames(cohort_hallmark)) next
    pathway_raw <- cohort_hallmark[cohort_samples$Sample, pathway_col]
    model_data <- data.frame(RD = cohort_samples$RD_Binary, Pathway_z = safe_zscore(pathway_raw), RL_Sig65_z = cohort_samples$RL_Sig65_z)
    fit <- fit_bias_reduced_logistic(model_data, RD ~ RL_Sig65_z + Pathway_z)
    row <- extract_model_term(fit, "Pathway_z")
    row$Cohort <- cohort; row$Pathway_ID <- pathway_id; row$Pathway <- clean_hallmark_name(pathway_id); row$N <- nrow(model_data); row$N_RD <- sum(model_data$RD == 1); row$N_pCR <- sum(model_data$RD == 0)
    pathway_association_rows[[paste(cohort, pathway_id)]] <- row %>% select(Cohort, Pathway_ID, Pathway, N, N_RD, N_pCR, beta, SE, OR, CI_low, CI_high, P_value)
  }
}
pathway_association <- bind_rows(pathway_association_rows)
pathway_meta_list <- list()
for (pathway_id in selected_pathway_ids) {
  pathway_rows <- pathway_association %>% filter(Pathway_ID == pathway_id)
  meta_row <- fit_random_effects_meta(pathway_rows %>% transmute(Cohort, Effect = Pathway, N, N_pCR, N_RD, beta, SE, OR, CI_low, CI_high, P_value), clean_hallmark_name(pathway_id))
  meta_row$Pathway_ID <- pathway_id; meta_row$Pathway <- clean_hallmark_name(pathway_id)
  pathway_meta_list[[pathway_id]] <- meta_row
}
pathway_meta <- bind_rows(pathway_meta_list) %>% mutate(FDR = p.adjust(P_value, method = "BH"))

COL_RLSIG65 <- "#B2476B"
COL_RLOOP <- "#3E86B8"
COL_MKI67 <- "#D79A2E"
COL_COMBINED <- "#4B4B4B"
COL_PCR <- "#B2476B"
COL_RD <- "#8FB5CE"
model_colors <- c("RL-Sig65" = COL_RLSIG65, "Original R-loop-related score" = COL_RLOOP, "MKI67 expression" = COL_MKI67, "RL-Sig65 + MKI67" = COL_COMBINED)

score_long <- sample_scores %>% select(Sample, Cohort, Response, RL_Sig65_z, Original_Rloop_z, MKI67_z) %>% pivot_longer(cols = c(RL_Sig65_z, Original_Rloop_z, MKI67_z), names_to = "Score_key", values_to = "Score_z") %>% mutate(Score = recode(Score_key, RL_Sig65_z = "RL-Sig65", Original_Rloop_z = "Original R-loop-related score", MKI67_z = "MKI67 expression"), Score = factor(Score, levels = c("RL-Sig65", "Original R-loop-related score", "MKI67 expression")), Cohort = factor(Cohort, levels = cohort_order), Response = factor(Response, levels = c("RD", "pCR")))

main_wilcox <- score_long %>% group_by(Score) %>% summarise(P_value = wilcox.test(Score_z ~ Response)$p.value, .groups = "drop") %>% mutate(FDR = p.adjust(P_value, method = "BH"), star = fdr_to_star(FDR), x1 = 1, x2 = 2, y = c(3.4, 3.2, 3.8), Score = factor(Score, levels = levels(score_long$Score)))

pA <- ggplot(score_long, aes(x = Response, y = Score_z, fill = Response)) +
  geom_violin(trim = FALSE, alpha = 0.78, color = "#404040", linewidth = 0.55) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "#404040", linewidth = 0.45) +
  geom_jitter(width = 0.06, size = 0.6, alpha = 0.35, color = "#3A3A3A") +
  geom_segment(data = main_wilcox, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = main_wilcox, aes(x = x1, xend = x1, y = y - 0.08, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(data = main_wilcox, aes(x = x2, xend = x2, y = y - 0.08, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = main_wilcox, aes(x = 1.5, y = y + 0.12, label = star), inherit.aes = FALSE, size = 4.8, fontface = "bold") +
  facet_wrap(~ Score, nrow = 1) +
  scale_fill_manual(values = c("RD" = COL_RD, "pCR" = COL_PCR)) +
  scale_x_discrete(labels = c("RD\n(n=356)", "pCR\n(n=40)")) +
  coord_cartesian(ylim = c(-4.2, 4.6), clip = "off") +
  labs(x = NULL, y = "Within-cohort standardized score") +
  base_theme +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 13), axis.text = element_text(size = 11, color = "black"), panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.35), panel.grid.minor = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 5.5))

forest_levels <- rev(c("Random-effects pooled", cohort_order))
univ_plot_df <- association_all %>% mutate(Cohort = factor(Cohort, levels = forest_levels), Effect = factor(Effect, levels = c("RL-Sig65", "Original R-loop-related score", "MKI67 expression")))
pB <- ggplot(univ_plot_df, aes(x = OR, y = Cohort, xmin = CI_low, xmax = CI_high, color = Effect, shape = Effect)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, color = "#666666") +
  geom_errorbarh(height = 0.16, linewidth = 0.7, position = position_dodge(width = 0.56)) +
  geom_point(size = 2.9, position = position_dodge(width = 0.56)) +
  scale_x_log10(limits = c(0.35, 64), breaks = c(0.5, 1, 2, 4, 8, 16, 32), labels = c("0.5", "1", "2", "4", "8", "16", "32")) +
  scale_color_manual(values = model_colors[c("RL-Sig65", "Original R-loop-related score", "MKI67 expression")]) +
  labs(title = "Univariable associations", x = "OR for pCR per 1-SD increase", y = NULL, color = NULL, shape = NULL) +
  base_theme +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title.x = element_text(face = "bold", size = 12), axis.text = element_text(size = 11, color = "black"), legend.position = "none", panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35), plot.margin = margin(5.5, 5.5, 5.5, 5.5))

adj_plot_df <- adjusted_all %>% mutate(Cohort = factor(Cohort, levels = forest_levels))
pC <- ggplot(adj_plot_df, aes(x = OR, y = Cohort, xmin = CI_low, xmax = CI_high)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, color = "#666666") +
  geom_errorbarh(height = 0.16, linewidth = 0.78, color = COL_RLSIG65) +
  geom_point(size = 3.0, color = COL_RLSIG65) +
  scale_x_log10(limits = c(0.35, 64), breaks = c(0.5, 1, 2, 4, 8, 16, 32), labels = c("0.5", "1", "2", "4", "8", "16", "32")) +
  labs(title = "RL-Sig65 adjusted for MKI67", x = "Adjusted OR for pCR per 1-SD increase", y = NULL) +
  base_theme +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title.x = element_text(face = "bold", size = 12), axis.text = element_text(size = 11, color = "black"), panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35), plot.margin = margin(5.5, 5.5, 5.5, 5.5))

legend_bc <- cowplot::get_legend(
  pB + theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 10), legend.key.width = unit(1.4, "lines"), legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0))
)

loco_labels <- loco_roc_summary %>% mutate(short = c("RL-Sig65", "Original score", "MKI67", "RL-Sig65 + MKI67"), lab = paste0(short, " (AUC ", sprintf("%.3f", AUC), ")"))
comp_1 <- loco_delong %>% filter(Model_1 == "RL-Sig65", Model_2 == "Original R-loop-related score")
comp_2 <- loco_delong %>% filter(Model_1 == "RL-Sig65", Model_2 == "MKI67 expression")
roc_note <- paste0("RL-Sig65 vs original score: ΔAUC = ", sprintf("%.3f", comp_1$Delta_AUC_1_minus_2), ", BH-adjusted ", format_p(comp_1$DeLong_BH), "\n",
                   "RL-Sig65 vs MKI67: ΔAUC = ", sprintf("%.3f", comp_2$Delta_AUC_1_minus_2), ", BH-adjusted ", format_p(comp_2$DeLong_BH))

pD <- ggplot(loco_roc_curves, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.5, color = "#9A9A9A") +
  geom_path(linewidth = 1.0) +
  annotate("text", x = 0.50, y = 0.07, label = roc_note, hjust = 0.5, vjust = 0, size = 3.9, lineheight = 1.15) +
  scale_color_manual(values = model_colors) +
  coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
  labs(x = "1 - Specificity", y = "Sensitivity", color = NULL) +
  base_theme +
  theme(axis.title = element_text(face = "bold", size = 13), axis.text = element_text(size = 11, color = "black"), legend.position = "none", panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35), plot.margin = margin(5.5, 5.5, 5.5, 5.5))

legend_d <- cowplot::get_legend(
  pD + theme(legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 10), legend.key.width = unit(1.4, "lines"))
)

gsea_plot_df <- selected_pathways %>% mutate(Pathway = factor(Pathway, levels = rev(Pathway[order(NES)])), Dir2 = ifelse(Direction == "Enriched in pCR", "Enriched in pCR", "Enriched in RD"))
pE <- ggplot(gsea_plot_df, aes(x = NES, y = Pathway, fill = Dir2)) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "#666666") +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("Enriched in pCR" = COL_PCR, "Enriched in RD" = COL_RD)) +
  labs(title = "Meta-ranked Hallmark GSEA", x = "Normalized enrichment score", y = NULL, fill = NULL) +
  base_theme +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title.x = element_text(face = "bold", size = 12), axis.text = element_text(size = 10.5, color = "black"), legend.position = "none", plot.margin = margin(5.5, 5.5, 5.5, 5.5))

pathway_meta_plot_df <- pathway_meta %>% left_join(selected_pathways %>% select(Pathway_ID, NES), by = "Pathway_ID") %>% mutate(Pathway = factor(Pathway, levels = rev(levels(gsea_plot_df$Pathway))))
pF <- ggplot(pathway_meta_plot_df, aes(x = OR, y = Pathway, xmin = CI_low, xmax = CI_high)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, color = "#666666") +
  geom_errorbarh(height = 0.16, linewidth = 0.75, color = "#505050") +
  geom_point(aes(fill = NES), shape = 21, size = 3.2, color = "black", stroke = 0.4) +
  scale_x_log10(limits = c(0.35, 3.5), breaks = c(0.5, 1, 2, 3), labels = c("0.5", "1", "2", "3")) +
  scale_fill_gradient2(low = COL_PCR, mid = "white", high = COL_RD, midpoint = 0, name = "GSEA NES") +
  labs(title = "Pathway association with residual disease", x = "Adjusted OR for RD per 1-SD ssGSEA score", y = NULL) +
  base_theme +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5), axis.title.x = element_text(face = "bold", size = 12), axis.text = element_text(size = 10.5, color = "black"), legend.position = "none", panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35), plot.margin = margin(5.5, 5.5, 5.5, 5.5))

legend_ef1 <- cowplot::get_legend(pE + theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.key.width = unit(1.4, "lines")))
legend_ef2 <- cowplot::get_legend(pF + theme(legend.position = "bottom", legend.text = element_text(size = 10), legend.key.width = unit(1.2, "lines")))
legend_ef <- cowplot::plot_grid(legend_ef1, legend_ef2, nrow = 1, rel_widths = c(0.52, 0.48))

cohort_wilcox <- score_long %>% group_by(Cohort, Score) %>% summarise(P_value = wilcox.test(Score_z ~ Response)$p.value, .groups = "drop") %>% group_by(Score) %>% mutate(FDR = p.adjust(P_value, method = "BH"), star = fdr_to_star(FDR), label = paste0(star, "\n", format_fdr(FDR)), x = 1.5, y = 3.4) %>% ungroup()

pS4A <- ggplot(score_long, aes(x = Response, y = Score_z, fill = Response)) +
  geom_violin(trim = FALSE, alpha = 0.78, color = "#404040", linewidth = 0.4) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "#404040", linewidth = 0.35) +
  geom_jitter(width = 0.05, size = 0.4, alpha = 0.28, color = "#3A3A3A") +
  geom_text(data = cohort_wilcox, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.5, vjust = 0) +
  facet_grid(rows = vars(Score), cols = vars(Cohort)) +
  scale_fill_manual(values = c("RD" = COL_RD, "pCR" = COL_PCR)) +
  coord_cartesian(ylim = c(-4.1, 4.3), clip = "off") +
  labs(x = NULL, y = "Within-cohort standardized score") +
  base_theme +
  theme(legend.position = "none", strip.background = element_rect(fill = "#F2F2F2", color = "#CFCFCF"), strip.text = element_text(face = "bold", size = 9), axis.text = element_text(size = 8.4, color = "black"), axis.title = element_text(face = "bold", size = 10.5), panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.25))

auc_plot_df <- cohort_roc_summary %>% mutate(Cohort = factor(Cohort, levels = rev(cohort_order)), Model = factor(Model, levels = c("RL-Sig65", "Original R-loop-related score", "MKI67 expression")))
pS4B <- ggplot(auc_plot_df, aes(x = AUC, y = Cohort, xmin = AUC_CI_low, xmax = AUC_CI_high, color = Model, shape = Model)) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.45, color = "#666666") +
  geom_errorbarh(height = 0.16, linewidth = 0.72, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.9, position = position_dodge(width = 0.5)) +
  scale_x_continuous(limits = c(0.35, 1), breaks = seq(0.4, 1, 0.1)) +
  scale_color_manual(values = model_colors[c("RL-Sig65", "Original R-loop-related score", "MKI67 expression")]) +
  labs(x = "AUC for pCR discrimination (95% CI)", y = NULL, color = NULL, shape = NULL) +
  base_theme + theme(axis.title.x = element_text(face = "bold"), axis.text = element_text(color = "black"), legend.position = "bottom", panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35))

cohort_auc_labels <- cohort_roc_summary %>% mutate(short_model = recode(Model, "RL-Sig65" = "RL-Sig65", "Original R-loop-related score" = "Original", "MKI67 expression" = "MKI67"), line = paste0(short_model, " AUC ", sprintf("%.2f", AUC))) %>% group_by(Cohort) %>% summarise(label = paste(line, collapse = "\n"), .groups = "drop") %>% mutate(x = 0.97, y = 0.04)
pS4C <- ggplot(cohort_roc_curves, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, linewidth = 0.45, color = "#666666") +
  geom_path(linewidth = 0.9) +
  geom_text(data = cohort_auc_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1, vjust = 0, size = 2.5) +
  facet_wrap(~ Cohort, ncol = 3) +
  scale_color_manual(values = model_colors[c("RL-Sig65", "Original R-loop-related score", "MKI67 expression")]) +
  coord_equal() + labs(x = "1 - Specificity", y = "Sensitivity", color = NULL) + base_theme + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"), legend.position = "bottom", legend.text = element_text(size = 8))

pcr_baseline <- mean(sample_scores$Response_Binary)
pS4D <- ggplot(loco_pr_curves, aes(x = Recall, y = Precision, color = Model)) +
  geom_hline(yintercept = pcr_baseline, linetype = 2, linewidth = 0.45, color = "#666666") +
  geom_path(linewidth = 1) +
  scale_color_manual(values = model_colors) + coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) + labs(x = "Recall", y = "Precision", color = NULL) + base_theme + theme(axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"), legend.position = "bottom", panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35))

pS4E <- ggplot(cohort_gsea_all %>% mutate(Pathway = factor(Pathway, levels = rev(unique(meta_gsea$Pathway))), Cohort = factor(Cohort, levels = cohort_order), NegLog10FDR = -log10(pmax(FDR, 1e-300))), aes(x = Cohort, y = Pathway, color = NES, size = NegLog10FDR)) +
  geom_point(alpha = 0.9) + scale_color_gradient2(low = COL_PCR, mid = "white", high = COL_RD, midpoint = 0, name = "NES") + scale_size_continuous(range = c(0.7, 5.0), name = expression(-log[10](FDR))) + labs(x = NULL, y = NULL) + base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"), axis.text.y = element_text(size = 7.2, color = "black"), legend.position = "right")

qc_corr <- sample_scores %>% group_by(Cohort) %>% summarise(Spearman_rho = suppressWarnings(cor(RL_Sig65_z, MKI67_z, method = "spearman", use = "complete.obs")), P_value = suppressWarnings(cor.test(RL_Sig65_z, MKI67_z, method = "spearman")$p.value), .groups = "drop") %>% mutate(FDR = p.adjust(P_value, method = "BH"), label = paste0("ρ = ", sprintf("%.2f", Spearman_rho), "\n", format_fdr(FDR)), x = Inf, y = Inf)
pQC <- ggplot(sample_scores, aes(x = MKI67_z, y = RL_Sig65_z, color = Response)) +
  geom_point(size = 1.1, alpha = 0.65) + geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "black") + geom_text(data = qc_corr, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 3.2) + facet_wrap(~ Cohort, ncol = 3) + scale_color_manual(values = c("RD" = COL_RD, "pCR" = COL_PCR)) + labs(x = "Within-cohort standardized MKI67 expression", y = "Within-cohort standardized RL-Sig65 score", color = NULL) + base_theme + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"), legend.position = "bottom")

save_pdf("F3A.NAC_Score_Distribution.pdf", pA, width = 8.0, height = 4.8)
save_pdf("F3B.NAC_Univariable_OR_Forest.pdf", pB, width = 5.2, height = 4.6)
save_pdf("F3C.NAC_Adjusted_OR_Forest.pdf", pC, width = 5.0, height = 4.6)
save_pdf("F3D.NAC_LOCO_ROC.pdf", pD + theme(legend.position = "bottom"), width = 8.0, height = 4.8)
save_pdf("F3E.NAC_Meta_Hallmark_GSEA.pdf", pE, width = 4.0, height = 4.6)
save_pdf("F3F.NAC_Pathway_RD_Association.pdf", pF, width = 4.6, height = 4.6)

save_pdf("S4A.NAC_Cohort_Score_Distributions.pdf", pS4A, width = 13, height = 8)
save_pdf("S4B.NAC_Cohort_AUC_Forest.pdf", pS4B, width = 7.4, height = 5.4)
save_pdf("S4C.NAC_Cohort_ROC_Curves.pdf", pS4C, width = 10.8, height = ceiling(length(cohort_order) / 3) * 3.6 + 0.8)
save_pdf("S4D.NAC_LOCO_PR_Curves.pdf", pS4D, width = 6.2, height = 5.2)
save_pdf("S4E.NAC_Residual_Disease_Hallmark_by_Cohort.pdf", pS4E, width = 8.5, height = 12.5)
save_pdf("01.QC_RL-Sig65_MKI67_Correlation.pdf", pQC, width = 10, height = 6.8)

left_top <- pA + theme(plot.margin = margin(0, 4, 0, 0))
left_bottom <- cowplot::plot_grid(
  pB + theme(plot.margin = margin(0, 0, 0, 0)),
  pC + theme(plot.margin = margin(0, 0, 0, 0)),
  nrow = 1,
  rel_widths = c(1, 1)
)
left_body <- cowplot::plot_grid(left_top, left_bottom, ncol = 1, rel_heights = c(1, 1))

right_top <- pD + theme(plot.margin = margin(0, 0, 0, 4))
right_bottom <- cowplot::plot_grid(
  pE + theme(plot.margin = margin(0, 0, 0, 0)),
  pF + theme(plot.margin = margin(0, 0, 0, 0)),
  nrow = 1,
  rel_widths = c(1, 1)
)
right_body <- cowplot::plot_grid(right_top, right_bottom, ncol = 1, rel_heights = c(1, 1))

main_body <- cowplot::plot_grid(left_body, right_body, nrow = 1, rel_widths = c(1.62, 1.0), align = 'h', axis = 'tb')

left_legend_block <- legend_bc
right_legend_block <- cowplot::plot_grid(legend_d, legend_ef, ncol = 1, rel_heights = c(0.50, 0.50))
legend_row <- cowplot::plot_grid(left_legend_block, right_legend_block, nrow = 1, rel_widths = c(1.62, 1.0), align = 'h', axis = 'tb')

final_fig <- cowplot::plot_grid(main_body, legend_row, ncol = 1, rel_heights = c(1, 0.13))

final_fig <- cowplot::ggdraw(final_fig) +
  cowplot::draw_plot_label(
    label = c("A", "B", "C", "D", "E", "F"),
    x = c(0.025, 0.025, 0.315, 0.64, 0.64, 0.84),
    y = c(0.985, 0.50, 0.50, 0.985, 0.41, 0.41),
    size = c(20, 20, 20, 20, 20, 20),
    fontface = "bold",
    family = "Arial"
  )

save_pdf("F3.NAC_Main_Figure.pdf", final_fig, width = 15.5, height = 10.0)

xlsx_file <- file.path(OUTPUT_DIR, "T_NAC.NAC_Comprehensive_Results.xlsx")
wb <- openxlsx::createWorkbook()
header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAF7", halign = "center", valign = "center", wrapText = TRUE, border = "Bottom", borderColour = "#808080")
write_result_sheet <- function(wb, sheet_name, data, widths = "auto") {
  openxlsx::addWorksheet(wb, sheet_name, gridLines = FALSE)
  openxlsx::writeData(wb, sheet = sheet_name, x = data, startRow = 1, startCol = 1, colNames = TRUE, rowNames = FALSE, keepNA = FALSE, headerStyle = header_style, borders = "rows")
  openxlsx::addFilter(wb, sheet = sheet_name, row = 1, cols = seq_len(ncol(data)))
  openxlsx::freezePane(wb, sheet = sheet_name, firstActiveRow = 2, firstActiveCol = 2)
  openxlsx::setRowHeights(wb, sheet = sheet_name, rows = 1, heights = 38)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(data)), widths = widths)
}

write_result_sheet(wb, "Cohort_summary", cohort_summary)
write_result_sheet(wb, "Sample_scores", sample_scores %>% mutate(Cohort = as.character(Cohort), Response = as.character(Response)) %>% arrange(Cohort, Sample))
write_result_sheet(wb, "Association_meta", bind_rows(association_all, adjusted_all) %>% arrange(Effect, Cohort))
write_result_sheet(wb, "Cohort_ROC", cohort_roc_summary)
write_result_sheet(wb, "LOCO_performance", loco_roc_summary %>% left_join(loco_pr_summary %>% select(Model, AUPRC, pCR_prevalence), by = "Model"))
write_result_sheet(wb, "LOCO_DeLong", loco_delong)
write_result_sheet(wb, "Residual_DE_meta", residual_meta)
write_result_sheet(wb, "Hallmark_GSEA_meta", meta_gsea)
write_result_sheet(wb, "Hallmark_GSEA_cohort", cohort_gsea_all)
write_result_sheet(wb, "Pathway_RD_meta", bind_rows(pathway_association, pathway_meta))
write_result_sheet(wb, "LOCO_predictions", loco_predictions %>% mutate(Cohort = as.character(Cohort), Model = as.character(Model)) %>% arrange(Cohort, Sample, Model))
write_result_sheet(wb, "Score_Wilcoxon", main_wilcox %>% select(Score, P_value, FDR, star))
openxlsx::saveWorkbook(wb, file = xlsx_file, overwrite = TRUE)
