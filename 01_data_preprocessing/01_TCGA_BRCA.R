# TCGA-BRCA preprocessing

PROJECT_DIR <- "."

data_dir <- file.path(PROJECT_DIR, "data", "raw", "TCGA_BRCA")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
setwd(data_dir)

suppressPackageStartupMessages({
  library(tidyverse)
  library(biomaRt)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dendextend)
  library(dendsort)
  library(ggrepel)
})

counts1<-read.table(file="TCGA-BRCA.star_counts.tsv",sep="\t",header=T)
rownames(counts1)<-counts1[ ,1]
counts1<-counts1[ ,-1]

table(substr(colnames(counts1),14,16))

counts1<-counts1[ ,substr(colnames(counts1),14,16)
                  %in%c("01A","11A")]

table(substr(colnames(counts1),14,16))

counts1 <- counts1[!duplicated(substr(rownames(counts1),1,15)),]
rownames(counts1) <- substr(rownames(counts1),1,15)

counts<-ceiling(2^(counts1)-1)

dir.create("count", showWarnings = FALSE)
write.table(counts, file = "count/BRCA_ENSgenename.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

library(org.Hs.eg.db)
library(AnnotationDbi)

counts <- read.table("count/BRCA_ENSgenename.txt",
                     sep = "\t",
                     row.names = 1,
                     check.names = FALSE,
                     stringsAsFactors = FALSE,
                     header = TRUE)

ensembl_ids <- rownames(counts)

annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  columns = c("SYMBOL", "GENETYPE"),
  keytype = "ENSEMBL")
colnames(annotations) <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
annotations <- annotations[!duplicated(annotations$ensembl_gene_id), ]

counts_df <- data.frame(ensembl_gene_id = rownames(counts), counts, stringsAsFactors = FALSE)
counts_annot <- merge(counts_df, annotations, by = "ensembl_gene_id", all.x = FALSE)
counts_annot <- counts_annot[!(is.na(counts_annot$hgnc_symbol) | counts_annot$hgnc_symbol == ""), ]

table(counts_annot$gene_biotype)

counts_annot <- counts_annot[, c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
                                 setdiff(colnames(counts_annot),
                                         c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")))]

convert_and_sum <- function(df, biotype_filter){
  df_sub <- df[df$gene_biotype %in% biotype_filter, ]
  summed <- df_sub %>%
    group_by(hgnc_symbol) %>%
    summarise(across(where(is.numeric), sum))
  summed_df <- as.data.frame(summed)
  rownames(summed_df) <- summed_df$hgnc_symbol
  summed_df <- summed_df[, -1]
  return(summed_df)
}

count_mRNA <- convert_and_sum(counts_annot, biotype_filter = c("protein-coding"))

write.table(count_mRNA, file = "count/BRCA_mRNA_symbolgenenames.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

extract_tumor_samples <- function(count_matrix) {
  tumor_samples <- colnames(count_matrix)[substr(colnames(count_matrix), 14, 16) == "01A"]
  count_matrix[, tumor_samples, drop = FALSE]
}
count_mRNA_01A <- extract_tumor_samples(count_mRNA)

write.table(count_mRNA_01A, file = "count/BRCA_01A_mRNA_symbolgenenames.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

fpkm1<-read.table(file="TCGA-BRCA.star_fpkm.tsv",sep="\t",header=T)
rownames(fpkm1)<-fpkm1[,1]
fpkm1<-fpkm1[,-1]

table(substr(colnames(fpkm1),14,16))
fpkm1<-fpkm1[ ,substr(colnames(fpkm1),14,16)
              %in%c("01A","11A")]
table(substr(colnames(fpkm1),14,16))

fpkm1 <- fpkm1[!duplicated(substr(rownames(fpkm1),1,15)),]
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm<-fpkm1

dir.create("fpkm", showWarnings = FALSE)
write.table(fpkm, file = "fpkm/BRCA_ENSgenename.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

fpkm <- read.table(
  file = "fpkm/BRCA_ENSgenename.txt",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE)

ensembl_ids <- rownames(fpkm)

annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  columns = c("SYMBOL", "GENETYPE"),
  keytype = "ENSEMBL")
colnames(annotations) <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
annotations <- annotations[!duplicated(annotations$ensembl_gene_id), ]

fpkm_df <- data.frame(ensembl_gene_id = rownames(fpkm), fpkm, stringsAsFactors = FALSE)
fpkm_annot <- merge(fpkm_df, annotations, by = "ensembl_gene_id", all.x = FALSE)
fpkm_annot <- fpkm_annot[!(is.na(fpkm_annot$hgnc_symbol) | fpkm_annot$hgnc_symbol == ""), ]

table(fpkm_annot$gene_biotype)

fpkm_annot <- fpkm_annot[, c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
                             setdiff(colnames(fpkm_annot),
                                     c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")))]

convert_and_sum <- function(df, biotype_filter){
  df_sub <- df[df$gene_biotype %in% biotype_filter, ]
  summed <- df_sub %>%
    group_by(hgnc_symbol) %>%
    summarise(across(where(is.numeric), sum))
  summed_df <- as.data.frame(summed)
  rownames(summed_df) <- summed_df$hgnc_symbol
  summed_df <- summed_df[, -1]
  return(summed_df)
}

fpkm_mRNA <- convert_and_sum(fpkm_annot, biotype_filter = c("protein-coding"))

write.table(fpkm_mRNA, file = "fpkm/BRCA_mRNA_symbolgenenames.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

extract_tumor_samples <- function(expr_mat) {
  tumor_samples <- colnames(expr_mat)[substr(colnames(expr_mat), 14, 16) == "01A"]
  expr_mat[, tumor_samples, drop = FALSE]}
fpkm_mRNA_01A <- extract_tumor_samples(fpkm_mRNA)

write.table(fpkm_mRNA_01A, file = "fpkm/BRCA_01A_mRNA_symbolgenenames.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

dir.create("tpm", showWarnings = FALSE)

fpkm_files <- list.files("fpkm", pattern = "\\.txt$", full.names = TRUE)

log2fpkm_to_tpm <- function(fpkm_log2) {
  fpkm <- 2^(fpkm_log2) - 1
  fpkm[fpkm == 0] <- 1e-6
  tpm <- apply(fpkm, 2, function(col) {
    exp(log(col) - log(sum(col)) + log(1e6))
  })
  as.data.frame(tpm)
}

for (file_path in fpkm_files) {

  file_name <- basename(file_path)

  fpkm_log2 <- read.table(file_path,
                          header = TRUE, sep = "\t",
                          row.names = 1, stringsAsFactors = FALSE)

  tpm <- log2fpkm_to_tpm(fpkm_log2)

  tpm_rounded <- round(tpm, 2)

  output_path <- file.path("tpm", paste0("tpm_", file_name))

  write.table(tpm_rounded, output_path,
              sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = NA)
}

pd <- read.delim("TCGA.BRCA.sampleMap_BRCA_clinicalMatrix")

pd2 <- pd[, c("ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012")]
rownames(pd2) <- pd$sampleID
rownames(pd2) <- paste0(gsub("-", ".", rownames(pd2)), "A")

tpm <- read.table(file = "tpm/tpm_BRCA_mRNA_symbolgenenames.txt", sep = "\t",
                  header = TRUE, row.names = 1, stringsAsFactors = FALSE)
counts <- read.table(file = "count/BRCA_mRNA_symbolgenenames.txt", sep = "\t",
                     header = TRUE, row.names = 1, stringsAsFactors = FALSE)

tpm_11A_samples <- colnames(tpm)[grepl("11A$", colnames(tpm))]
counts_11A_samples <- colnames(counts)[grepl("11A$", colnames(counts))]

subtype_filters <- list(
  HER2pos = function(x) x["HER2_Final_Status_nature2012"] == "Positive",
  HRpos_HER2neg = function(x) ((x["ER_Status_nature2012"] == "Positive" | x["PR_Status_nature2012"] == "Positive") & x["HER2_Final_Status_nature2012"] == "Negative"),
  TNBC = function(x) (x["ER_Status_nature2012"] == "Negative" & x["PR_Status_nature2012"] == "Negative" & x["HER2_Final_Status_nature2012"] == "Negative"))

sample_counts <- list()
for (subtype_name in names(subtype_filters)) {

  dir.create(file.path(subtype_name, "tpm"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(subtype_name, "count"), recursive = TRUE, showWarnings = FALSE)

  k <- apply(pd2, 1, subtype_filters[[subtype_name]])
  pd2_subtype <- pd2[k, ]
  pd2_subtype_01A <- rownames(pd2_subtype)[grepl("01A$", rownames(pd2_subtype))]
  tpm_subtype_01A <- intersect(pd2_subtype_01A, colnames(tpm))
  counts_subtype_01A <- intersect(pd2_subtype_01A, colnames(counts))

  tpm_subtype_all <- tpm[, c(tpm_subtype_01A, tpm_11A_samples), drop = FALSE]
  write.table(tpm_subtype_all, file = file.path(subtype_name, "tpm", paste0(subtype_name, ".txt")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  tpm_subtype_tumor <- tpm[, tpm_subtype_01A, drop = FALSE]
  write.table(tpm_subtype_tumor, file = file.path(subtype_name, "tpm", paste0(subtype_name, "_tumor_only.txt")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  counts_subtype_all <- counts[, c(counts_subtype_01A, counts_11A_samples), drop = FALSE]
  write.table(counts_subtype_all, file = file.path(subtype_name, "count", paste0(subtype_name, ".txt")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  counts_subtype_tumor <- counts[, counts_subtype_01A, drop = FALSE]
  write.table(counts_subtype_tumor, file = file.path(subtype_name, "count", paste0(subtype_name, "_tumor_only.txt")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  sample_counts[[subtype_name]] <- list(tpm = length(tpm_subtype_01A), counts = length(counts_subtype_01A))
}

library(tidyverse)

surv <- read.table(file = 'TCGA-BRCA.survival.tsv', sep = '\t', header = TRUE)
surv$sample <- gsub("-", ".", surv$sample)
rownames(surv) <- surv$sample
surv <- surv[, c("OS", "OS.time")]

clini <- read.table(file = 'TCGA-BRCA.clinical.tsv', sep = '\t', header = TRUE)
clini$sample <- gsub("-", ".", clini$sample)

colnames(clini)
table(clini$progression_or_recurrence.diagnoses)

table(clini$treatment_type.treatments.diagnoses)
table(clini$treatment_or_therapy.treatments.diagnoses)

keep_cols <- c(
  "sample",
  "age_at_diagnosis.diagnoses",
  "gender.demographic",
  "ajcc_pathologic_stage.diagnoses",
  "ajcc_pathologic_t.diagnoses",
  "ajcc_pathologic_n.diagnoses",
  "ajcc_pathologic_m.diagnoses",
  "sample_type.samples",
  "prior_malignancy.diagnoses",
  "synchronous_malignancy.diagnoses"
)

keep_cols_in_data <- intersect(keep_cols, colnames(clini))
clini_sub <- clini[, keep_cols_in_data, drop = FALSE]

colnames(clini_sub) <- c(
  "sample",
  "age",
  "gender",
  "stage",
  "T",
  "N",
  "M",
  "sample_type",
  "prior_malignancy",
  "synchronous_malignancy"
)

clini_sub <- clini_sub[clini_sub$gender == "female", ]
clini_sub <- clini_sub[clini_sub$sample_type == "Primary Tumor", ]
clini_sub <- clini_sub[!(clini_sub$prior_malignancy %in% c("yes", "not reported")), ]
clini_sub <- clini_sub[clini_sub$synchronous_malignancy == "No", ]
clini_sub <- clini_sub[, !(colnames(clini_sub) %in%
                             c("gender", "sample_type", "prior_malignancy", "synchronous_malignancy"))]

clini_sub$age <- round(as.numeric(clini_sub$age) / 365, 1)

clini_sub$stage <- gsub("Stage ", "", clini_sub$stage)
clini_sub$stage <- case_when(
  grepl("^IV",  clini_sub$stage) ~ 4,
  grepl("^III", clini_sub$stage) ~ 3,
  grepl("^II",  clini_sub$stage) ~ 2,
  grepl("^I",   clini_sub$stage) ~ 1,
  TRUE ~ NA_real_
)

clini_sub$T <- case_when(
  grepl("^Tis", clini_sub$T, ignore.case = TRUE) ~ 0,
  grepl("^T1",  clini_sub$T, ignore.case = TRUE) ~ 1,
  grepl("^T2",  clini_sub$T, ignore.case = TRUE) ~ 2,
  grepl("^T3",  clini_sub$T, ignore.case = TRUE) ~ 3,
  grepl("^T4",  clini_sub$T, ignore.case = TRUE) ~ 4,
  TRUE ~ NA_real_
)

clini_sub$N <- case_when(
  grepl("^N0", clini_sub$N, ignore.case = TRUE) ~ 0,
  grepl("^N1", clini_sub$N, ignore.case = TRUE) ~ 1,
  grepl("^N2", clini_sub$N, ignore.case = TRUE) ~ 2,
  grepl("^N3", clini_sub$N, ignore.case = TRUE) ~ 3,
  TRUE ~ NA_real_
)

clini_sub$M <- case_when(
  grepl("^M0\\(i\\+\\)", clini_sub$M, ignore.case = TRUE) ~ NA_real_,
  grepl("^M0", clini_sub$M, ignore.case = TRUE) ~ 0,
  grepl("^M1", clini_sub$M, ignore.case = TRUE) ~ 1,
  TRUE ~ NA_real_
)

clini_sub <- clini_sub[complete.cases(clini_sub), ]

table(clini_sub$stage, useNA = "ifany")
table(clini_sub$T,     useNA = "ifany")
table(clini_sub$N,     useNA = "ifany")
table(clini_sub$M,     useNA = "ifany")

common_samples <- intersect(clini_sub$sample, rownames(surv))
clini_common   <- clini_sub[clini_sub$sample %in% common_samples, ]
surv_common    <- surv[common_samples, ]
clini_surv     <- cbind(clini_common, surv_common[clini_common$sample, ])
rownames(clini_surv) <- clini_surv$sample

subtypes <- c("HER2pos", "HRpos_HER2neg", "TNBC")

for (subtype_name in subtypes) {

  tpm_tumor <- read.table(
    file = file.path(subtype_name, "tpm", paste0(subtype_name, "_tumor_only.txt")),
    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  tpm_tumor_t <- as.data.frame(t(tpm_tumor))
  common_tumor <- intersect(rownames(tpm_tumor_t), rownames(clini_surv))
  tpm_tumor_common    <- tpm_tumor_t[common_tumor, , drop = FALSE]
  clini_surv_common   <- clini_surv[common_tumor, ]
  merged_tumor        <- cbind(clini_surv_common, tpm_tumor_common)
  out_tumor <- file.path(subtype_name, paste0(subtype_name, "_clini_surv.txt"))
  write.table(merged_tumor, file = out_tumor,
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  tpm_all <- read.table(
    file = file.path(subtype_name, "tpm", paste0(subtype_name, ".txt")),
    sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
  tpm_all_t <- as.data.frame(t(tpm_all))

  all_samples  <- rownames(tpm_all_t)
  tumor_samples  <- all_samples[grepl("01A$", all_samples)]
  normal_samples <- all_samples[grepl("11A$", all_samples)]

  common_tumor_all <- intersect(tumor_samples, rownames(clini_surv))

  tpm_tumor_part  <- tpm_all_t[common_tumor_all, , drop = FALSE]
  clini_surv_part <- clini_surv[common_tumor_all, ]
  merged_tumor_part <- cbind(clini_surv_part, tpm_tumor_part)

  tpm_normal_part   <- tpm_all_t[normal_samples, , drop = FALSE]
  na_clini          <- as.data.frame(
    matrix(NA, nrow = length(normal_samples),
           ncol = ncol(clini_surv),
           dimnames = list(normal_samples, colnames(clini_surv))))
  merged_normal_part <- cbind(na_clini, tpm_normal_part)

  merged_all <- rbind(merged_tumor_part, merged_normal_part)
  out_all <- file.path(subtype_name, paste0(subtype_name, "_all_clini_surv.txt"))
  write.table(merged_all, file = out_all,
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}

subtypes <- c("HER2pos", "HRpos_HER2neg", "TNBC")
for (subtype_name in subtypes) {

  file_path <- file.path(subtype_name, paste0(subtype_name, "_all_clini_surv.txt"))
  data <- read.table(file = file_path,
                     sep = "\t",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)

  expr_data <- data[, -c(1:8)]

  expr_data_t <- t(expr_data)

  out_file <- file.path(subtype_name, paste0(subtype_name, "_all_clini_surv_exp.txt"))
  write.table(expr_data_t,
              file = out_file,
              sep = "\t",
              quote = FALSE,
              row.names = TRUE,
              col.names = NA)

  all_samples <- colnames(expr_data_t)
  tumor_samples <- sum(grepl("01A$", all_samples))
  normal_samples <- sum(grepl("11A$", all_samples))

}
