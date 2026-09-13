# GSE96058 preprocessing

PROJECT_DIR <- "."

raw_dir <- file.path(PROJECT_DIR, "data", "raw", "GSE96058")
out_dir <- file.path(PROJECT_DIR, "results", "bulk", "GSE96058")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(affy)
library(GEOquery)
library(tidyverse)
library(limma)
library(rtracklayer)
library(org.Hs.eg.db)
library(AnnotationDbi)

gset <- getGEO(filename = file.path(raw_dir, "GSE96058-GPL11154_series_matrix.txt.gz"),
               destdir = raw_dir, AnnotGPL = FALSE, getGPL = FALSE)

pdata <- pData(gset)

columns_needed <- c(
  "title",
  "geo_accession",
  "age at diagnosis:ch1",
  "er status:ch1",
  "pgr status:ch1",
  "her2 status:ch1",
  "ki67 status:ch1",
  "pam50 subtype:ch1",
  "chemo treated:ch1",
  "endocrine treated:ch1",
  "overall survival days:ch1",
  "overall survival event:ch1",
  "characteristics_ch1.10"
)

pdata_subset <- pdata[, columns_needed]

pdata_subset <- as.data.frame(lapply(pdata_subset, function(x) {
  x <- as.character(x)
  x[x == "NA"] <- NA
  return(x)
}))
pdata_subset <- na.omit(pdata_subset)

cols_to_numeric <- c("er.status.ch1", "pgr.status.ch1", "her2.status.ch1")
pdata_subset[, cols_to_numeric] <- lapply(pdata_subset[, cols_to_numeric],
                                          function(x) as.numeric(as.character(x)))

pdata_subset <- pdata_subset[
  ((pdata_subset$er.status.ch1 == 1) | (pdata_subset$pgr.status.ch1 == 1)) &
    (pdata_subset$her2.status.ch1 == 0),
]

pdata_subset$title <- as.character(pdata_subset$title)
pdata_subset <- pdata_subset[!grepl("repl", pdata_subset$title, ignore.case = TRUE), ]

pdata_subset <- pdata_subset[pdata_subset$pam50.subtype.ch1 %in% c("LumA", "LumB"), ]

cols_treat <- c("chemo.treated.ch1", "endocrine.treated.ch1")
pdata_subset[, cols_treat] <- lapply(pdata_subset[, cols_treat],
                                     function(x) as.numeric(as.character(x)))
pdata_subset <- pdata_subset[pdata_subset$endocrine.treated.ch1 != 0, ]

colnames(pdata_subset) <- c(
  "sample",
  "gsm_sample",
  "age",
  "ER",
  "PR",
  "HER2",
  "ki67",
  "PAM50",
  "chemo_tre",
  "endocrine_tre",
  "OS days",
  "OS event",
  "grade"
)

char_cols <- c("sample", "gsm_sample", "PAM50", "grade")

pdata_subset[] <- lapply(pdata_subset, as.character)

pdata_subset <- lapply(names(pdata_subset), function(col) {
  if (col %in% char_cols) {
    pdata_subset[[col]]
  } else {
    as.numeric(pdata_subset[[col]])
  }
})

pdata_subset <- as.data.frame(pdata_subset)
colnames(pdata_subset) <- c(
  "sample",
  "gsm_sample",
  "age",
  "ER",
  "PR",
  "HER2",
  "ki67",
  "PAM50",
  "chemo_tre",
  "endocrine_tre",
  "OS days",
  "OS event",
  "grade"
)

pdata_subset$grade <- sub(".*?:\\s*", "", pdata_subset$grade)

colnames(pdata_subset)[colnames(pdata_subset) == "OS days"] <- "OS years"

pdata_subset$`OS years` <- round(as.numeric(pdata_subset$`OS years`) / 365, 2)

gene_expr_data <- read.csv(file.path(raw_dir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv"),
                           header = TRUE, stringsAsFactors = FALSE, row.names = 1)

common_samples <- intersect(colnames(gene_expr_data), pdata_subset$sample)

gene_expr_data <- gene_expr_data[, common_samples]

gene_expr_data <- gene_expr_data[, pdata_subset$sample]

all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
current_ids <- rownames(gene_expr_data)

is_standard <- current_ids %in% all_symbols
n_standard <- sum(is_standard)
n_nonstandard <- sum(!is_standard)

unmapped_ids <- current_ids[!is_standard]
new_rownames <- current_ids

if(length(unmapped_ids) > 0) {
  alias_map <- tryCatch({
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unmapped_ids,
      keytype = "ALIAS",
      columns = "SYMBOL"
    )
  }, error = function(e) {
    data.frame(ALIAS = character(0), SYMBOL = character(0))
  })

  alias_map <- alias_map[!is.na(alias_map$SYMBOL), ]
  alias_map <- alias_map[!duplicated(alias_map$ALIAS), ]

  if(nrow(alias_map) > 0) {
    matched_idx <- match(alias_map$ALIAS, new_rownames)
    new_rownames[matched_idx] <- alias_map$SYMBOL
  }
}

is_standard_now <- new_rownames %in% all_symbols

unmapped_ids <- new_rownames[!is_standard_now]

if(length(unmapped_ids) > 0) {
  gtf_data <- import(file.path(raw_dir, "UCSC_hg38_knownGenes_22sep2014.gtf"))
  gtf_df <- as.data.frame(gtf_data)
  gtf_df_sub <- gtf_df[, c("mRNA", "geneSymbol")]

  gtf_df_sub <- gtf_df_sub[!is.na(gtf_df_sub$geneSymbol) & gtf_df_sub$geneSymbol != "", ]

  gtf_df_sub <- gtf_df_sub[!duplicated(gtf_df_sub$mRNA), ]

  gtf_map <- setNames(gtf_df_sub$geneSymbol, gtf_df_sub$mRNA)

  matched_idx <- new_rownames %in% names(gtf_map)
  if(sum(matched_idx) > 0) {
    new_rownames[matched_idx] <- gtf_map[new_rownames[matched_idx]]
  }
}

is_standard_now <- new_rownames %in% all_symbols

keep_idx <- !duplicated(new_rownames)
gene_expr_data <- gene_expr_data[keep_idx, ]
new_rownames <- new_rownames[keep_idx]

is_standard_final <- new_rownames %in% all_symbols

rownames(gene_expr_data) <- new_rownames
gene_expr_final <- gene_expr_data

if(sum(!is_standard_final) > 0) {
  print(head(new_rownames[!is_standard_final], 30))
}

gene_expr_final <- gene_expr_final[is_standard_final, ]

write.csv(gene_expr_final, file = file.path(out_dir, "GSE96058_exp.csv"), quote = FALSE)
write.csv(pdata_subset, file = file.path(out_dir, "GSE96058_clinisurv.csv"), row.names = FALSE, quote = FALSE)
