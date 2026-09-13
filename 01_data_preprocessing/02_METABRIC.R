# METABRIC preprocessing

PROJECT_DIR <- "."

rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)

data_dir <- file.path(PROJECT_DIR, "data", "raw", "METABRIC")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
setwd(data_dir)

cl <- data.table::fread("./brca_metabric/data_clinical_sample.txt", data.table = FALSE, skip = 4)

keep_cols <- c(
  "PATIENT_ID", "ER_STATUS", "HER2_STATUS", "PR_STATUS",
  "SAMPLE_TYPE", "TUMOR_SIZE", "TUMOR_STAGE", "TMB_NONSYNONYMOUS", "GRADE"
)
cl <- cl[, intersect(keep_cols, colnames(cl)), drop = FALSE]

cl <- cl %>%
  filter(ER_STATUS != "", HER2_STATUS != "", PR_STATUS != "")

cl$Type <- ifelse(
  (cl$ER_STATUS == "Positive" | cl$PR_STATUS == "Positive") & cl$HER2_STATUS == "Negative", "HR+/HER2-",
  ifelse(
    cl$HER2_STATUS == "Positive", "HER2+",
    ifelse(
      cl$ER_STATUS == "Negative" & cl$PR_STATUS == "Negative" & cl$HER2_STATUS == "Negative", "TNBC",
      "Other"
    )
  )
)

OS <- data.table::fread("./brca_metabric/data_clinical_patient.txt", data.table = FALSE, skip = 4)
colnames(OS)

keep_cols2 <- c(
  "PATIENT_ID", "AGE_AT_DIAGNOSIS", "SEX",
  "OS_MONTHS", "OS_STATUS", "RFS_MONTHS", "RFS_STATUS",
  "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY",
  "LYMPH_NODES_EXAMINED_POSITIVE", "NPI", "INFERRED_MENOPAUSAL_STATE"
)
OS <- OS[, intersect(keep_cols2, colnames(OS)), drop = FALSE]

OS <- OS[apply(OS, 1, function(x) all(!is.na(x) & x != "")), ]

pd <- cl %>% inner_join(OS, by = "PATIENT_ID")

old_names <- c(
  "PATIENT_ID", "ER_STATUS", "HER2_STATUS", "PR_STATUS", "SAMPLE_TYPE",
  "TUMOR_SIZE", "TUMOR_STAGE", "TMB_NONSYNONYMOUS", "GRADE", "Type",
  "AGE_AT_DIAGNOSIS", "SEX", "OS_MONTHS", "OS_STATUS", "RFS_MONTHS",
  "RFS_STATUS", "CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY",
  "LYMPH_NODES_EXAMINED_POSITIVE", "NPI", "INFERRED_MENOPAUSAL_STATE"
)
new_names <- c(
  "ID", "ER", "HER2", "PR", "SampleType",
  "TumorSize", "Stage", "TMB", "Grade", "Subtype",
  "Age", "Sex", "OS_mo", "OS_evt", "RFS_mo",
  "RFS_evt", "Chemo", "Hormone", "Radio",
  "LN_pos", "NPI", "Menopause"
)
colnames(pd) <- new_names

pd$OS_evt  <- as.integer(sub(":.*", "", as.character(pd$OS_evt)))
pd$RFS_evt <- as.integer(sub(":.*", "", as.character(pd$RFS_evt)))

exprSet <- data.table::fread(
  "./brca_metabric/data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt",
  data.table = FALSE
)
exprSet <- exprSet[, -2]

exprSet <- limma::avereps(exprSet[, -1], ID = exprSet$Hugo_Symbol)
exprSet  <- as.data.frame(exprSet)

df <- as.data.frame(t(exprSet))

rownames(pd) <- pd$ID
pd$ID <- NULL

common_samples <- intersect(rownames(pd), rownames(df))

pd <- pd[common_samples, ]
df <- df[common_samples, ]

df_final <- cbind(pd, df)

table(df_final$SampleType)
df_final$SampleType <- NULL

df_HRpos_HER2neg <- df_final[df_final$Subtype == "HR+/HER2-", ]
df_HER2pos       <- df_final[df_final$Subtype == "HER2+", ]
df_TNBC          <- df_final[df_final$Subtype == "TNBC", ]

exp_HRpos_HER2neg <- df_HRpos_HER2neg[, -(1:20)]
exp_HER2pos       <- df_HER2pos[, -(1:20)]
exp_TNBC          <- df_TNBC[, -(1:20)]

dir.create("HRpos_HER2neg", showWarnings = FALSE)
write.table(df_HRpos_HER2neg,  file = "HRpos_HER2neg/HRpos_HER2neg_clinisurv_exp.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(exp_HRpos_HER2neg, file = "HRpos_HER2neg/exp_HRpos_HER2neg.txt",           sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

dir.create("HER2pos", showWarnings = FALSE)
write.table(df_HER2pos,  file = "HER2pos/HER2pos_clinisurv_exp.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(exp_HER2pos, file = "HER2pos/exp_HER2pos.txt",           sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

dir.create("TNBC", showWarnings = FALSE)
write.table(df_TNBC,  file = "TNBC/TNBC_clinisurv_exp.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(exp_TNBC, file = "TNBC/exp_TNBC.txt",           sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
