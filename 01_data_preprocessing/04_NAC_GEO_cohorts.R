# NAC GEO cohort preprocessing

PROJECT_DIR <- "."

NAC_RAW_DIR <- file.path(PROJECT_DIR, "data/raw/NAC")
dir.create(NAC_RAW_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(NAC_RAW_DIR)
NAC_OUT_DIR <- file.path(PROJECT_DIR, "data/processed/NAC")
dir.create(NAC_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

library(affy)
library(GEOquery)
library(tidyverse)
library(limma)

gpl <- getGEO("GPL96", destdir = ".", AnnotGPL = FALSE)

cel_files <- list.celfiles(path = "./GSE20271_RAW", full.names = TRUE)
raw_data <- ReadAffy(filenames = cel_files)

norm_data <- rma(raw_data)

exp <- exprs(norm_data)

gset <- getGEO("GSE20271", destdir = ".", AnnotGPL=F, getGPL=F)
pdata <- pData(gset[[1]])

pdata<-pdata[,c(2,10,22,24,25)]

extract_after_colon <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    sub("^[^:]*:(.*)$", "\\1", x)
  )
}

pdata_cleaned <- pdata %>%
  mutate(across(-geo_accession, extract_after_colon))

pdata <- pdata_cleaned

colnames(pdata) <- c("sample", "group","PR","ER","HER2")

pdata <- pdata %>%

  mutate(
    PR = toupper(trimws(PR)),
    ER = toupper(trimws(ER)),
    HER2 = toupper(trimws(HER2)),
    group = trimws(group)
  ) %>%

  mutate(
    Type = case_when(
      (ER == "P" | PR == "P") & HER2 == "N" ~ "HR+/HER2-",
      HER2 == "P" ~ "HER2+",
      ER == "N" & PR == "N" & HER2 == "N" ~ "TNBC",
      TRUE ~ "Other"
    )
  ) %>%

  mutate(
    subgroup = case_when(
      group == "pCR" & Type == "HR+/HER2-" ~ "pCR_HR+/HER2-",
      group == "RD" & Type == "HR+/HER2-" ~ "RD_HR+/HER2-",
      TRUE ~ "Other"
    )
  )

library(stringr)

colnames(exp) <- toupper(colnames(exp)) %>%
  str_replace("\\.CEL$", "") %>%
  str_extract("GSM\\d+")

pdata$sample_upper <- toupper(pdata$sample)

filtered_pdata <- pdata %>%
  filter(
    (group == "pCR" & Type == "HR+/HER2-") |
      (group == "RD" & Type == "HR+/HER2-")
  )

samples_to_keep <- filtered_pdata$sample_upper

exp_filtered <- exp[, colnames(exp) %in% samples_to_keep]

pdata_filtered <- filtered_pdata

sample_to_group <- setNames(pdata_filtered$group, pdata_filtered$sample_upper)
new_colnames <- paste0(colnames(exp_filtered), "_", sample_to_group[colnames(exp_filtered)])
colnames(exp_filtered) <- new_colnames

ids <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(ids) <- c("probe_id", "symbol")
exp_df <- as.data.frame(exp_filtered)
exp_df$probe_id <- rownames(exp_df)
exp_annot <- inner_join(ids, exp_df, by = "probe_id")

exp_annot <- exp_annot[!duplicated(exp_annot$symbol), ]
rownames(exp_annot) <- exp_annot$symbol

exp_annot <- exp_annot[, !(colnames(exp_annot) %in% c("probe_id", "symbol"))]

library(tidyverse)
exp_annot_df <- as.data.frame(exp_annot)

exp_annot_df <- exp_annot_df %>% rownames_to_column(var = "Gene")

exp_annot_df <- exp_annot_df %>%
  separate_rows(Gene, sep = " /// ") %>%
  mutate(Gene = str_trim(Gene)) %>%
  distinct(Gene, .keep_all = TRUE)

exp_annot_df <- exp_annot_df %>% column_to_rownames(var = "Gene")
exp_annot <- exp_annot_df

table(sub("^.*?_(.*)$", "\\1", colnames(exp_annot)))

write.table(exp_annot, file = file.path(NAC_OUT_DIR, "GSE20271_RMA_normalized_exp.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

library(affy)
library(GEOquery)
library(tidyverse)
library(limma)

gpl <- getGEO("GPL96", destdir = ".", AnnotGPL = FALSE)

cel_files <- list.celfiles(path = "./GSE22093_RAW", full.names = TRUE)
raw_data <- ReadAffy(filenames = cel_files)

norm_data <- rma(raw_data)

exp <- exprs(norm_data)

gset <- getGEO("GSE22093", destdir = ".", AnnotGPL=F, getGPL=F)
pdata <- pData(gset[[1]])

pdata<-pdata[,c(2,11,54)]

extract_after_colon <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    sub("^[^:]*:(.*)$", "\\1", x)
  )
}

pdata_cleaned <- pdata %>%
  mutate(across(-geo_accession, extract_after_colon))

pdata <- pdata_cleaned

colnames(pdata) <- c("sample", "group","ER")

pdata$ER <- ifelse(grepl("neg", pdata$ER, ignore.case = TRUE), "N",
                   ifelse(grepl("pos", pdata$ER, ignore.case = TRUE), "P", "other"))

head(pdata)

pdata <- pdata %>%

  mutate(
    ER = toupper(trimws(ER)),
    group = trimws(group)
  ) %>%

  mutate(
    Type = case_when(
      ER == "P"  ~ "HR+/HER2-",
      TRUE ~ "Other"
    )
  ) %>%

  mutate(
    subgroup = case_when(
      group == "pCR" & Type == "HR+/HER2-" ~ "pCR_HR+/HER2-",
      group == "RD" & Type == "HR+/HER2-" ~ "RD_HR+/HER2-",
      TRUE ~ "Other"
    )
  )

library(stringr)

colnames(exp) <- toupper(colnames(exp)) %>%
  str_replace("\\.CEL$", "") %>%
  str_extract("GSM\\d+")

pdata$sample_upper <- toupper(pdata$sample)

filtered_pdata <- pdata %>%
  filter(
    (group == "pCR" & Type == "HR+/HER2-") |
      (group == "RD" & Type == "HR+/HER2-")
  )

samples_to_keep <- filtered_pdata$sample_upper

exp_filtered <- exp[, colnames(exp) %in% samples_to_keep]

pdata_filtered <- filtered_pdata

sample_to_group <- setNames(pdata_filtered$group, pdata_filtered$sample_upper)
new_colnames <- paste0(colnames(exp_filtered), "_", sample_to_group[colnames(exp_filtered)])
colnames(exp_filtered) <- new_colnames

ids <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(ids) <- c("probe_id", "symbol")
exp_df <- as.data.frame(exp_filtered)
exp_df$probe_id <- rownames(exp_df)
exp_annot <- inner_join(ids, exp_df, by = "probe_id")

exp_annot <- exp_annot[!duplicated(exp_annot$symbol), ]
rownames(exp_annot) <- exp_annot$symbol

exp_annot <- exp_annot[, !(colnames(exp_annot) %in% c("probe_id", "symbol"))]

library(tidyverse)
exp_annot_df <- as.data.frame(exp_annot)

exp_annot_df <- exp_annot_df %>% rownames_to_column(var = "Gene")

exp_annot_df <- exp_annot_df %>%
  separate_rows(Gene, sep = " /// ") %>%
  mutate(Gene = str_trim(Gene)) %>%
  distinct(Gene, .keep_all = TRUE)

exp_annot_df <- exp_annot_df %>% column_to_rownames(var = "Gene")
exp_annot <- exp_annot_df

table(sub("^.*?_(.*)$", "\\1", colnames(exp_annot)))

write.table(exp_annot, file = file.path(NAC_OUT_DIR, "GSE22093_RMA_normalized_exp.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

library(affy)
library(GEOquery)
library(tidyverse)
library(limma)

gpl <- getGEO("GPL96", destdir = ".", AnnotGPL = FALSE)

cel_files <- list.celfiles(path = "./GSE23988_RAW", full.names = TRUE)
raw_data <- ReadAffy(filenames = cel_files)

norm_data <- rma(raw_data)

exp <- exprs(norm_data)

gset <- getGEO("GSE23988", destdir = ".", AnnotGPL=F, getGPL=F)
pdata <- pData(gset[[1]])

pdata<-pdata[,c(2,10,43)]

extract_after_colon <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    sub("^[^:]*:(.*)$", "\\1", x)
  )
}

pdata_cleaned <- pdata %>%
  mutate(across(-geo_accession, extract_after_colon))

pdata <- pdata_cleaned

colnames(pdata) <- c("sample", "group","ER")

pdata$ER <- ifelse(grepl("neg", pdata$ER, ignore.case = TRUE), "N",
                   ifelse(grepl("pos", pdata$ER, ignore.case = TRUE), "P", "other"))

head(pdata)

pdata <- pdata %>%

  mutate(
    ER = toupper(trimws(ER)),
    group = trimws(group)
  ) %>%

  mutate(
    Type = case_when(
      ER == "P"  ~ "HR+/HER2-",
      TRUE ~ "Other"
    )
  ) %>%

  mutate(
    subgroup = case_when(
      group == "pCR" & Type == "HR+/HER2-" ~ "pCR_HR+/HER2-",
      group == "RD" & Type == "HR+/HER2-" ~ "RD_HR+/HER2-",
      TRUE ~ "Other"
    )
  )

library(stringr)

colnames(exp) <- toupper(colnames(exp)) %>%
  str_replace("\\.CEL$", "") %>%
  str_extract("GSM\\d+")

pdata$sample_upper <- toupper(pdata$sample)

filtered_pdata <- pdata %>%
  filter(
    (group == "pCR" & Type == "HR+/HER2-") |
      (group == "RD" & Type == "HR+/HER2-")
  )

samples_to_keep <- filtered_pdata$sample_upper

exp_filtered <- exp[, colnames(exp) %in% samples_to_keep]

pdata_filtered <- filtered_pdata

sample_to_group <- setNames(pdata_filtered$group, pdata_filtered$sample_upper)
new_colnames <- paste0(colnames(exp_filtered), "_", sample_to_group[colnames(exp_filtered)])
colnames(exp_filtered) <- new_colnames

ids <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(ids) <- c("probe_id", "symbol")
exp_df <- as.data.frame(exp_filtered)
exp_df$probe_id <- rownames(exp_df)
exp_annot <- inner_join(ids, exp_df, by = "probe_id")

exp_annot <- exp_annot[!duplicated(exp_annot$symbol), ]
rownames(exp_annot) <- exp_annot$symbol

exp_annot <- exp_annot[, !(colnames(exp_annot) %in% c("probe_id", "symbol"))]

library(tidyverse)
exp_annot_df <- as.data.frame(exp_annot)

exp_annot_df <- exp_annot_df %>% rownames_to_column(var = "Gene")

exp_annot_df <- exp_annot_df %>%
  separate_rows(Gene, sep = " /// ") %>%
  mutate(Gene = str_trim(Gene)) %>%
  distinct(Gene, .keep_all = TRUE)

exp_annot_df <- exp_annot_df %>% column_to_rownames(var = "Gene")
exp_annot <- exp_annot_df

table(sub("^.*?_(.*)$", "\\1", colnames(exp_annot)))

write.table(exp_annot, file = file.path(NAC_OUT_DIR, "GSE23988_RMA_normalized_exp.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

setwd(file.path(PROJECT_DIR, "data/raw/NAC"))

library(affy)
library(GEOquery)
library(tidyverse)
library(limma)

gpl <- getGEO("GPL570", destdir = ".", AnnotGPL = FALSE)

cel_files <- list.celfiles(path = "./GSE32646_RAW", full.names = TRUE)
raw_data <- ReadAffy(filenames = cel_files)

norm_data <- rma(raw_data)

exp <- exprs(norm_data)

gset <- getGEO("GSE32646", destdir = ".", AnnotGPL=F, getGPL=F)
pdata <- pData(gset[[1]])

pdata<-pdata[,c(2,47,48,51,52)]

extract_after_colon <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    sub("^[^:]*:(.*)$", "\\1", x)
  )
}

pdata_cleaned <- pdata %>%
  mutate(across(-geo_accession, extract_after_colon))

pdata <- pdata_cleaned

colnames(pdata) <- c("sample","ER","HER2","group","PR")

pdata$ER <- ifelse(grepl("neg", pdata$ER, ignore.case = TRUE), "N",
                   ifelse(grepl("pos", pdata$ER, ignore.case = TRUE), "P", "other"))
pdata$HER2 <- ifelse(grepl("neg", pdata$HER2, ignore.case = TRUE), "N",
                     ifelse(grepl("pos", pdata$HER2, ignore.case = TRUE), "P", "other"))
pdata$PR <- ifelse(grepl("neg", pdata$PR, ignore.case = TRUE), "N",
                   ifelse(grepl("pos", pdata$PR, ignore.case = TRUE), "P", "other"))
pdata$group <- ifelse(grepl("n", pdata$group, ignore.case = TRUE), "RD",
                      ifelse(grepl("p", pdata$group, ignore.case = TRUE), "pCR", "other"))

head(pdata)

pdata <- pdata %>%

  mutate(
    PR = toupper(trimws(PR)),
    ER = toupper(trimws(ER)),
    HER2 = toupper(trimws(HER2)),
    group = trimws(group)
  ) %>%

  mutate(
    Type = case_when(
      (ER == "P" | PR == "P") & HER2 == "N" ~ "HR+/HER2-",
      HER2 == "P" ~ "HER2+",
      ER == "N" & PR == "N" & HER2 == "N" ~ "TNBC",
      TRUE ~ "Other"
    )
  ) %>%

  mutate(
    subgroup = case_when(
      group == "pCR" & Type == "HR+/HER2-" ~ "pCR_HR+/HER2-",
      group == "RD" & Type == "HR+/HER2-" ~ "RD_HR+/HER2-",
      TRUE ~ "Other"
    )
  )

library(stringr)

colnames(exp) <- toupper(colnames(exp)) %>%
  str_replace("\\.CEL$", "") %>%
  str_extract("GSM\\d+")

pdata$sample_upper <- toupper(pdata$sample)

filtered_pdata <- pdata %>%
  filter(
    (group == "pCR" & Type == "HR+/HER2-") |
      (group == "RD" & Type == "HR+/HER2-")
  )

samples_to_keep <- filtered_pdata$sample_upper

exp_filtered <- exp[, colnames(exp) %in% samples_to_keep]

pdata_filtered <- filtered_pdata

sample_to_group <- setNames(pdata_filtered$group, pdata_filtered$sample_upper)
new_colnames <- paste0(colnames(exp_filtered), "_", sample_to_group[colnames(exp_filtered)])
colnames(exp_filtered) <- new_colnames

ids <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(ids) <- c("probe_id", "symbol")
exp_df <- as.data.frame(exp_filtered)
exp_df$probe_id <- rownames(exp_df)
exp_annot <- inner_join(ids, exp_df, by = "probe_id")

exp_annot <- exp_annot[!duplicated(exp_annot$symbol), ]
rownames(exp_annot) <- exp_annot$symbol

exp_annot <- exp_annot[, !(colnames(exp_annot) %in% c("probe_id", "symbol"))]

library(tidyverse)
exp_annot_df <- as.data.frame(exp_annot)

exp_annot_df <- exp_annot_df %>% rownames_to_column(var = "Gene")

exp_annot_df <- exp_annot_df %>%
  separate_rows(Gene, sep = " /// ") %>%
  mutate(Gene = str_trim(Gene)) %>%
  distinct(Gene, .keep_all = TRUE)

exp_annot_df <- exp_annot_df %>% column_to_rownames(var = "Gene")
exp_annot <- exp_annot_df

table(sub("^.*?_(.*)$", "\\1", colnames(exp_annot)))

write.table(exp_annot, file = file.path(NAC_OUT_DIR, "GSE32646_RMA_normalized_exp.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

setwd(file.path(PROJECT_DIR, "data/raw/NAC"))

library(affy)
library(GEOquery)
library(tidyverse)
library(limma)

gpl <- getGEO("GPL96", destdir = ".", AnnotGPL = FALSE)

cel_files <- list.celfiles(path = "./GSE25066_RAW", full.names = TRUE)
raw_data <- ReadAffy(filenames = cel_files)

norm_data <- rma(raw_data)

exp <- exprs(norm_data)

gset <- getGEO("GSE25066", destdir = ".", AnnotGPL=F, getGPL=F)
pdata <- pData(gset[[1]])

pdata<-pdata[,c(2,13,14,15,21)]

extract_after_colon <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    sub("^[^:]*:(.*)$", "\\1", x)
  )
}

pdata_cleaned <- pdata %>%
  mutate(across(-geo_accession, extract_after_colon))

pdata <- pdata_cleaned

colnames(pdata) <- c("sample","ER","PR","HER2","group")

pdata <- pdata %>%

  mutate(
    PR = toupper(trimws(PR)),
    ER = toupper(trimws(ER)),
    HER2 = toupper(trimws(HER2)),
    group = trimws(group)
  ) %>%

  mutate(
    Type = case_when(
      (ER == "P" | PR == "P") & HER2 == "N" ~ "HR+/HER2-",
      HER2 == "P" ~ "HER2+",
      ER == "N" & PR == "N" & HER2 == "N" ~ "TNBC",
      TRUE ~ "Other"
    )
  ) %>%

  mutate(
    subgroup = case_when(
      group == "pCR" & Type == "HR+/HER2-" ~ "pCR_HR+/HER2-",
      group == "RD" & Type == "HR+/HER2-" ~ "RD_HR+/HER2-",
      TRUE ~ "Other"
    )
  )

library(stringr)

colnames(exp) <- toupper(colnames(exp)) %>%
  str_replace("\\.CEL$", "") %>%
  str_extract("GSM\\d+")

pdata$sample_upper <- toupper(pdata$sample)

filtered_pdata <- pdata %>%
  filter(
    (group == "pCR" & Type == "HR+/HER2-") |
      (group == "RD" & Type == "HR+/HER2-")
  )

samples_to_keep <- filtered_pdata$sample_upper

exp_filtered <- exp[, colnames(exp) %in% samples_to_keep]

pdata_filtered <- filtered_pdata

sample_to_group <- setNames(pdata_filtered$group, pdata_filtered$sample_upper)
new_colnames <- paste0(colnames(exp_filtered), "_", sample_to_group[colnames(exp_filtered)])
colnames(exp_filtered) <- new_colnames

ids <- Table(gpl)[, c("ID", "Gene Symbol")]
colnames(ids) <- c("probe_id", "symbol")
exp_df <- as.data.frame(exp_filtered)
exp_df$probe_id <- rownames(exp_df)
exp_annot <- inner_join(ids, exp_df, by = "probe_id")

exp_annot <- exp_annot[!duplicated(exp_annot$symbol), ]
rownames(exp_annot) <- exp_annot$symbol

exp_annot <- exp_annot[, !(colnames(exp_annot) %in% c("probe_id", "symbol"))]

library(tidyverse)
exp_annot_df <- as.data.frame(exp_annot)

exp_annot_df <- exp_annot_df %>% rownames_to_column(var = "Gene")

exp_annot_df <- exp_annot_df %>%
  separate_rows(Gene, sep = " /// ") %>%
  mutate(Gene = str_trim(Gene)) %>%
  distinct(Gene, .keep_all = TRUE)

exp_annot_df <- exp_annot_df %>% column_to_rownames(var = "Gene")
exp_annot <- exp_annot_df

table(sub("^.*?_(.*)$", "\\1", colnames(exp_annot)))

write.table(exp_annot, file = file.path(NAC_OUT_DIR, "GSE25066_RMA_normalized_exp.txt"),
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
