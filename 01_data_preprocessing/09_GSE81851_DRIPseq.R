# GSE81851 DRIP-seq preprocessing

PROJECT_DIR <- "."

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ChIPseeker)
  library(GenomicRanges)
  library(IRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
})

raw_dir <- file.path(PROJECT_DIR, "data", "raw", "GSE81851")
out_dir <- file.path(PROJECT_DIR, "results", "GSE81851", "processed")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

raw_files <- list.files(raw_dir, pattern = "^GSE81851_DESeq_differential_calls\\.tsv$",
                        recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

if (length(raw_files) > 0) {
  file_info <- file.info(raw_files)
  raw_files <- raw_files[!is.na(file_info$isdir) & !file_info$isdir]
}

if (length(raw_files) == 0) {
  stop("GSE81851_DESeq_differential_calls.tsv, check.")
}

raw_file <- raw_files[1]

x <- fread(raw_file)

required_cols <- c(
  "chr", "start", "end",
  "log2FoldChange.T0_T2_DRIP", "padj.T0_T2_DRIP",
  "log2FoldChange.T0_T24_DRIP", "padj.T0_T24_DRIP",
  "log2FoldChange.T2_DRIP_over_Input", "padj.T2_DRIP_over_Input",
  "log2FoldChange.T24_DRIP_over_Input", "padj.T24_DRIP_over_Input"
)

missing_cols <- setdiff(required_cols, names(x))
if (length(missing_cols) > 0) stop("missing: ", paste(missing_cols, collapse = ", "))

fdr_cut <- 0.10

x[, Original_2h := !is.na(padj.T0_T2_DRIP) &
    padj.T0_T2_DRIP < fdr_cut &
    log2FoldChange.T0_T2_DRIP > 0]

x[, Original_24h := !is.na(padj.T0_T24_DRIP) &
    padj.T0_T24_DRIP < fdr_cut &
    log2FoldChange.T0_T24_DRIP > 0]

x[, Original_persistent := Original_2h & Original_24h]

x[, Strict_2h := Original_2h &
    !is.na(padj.T2_DRIP_over_Input) &
    padj.T2_DRIP_over_Input < fdr_cut &
    log2FoldChange.T2_DRIP_over_Input > 0]

x[, Strict_24h := Original_24h &
    !is.na(padj.T24_DRIP_over_Input) &
    padj.T24_DRIP_over_Input < fdr_cut &
    log2FoldChange.T24_DRIP_over_Input > 0]

x[, Strict_persistent := Strict_2h & Strict_24h]
x[, Strict_any := Strict_2h | Strict_24h]

region_table <- x[Strict_any == TRUE]

fwrite(region_table,
       file.path(out_dir, "01_E2DRIPsummary.csv"))

peak_gr <- GRanges(
  seqnames = as.character(region_table$chr),
  ranges = IRanges(start = as.numeric(region_table$start) + 1,
                   end = as.numeric(region_table$end))
)

anno <- annotatePeak(
  peak_gr,
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Hs.eg.db",
  verbose = FALSE
)

anno_df <- as.data.frame(anno)
anno_df$start_raw <- anno_df$start - 1

anno_key <- paste(anno_df$seqnames, anno_df$start_raw, anno_df$end, sep = ":")
region_key <- paste(region_table$chr, region_table$start, region_table$end, sep = ":")
idx <- match(anno_key, region_key)

add_cols <- c(
  "Original_2h", "Original_24h", "Original_persistent",
  "Strict_2h", "Strict_24h", "Strict_persistent", "Strict_any",
  "log2FoldChange.T0_T2_DRIP", "padj.T0_T2_DRIP",
  "log2FoldChange.T2_DRIP_over_Input", "padj.T2_DRIP_over_Input",
  "log2FoldChange.T0_T24_DRIP", "padj.T0_T24_DRIP",
  "log2FoldChange.T24_DRIP_over_Input", "padj.T24_DRIP_over_Input"
)

for (nm in add_cols) anno_df[[nm]] <- region_table[[nm]][idx]

anno_df$Locus_class <- "Distal intergenic"

anno_df$Locus_class[
  grepl("^Promoter", anno_df$annotation)
] <- "Promoter"

anno_df$Locus_class[
  grepl("Exon|Intron|UTR", anno_df$annotation, ignore.case = TRUE)
] <- "Gene body"

anno_df$Locus_class[
  grepl("^Downstream", anno_df$annotation)
] <- "Downstream"

anno_df$Direct_locus <- anno_df$Locus_class %in% c("Promoter", "Gene body")

fwrite(anno_df,
       file.path(out_dir, "02_E2DRIPgeneannotationsummary.csv"))

direct_df <- as.data.table(anno_df)[
  Direct_locus == TRUE & !is.na(SYMBOL) & SYMBOL != ""
]

gene_summary <- direct_df[, .(
  DRIP_2h = any(Strict_2h),
  DRIP_24h = any(Strict_24h),
  DRIP_persistent = any(Strict_persistent),
  N_direct_peaks = .N,
  N_promoter_peaks = sum(Locus_class == "Promoter"),
  N_genebody_peaks = sum(Locus_class == "Gene body")
), by = SYMBOL]

gene_summary[, DRIP_any := DRIP_2h | DRIP_24h]
setorder(gene_summary, SYMBOL)

fwrite(gene_summary,
       file.path(out_dir, "03_E2DRIPgenesummary.csv"))

candidate_genes <- c("HMGB2", "H2AZ2", "H2AFV", "UBALD2")

candidate_hits <- anno_df[
  !is.na(anno_df$SYMBOL) & anno_df$SYMBOL %in% candidate_genes, ,
  drop = FALSE
]

candidate_cols <- c(
  "SYMBOL", "seqnames", "start", "end", "annotation",
  "Locus_class", "Direct_locus", "distanceToTSS",
  "Strict_2h", "Strict_24h", "Strict_persistent",
  "log2FoldChange.T0_T2_DRIP", "padj.T0_T2_DRIP",
  "log2FoldChange.T2_DRIP_over_Input", "padj.T2_DRIP_over_Input",
  "log2FoldChange.T0_T24_DRIP", "padj.T0_T24_DRIP",
  "log2FoldChange.T24_DRIP_over_Input", "padj.T24_DRIP_over_Input"
)

candidate_hits <- candidate_hits[, candidate_cols, drop = FALSE]

fwrite(candidate_hits,
       file.path(out_dir, "04_candidategene_DRIPsummary.csv"))

if (nrow(candidate_hits) == 0) {
} else {
  print(candidate_hits[, c(
    "SYMBOL", "annotation", "Locus_class", "distanceToTSS",
    "Strict_2h", "Strict_24h", "Strict_persistent",
    "log2FoldChange.T0_T24_DRIP", "padj.T0_T24_DRIP"
  ), drop = FALSE])
}

basic_summary <- data.table(
  Section = "genenumber",
  Group = c(
    "raw",
    "raw_2h",
    "raw_24h",
    "raw_",
    "_2h",
    "_24h",
    "_",
    "_",
    "gene_2h",
    "gene_24h",
    "gene_",
    "gene_"
  ),
  N = c(
    nrow(x),
    sum(x$Original_2h),
    sum(x$Original_24h),
    sum(x$Original_persistent),
    sum(x$Strict_2h),
    sum(x$Strict_24h),
    sum(x$Strict_persistent),
    sum(x$Strict_any),
    sum(gene_summary$DRIP_2h),
    sum(gene_summary$DRIP_24h),
    sum(gene_summary$DRIP_persistent),
    sum(gene_summary$DRIP_any)
  ),
  Percent = NA_real_
)

make_locus_summary <- function(flag_name, label) {

  tmp <- as.data.table(anno_df)
  tmp <- tmp[get(flag_name) == TRUE, .N, by = Locus_class]
  tmp[, Percent := N / sum(N) * 100]
  tmp[, Section := paste0("gene_", label)]
  tmp[, Group := Locus_class]

  return(tmp[, .(Section, Group, N, Percent)])
}

locus_summary <- rbind(
  make_locus_summary("Strict_2h", "2h"),
  make_locus_summary("Strict_24h", "24h"),
  make_locus_summary("Strict_persistent", "")
)

final_summary <- rbind(basic_summary, locus_summary, fill = TRUE)

fwrite(final_summary,
       file.path(out_dir, "05_GSE81851_DRIPstatisticssummary.csv"))

enrichment_genes <- sort(unique(
  gene_summary$SYMBOL[gene_summary$DRIP_any == TRUE]
))

writeLines(enrichment_genes,
           file.path(out_dir, "06_enrichment_DRIPgene.txt"))

method_note <- c(
  "GSE81851 DRIP-seq analysis",
  "Genome build: hg19 / GRCh37.",
  "",
  "Original E2-induced criterion:",
  "E2 versus Mock log2FoldChange > 0 and adjusted P < 0.10.",
  "",
  "Strict E2-induced criterion:",
  "E2 versus Mock log2FoldChange > 0 and adjusted P < 0.10,",
  "and DRIP versus Input at the corresponding time point",
  "log2FoldChange > 0 and adjusted P < 0.10.",
  "",
  "Persistent DRIP region:",
  "A region satisfying the strict criteria at both 2 h and 24 h.",
  "",
  "Coordinate conversion:",
  "BED-style start coordinate was converted from 0-based to 1-based",
  "before GRanges annotation.",
  "",
  "Direct gene-locus DRIP evidence:",
  "Peak located in promoter within 1 kb of TSS, exon, intron or UTR.",
  "Distal intergenic peaks were retained in the annotation table",
  "but were not classified as direct DRIP evidence for the nearest gene."
)

writeLines(method_note,
           file.path(out_dir, "07_GSE81851_DRIPanalysis.txt"))
