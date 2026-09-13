# Cross-network perturbation and Hallmark analysis

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(ggplot2)
})

# Paths
project_dir <- "."
out_dir <- file.path(project_dir, "results", "scTenifoldKnk", "cross_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

c4_file <- file.path(
  project_dir, "results", "scTenifoldKnk", "C4",
  "03_three_gene_differential_regulation.csv"
)
c5_file <- file.path(
  project_dir, "results", "scTenifoldKnk", "C5",
  "03_three_gene_differential_regulation.csv"
)
rl65_file <- file.path(project_dir, "gene_sets", "RL_Sig65.txt")
c5sig_file <- file.path(project_dir, "gene_sets", "C5_stable_up_genes.txt")
hallmark_gmt <- file.path(
  project_dir, "data", "gene_sets",
  "h.all.v2026.1.Hs.symbols.gmt"
)

required <- c(c4_file, c5_file, rl65_file, c5sig_file, hallmark_gmt)
if (any(!file.exists(required))) {
  stop("Missing required input files:\n",
       paste(required[!file.exists(required)], collapse = "\n"))
}

# Utilities
canon <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x == "H2AFV"] <- "H2AZ2"
  x
}

read_set <- function(file) {
  x <- canon(readLines(file, warn = FALSE))
  unique(x[!is.na(x) & nzchar(x)])
}

read_ko <- function(file, network) {
  d <- fread(file)
  if (!all(c("KO", "gene", "Z") %in% names(d))) {
    stop("KO, gene and Z columns are required in: ", file)
  }
  d[, gene := canon(gene)]
  d[, KO := canon(KO)]
  d[, Network := network]
  if (!"Self_KO" %in% names(d)) d[, Self_KO := gene == KO]
  d[Self_KO == FALSE & is.finite(Z)]
}

rl65 <- read_set(rl65_file)
c5sig <- read_set(c5sig_file)
gene_sets <- list(
  `RL-Sig65` = rl65,
  `C5 stable-up` = c5sig
)

all_dr <- rbindlist(list(
  read_ko(c5_file, "GSE306201 C5"),
  read_ko(c4_file, "GSE245601 C4")
))

# Program-level perturbation
program_test <- function(d, genes) {
  z_set <- d[gene %in% genes, Z]
  z_bg <- d[!gene %in% genes, Z]

  if (length(z_set) < 5 || length(z_bg) < 5) {
    return(data.table(
      N_in_network = length(z_set),
      Median_Z = NA_real_,
      Background_median_Z = NA_real_,
      Delta_median_Z = NA_real_,
      P = NA_real_
    ))
  }

  wt <- wilcox.test(
    z_set,
    z_bg,
    alternative = "greater",
    exact = FALSE
  )

  data.table(
    N_in_network = length(z_set),
    Median_Z = median(z_set, na.rm = TRUE),
    Background_median_Z = median(z_bg, na.rm = TRUE),
    Delta_median_Z = median(z_set, na.rm = TRUE) -
      median(z_bg, na.rm = TRUE),
    P = wt$p.value
  )
}

program_stats <- all_dr[, {
  rbindlist(lapply(names(gene_sets), function(nm) {
    z <- program_test(.SD, gene_sets[[nm]])
    z[, Program := nm]
    z
  }))
}, by = .(Network, KO)]

program_stats[, FDR := p.adjust(P, method = "BH")]
fwrite(
  program_stats,
  file.path(out_dir, "01_program_perturbation_statistics.csv")
)

# Hallmark GSEA
gmt <- strsplit(readLines(hallmark_gmt, warn = FALSE), "\t")
term2gene <- rbindlist(lapply(gmt, function(x) {
  data.table(
    term = x[1],
    gene = canon(x[-c(1, 2)])
  )
}))
term2gene <- unique(term2gene[nzchar(gene)])

run_gsea <- function(d) {
  ranks <- d[, .(Z = max(Z, na.rm = TRUE)), by = gene]
  ranks <- ranks[is.finite(Z)]
  setorder(ranks, -Z)

  gene_list <- ranks$Z
  names(gene_list) <- ranks$gene

  fit <- suppressWarnings(
    clusterProfiler::GSEA(
      gene_list,
      TERM2GENE = term2gene,
      pvalueCutoff = 1,
      verbose = FALSE,
      seed = TRUE
    )
  )

  as.data.table(fit@result)
}

hallmark <- all_dr[, {
  z <- run_gsea(.SD)
  if (!nrow(z)) data.table() else z
}, by = .(Network, KO)]

if (nrow(hallmark)) {
  hallmark[, p.adjust := as.numeric(p.adjust)]
  fwrite(
    hallmark,
    file.path(out_dir, "02_Hallmark_GSEA_all.csv")
  )

  concordant <- hallmark[
    p.adjust < 0.05 & NES > 0,
    .(
      Networks = uniqueN(Network),
      Min_NES = min(NES),
      Max_FDR = max(p.adjust)
    ),
    by = .(KO, ID, Description)
  ][Networks == 2]

  fwrite(
    concordant,
    file.path(out_dir, "03_Hallmark_positive_concordance.csv")
  )
}

# Figure
p <- ggplot(
  program_stats,
  aes(x = KO, y = Delta_median_Z, shape = Program)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_point(size = 2) +
  facet_wrap(~Network, nrow = 1) +
  theme_classic(base_size = 9) +
  labs(
    x = NULL,
    y = "Delta median regulatory perturbation Z"
  )

ggsave(
  file.path(out_dir, "04_program_delta_median_Z.pdf"),
  p,
  width = 7,
  height = 3.5
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "99_sessionInfo.txt")
)
