# E2-induced DRIP enrichment analysis

PROJECT_DIR <- "."

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ChIPseeker)
  library(GenomicRanges)
  library(IRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE81851")
drip_dir <- file.path(base_dir, "processed")
out_dir <- file.path(base_dir, "enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

project_dir <- PROJECT_DIR

rl65_file <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")

c5_file <- file.path(PROJECT_DIR, "gene_sets/C5_stable_up_genes.txt")

candidate24_file <- file.path(PROJECT_DIR, "gene_sets/candidate24_genes.txt")

drip_gene_file <- file.path(
  drip_dir,
  "03_E2DRIPgenesummary.csv"
)

required_files <- c(rl65_file, c5_file, candidate24_file, drip_gene_file)

if (any(!file.exists(required_files))) {
  stop(": \n",
       paste(required_files[!file.exists(required_files)], collapse = "\n"))
}

raw_files <- list.files(
  base_dir,
  pattern = "^GSE81851_DESeq_differential_calls\\.tsv$",
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

if (length(raw_files) > 0) {
  info <- file.info(raw_files)
  raw_files <- raw_files[!is.na(info$isdir) & !info$isdir]
}

if (length(raw_files) == 0) {
  stop("GSE81851_DESeq_differential_calls.tsv.")
}

raw_file <- raw_files[1]

clean_genes <- function(z) {

  z <- toupper(trimws(as.character(z)))
  z[z == "H2AFV"] <- "H2AZ2"
  z <- unique(z[!is.na(z) & nzchar(z)])

  return(z)
}

rl65_genes <- clean_genes(readLines(rl65_file, warn = FALSE))
c5_genes <- clean_genes(readLines(c5_file, warn = FALSE))
candidate24 <- clean_genes(readLines(candidate24_file, warn = FALSE))

final3 <- c("HMGB2", "H2AZ2", "UBALD2")

gene_sets <- list(
  RL_Sig65 = rl65_genes,
  C5_stable_up = c5_genes,
  Candidate_24 = candidate24
)

strict_df <- fread(drip_gene_file)
strict_df[, SYMBOL := clean_genes(SYMBOL)]

strict_genes <- clean_genes(
  strict_df$SYMBOL[strict_df$DRIP_any %in% TRUE]
)

x <- fread(raw_file)

needed <- c(
  "chr", "start", "end",
  "log2FoldChange.T0_T2_DRIP",
  "log2FoldChange.T0_T24_DRIP"
)

if (!all(needed %in% colnames(x))) {
  stop("GSE81851rawmissing.")
}

all_gr <- GRanges(
  seqnames = as.character(x$chr),
  ranges = IRanges(
    start = as.numeric(x$start) + 1,
    end = as.numeric(x$end)
  )
)

all_anno <- annotatePeak(
  all_gr,
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Hs.eg.db",
  verbose = FALSE
)

all_anno <- as.data.frame(all_anno)
all_anno$start_raw <- all_anno$start - 1

all_anno$Locus_class <- "Distal intergenic"
all_anno$Locus_class[grepl("^Promoter", all_anno$annotation)] <- "Promoter"
all_anno$Locus_class[grepl("Exon|Intron|UTR", all_anno$annotation,
                           ignore.case = TRUE)] <- "Gene body"
all_anno$Locus_class[grepl("^Downstream", all_anno$annotation)] <- "Downstream"

all_anno$Direct_locus <- all_anno$Locus_class %in% c("Promoter", "Gene body")

anno_key <- paste(all_anno$seqnames, all_anno$start_raw,
                  all_anno$end, sep = ":")

raw_key <- paste(x$chr, x$start, x$end, sep = ":")
idx <- match(anno_key, raw_key)

all_anno$LFC_2h <- x$log2FoldChange.T0_T2_DRIP[idx]
all_anno$LFC_24h <- x$log2FoldChange.T0_T24_DRIP[idx]

all_anno$SYMBOL <- toupper(all_anno$SYMBOL)
all_anno$SYMBOL[all_anno$SYMBOL == "H2AFV"] <- "H2AZ2"

direct_all <- as.data.table(all_anno)[
  Direct_locus == TRUE & !is.na(SYMBOL) & SYMBOL != ""
]

gene_rank <- direct_all[, .(
  N_DRIP_regions = .N,
  Mean_LFC_2h = mean(LFC_2h, na.rm = TRUE),
  Mean_LFC_24h = mean(LFC_24h, na.rm = TRUE)
), by = SYMBOL]

gene_rank[!is.finite(Mean_LFC_2h), Mean_LFC_2h := NA_real_]
gene_rank[!is.finite(Mean_LFC_24h), Mean_LFC_24h := NA_real_]

universe_genes <- unique(gene_rank$SYMBOL)
strict_genes <- intersect(strict_genes, universe_genes)

fisher_list <- lapply(names(gene_sets), function(set_name) {

  genes <- intersect(clean_genes(gene_sets[[set_name]]), universe_genes)
  overlap <- intersect(genes, strict_genes)

  a <- length(overlap)
  b <- length(genes) - a
  c <- length(strict_genes) - a
  d <- length(universe_genes) - a - b - c

  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  fit <- fisher.test(tab, alternative = "greater")

  data.frame(
    GeneSet = set_name,
    Input_size = length(gene_sets[[set_name]]),
    In_universe = length(genes),
    DRIP_overlap = a,
    Overlap_fraction = ifelse(length(genes) > 0, a / length(genes), NA),
    Background_fraction = length(strict_genes) / length(universe_genes),
    Odds_ratio = unname(fit$estimate),
    CI_low = fit$conf.int[1],
    CI_high = fit$conf.int[2],
    P_value = fit$p.value,
    stringsAsFactors = FALSE
  )
})

fisher_df <- rbindlist(fisher_list, fill = TRUE)
fisher_df[, FDR := p.adjust(P_value, method = "BH")]
fisher_df[, Analysis := "Fisher_ORA"]

run_rank_test <- function(set_name, genes, score_col, time_label) {

  genes <- intersect(clean_genes(genes), universe_genes)

  dat <- gene_rank[
    !is.na(get(score_col)) & is.finite(get(score_col))
  ]

  in_set <- dat$SYMBOL %in% genes
  score_set <- dat[[score_col]][in_set]
  score_bg <- dat[[score_col]][!in_set]

  if (length(score_set) < 3 || length(score_bg) < 3) {
    return(data.frame(
      GeneSet = set_name, Analysis = time_label,
      In_universe = length(score_set), Median_set = NA,
      Median_background = NA, Delta_median = NA, P_value = NA
    ))
  }

  fit <- wilcox.test(
    score_set,
    score_bg,
    alternative = "greater",
    exact = FALSE
  )

  data.frame(
    GeneSet = set_name,
    Analysis = time_label,
    In_universe = length(score_set),
    Median_set = median(score_set, na.rm = TRUE),
    Median_background = median(score_bg, na.rm = TRUE),
    Delta_median = median(score_set, na.rm = TRUE) -
      median(score_bg, na.rm = TRUE),
    P_value = fit$p.value,
    stringsAsFactors = FALSE
  )
}

rank_list <- list()

for (nm in names(gene_sets)) {

  rank_list[[paste0(nm, "_2h")]] <- run_rank_test(
    nm, gene_sets[[nm]], "Mean_LFC_2h", "Rank_2h"
  )

  rank_list[[paste0(nm, "_24h")]] <- run_rank_test(
    nm, gene_sets[[nm]], "Mean_LFC_24h", "Rank_24h"
  )
}

rank_df <- rbindlist(rank_list, fill = TRUE)

rank_df[, FDR := p.adjust(P_value, method = "BH"), by = Analysis]

fisher_out <- fisher_df[, .(
  GeneSet, Analysis, Input_size, In_universe,
  DRIP_overlap, Overlap_fraction, Background_fraction,
  Odds_ratio, CI_low, CI_high, P_value, FDR
)]

rank_out <- rank_df[, .(
  GeneSet, Analysis,
  Input_size = NA_integer_,
  In_universe,
  DRIP_overlap = NA_integer_,
  Overlap_fraction = NA_real_,
  Background_fraction = NA_real_,
  Odds_ratio = NA_real_,
  CI_low = NA_real_,
  CI_high = NA_real_,
  P_value, FDR,
  Median_set,
  Median_background,
  Delta_median
)]

result_all <- rbind(fisher_out, rank_out, fill = TRUE)

fwrite(
  result_all,
  file.path(out_dir, "01_gene_set_DRIP_enrichment_statistics.csv")
)

overlap_list <- lapply(names(gene_sets), function(nm) {

  hit <- sort(intersect(
    clean_genes(gene_sets[[nm]]),
    strict_genes
  ))

  data.frame(
    GeneSet = nm,
    Gene = hit,
    stringsAsFactors = FALSE
  )
})

final3_hit <- intersect(final3, strict_genes)

overlap_df <- rbindlist(
  c(overlap_list,
    list(data.frame(GeneSet = "Final_3",
                    Gene = final3_hit))),
  fill = TRUE
)

fwrite(
  overlap_df,
  file.path(out_dir, "02_gene_set_DRIP_overlap_genes.csv")
)

gene_rank[, Strict_E2_DRIP := SYMBOL %in% strict_genes]
gene_rank[, RL_Sig65 := SYMBOL %in% rl65_genes]
gene_rank[, C5_stable_up := SYMBOL %in% c5_genes]
gene_rank[, Candidate_24 := SYMBOL %in% candidate24]

setorder(gene_rank, -Mean_LFC_24h)

fwrite(
  gene_rank,
  file.path(out_dir, "03_gene_level_E2_DRIP_effects.csv")
)

plot_fisher <- copy(fisher_df)
plot_fisher[, Effect := log2(Odds_ratio)]
plot_fisher[, Label := paste0(
  DRIP_overlap, "/", In_universe
)]

plot_rank <- copy(rank_df)
plot_rank[, Effect := Delta_median]
plot_rank[, Label := sprintf("FDR=%.3g", FDR)]

plot_data <- rbind(
  plot_fisher[, .(
    GeneSet,
    Analysis = "Strict DRIP enrichment",
    Effect,
    FDR,
    Label
  )],
  plot_rank[, .(
    GeneSet,
    Analysis = ifelse(Analysis == "Rank_2h",
                      "E2-induced DRIP shift: 2 h",
                      "E2-induced DRIP shift: 24 h"),
    Effect,
    FDR,
    Label
  )]
)

plot_data$GeneSet <- factor(
  plot_data$GeneSet,
  levels = c("RL_Sig65", "C5_stable_up", "Candidate_24")
)

p <- ggplot(
  plot_data,
  aes(x = Effect, y = GeneSet)
) +
  geom_vline(xintercept = 0, linetype = 2,
             linewidth = 0.4, color = "grey60") +
  geom_point(aes(size = -log10(pmax(FDR, 1e-300))),
             shape = 21, fill = "white", stroke = 0.8) +
  geom_text(aes(label = Label), nudge_y = 0.22,
            size = 2.5, check_overlap = TRUE) +
  facet_wrap(~Analysis, scales = "free_x", ncol = 1) +
  labs(
    x = "Enrichment effect",
    y = NULL,
    size = expression(-log[10](FDR))
  ) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.text = element_text(size = 7),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "04_gene_set_DRIPenrichmentplot.pdf"),
  p,
  width = 8.5,
  height = 11,
  units = "cm"
)

method_text <- c(
  "GSE81851 direct DRIP enrichment analysis",
  "",
  paste0("DRIP background genes: ", length(universe_genes)),
  paste0("Strict E2-induced direct DRIP genes: ", length(strict_genes)),
  "",
  "Direct locus definition: promoter <=1 kb from TSS, exon, intron or UTR.",
  "Distal intergenic regions were excluded from direct gene-locus enrichment.",
  "",
  "Fisher analysis:",
  "One-sided Fisher exact test using genes directly annotatable in the GSE81851 DRIP-seq processed region universe as background.",
  "P values were adjusted by BH across the three tested gene sets.",
  "",
  "Rank-based analysis:",
  "For each directly annotated gene, mean E2-induced DRIP log2FC was calculated across its promoter/gene-body regions.",
  "One-sided Wilcoxon rank-sum tests evaluated whether each gene set showed higher DRIP changes than the remaining background genes.",
  "2 h and 24 h were tested separately, with BH correction within each time point.",
  "",
  "Final three genes were evaluated descriptively and were not subjected to gene-set enrichment because n=3."
)

writeLines(
  method_text,
  file.path(out_dir, "05_DRIPenrichmentanalysis.txt")
)

print(fisher_df[, .(
  GeneSet, In_universe, DRIP_overlap,
  Odds_ratio, P_value, FDR
)])

print(rank_df[, .(
  GeneSet, Analysis, Median_set,
  Median_background, Delta_median, P_value, FDR
)])

options(stringsAsFactors=FALSE)

if(!requireNamespace("patchwork",quietly=TRUE)) install.packages("patchwork")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE81851")
enrich_dir <- file.path(base_dir,"enrichment")
drip_dir <- file.path(base_dir,"processed")
project_dir <- PROJECT_DIR

stat_file <- file.path(enrich_dir,"01_gene_set_DRIP_enrichment_statistics.csv")
rank_file <- file.path(enrich_dir,"03_gene_level_E2_DRIP_effects.csv")
direct_file <- file.path(drip_dir,"03_E2DRIPgenesummary.csv")
candidate_file <- file.path(PROJECT_DIR, "gene_sets", "candidate24_genes.txt")

need <- c(stat_file,rank_file,direct_file,candidate_file)
if(any(!file.exists(need))) stop(": \n",paste(need[!file.exists(need)],collapse="\n"))

res <- fread(stat_file)
gene_rank <- fread(rank_file)
direct <- fread(direct_file)

norm_gene <- function(x){
  x <- toupper(trimws(as.character(x)))
  x[x=="H2AFV"] <- "H2AZ2"
  x
}

gene_rank[,SYMBOL:=norm_gene(SYMBOL)]
direct[,SYMBOL:=norm_gene(SYMBOL)]
candidate24 <- unique(norm_gene(readLines(candidate_file,warn=FALSE)))
candidate24 <- candidate24[!is.na(candidate24) & nzchar(candidate24)]
final3 <- c("HMGB2","H2AZ2","UBALD2")

col_set <- c(
  "RL-Sig65"="#b43665",
  "C5 stable-up"="#6fa6cf",
  "24 candidates"="#8560af"
)

fmt_fdr <- function(x){
  out <- rep("NA",length(x)); ok <- !is.na(x)
  out[ok & x<1e-99] <- "<1e-99"
  out[ok & x>=1e-99 & x<0.001] <- formatC(x[ok & x>=1e-99 & x<0.001],format="e",digits=1)
  out[ok & x>=0.001] <- sprintf("%.3f",x[ok & x>=0.001])
  out
}

theme_pub <- theme_classic(base_size=8) +
  theme(
    axis.text=element_text(size=7,colour="black"),
    axis.title=element_text(size=8,colour="black"),
    axis.line=element_line(linewidth=.4,colour="black"),
    axis.ticks=element_line(linewidth=.35,colour="black"),
    legend.text=element_text(size=7),
    legend.title=element_text(size=7),
    strip.background=element_blank(),
    strip.text=element_text(size=8,face="bold"),
    plot.title=element_text(size=9,face="bold",hjust=.5),
    plot.subtitle=element_text(size=6.8,colour="grey35",hjust=.5),
    plot.tag=element_text(size=10,face="bold"),
    plot.margin=margin(5,6,5,6)
  )

name_map <- c(
  RL_Sig65="RL-Sig65",
  C5_stable_up="C5 stable-up",
  Candidate_24="24 candidates"
)

fisher <- res[Analysis=="Fisher_ORA"]
fisher[,GeneSet_label:=name_map[GeneSet]]
fisher[,Support_pct:=100*DRIP_overlap/In_universe]
fisher[,Inside_label:=paste0(sprintf("%.1f",Support_pct),"%  ",
                             "(",DRIP_overlap,"/",In_universe,")")]
fisher[,Stat_label:=paste0("OR ",sprintf("%.2f",Odds_ratio),
                           "   FDR ",fmt_fdr(FDR))]

fisher$GeneSet_label <- factor(
  fisher$GeneSet_label,
  levels=c("24 candidates","C5 stable-up","RL-Sig65")
)

bg_pct <- 100*fisher$Background_fraction[1]

pA <- ggplot(fisher,aes(x=GeneSet_label)) +
  geom_col(aes(y=100),width=.62,fill="#e9e9e9") +
  geom_col(aes(y=Support_pct,fill=GeneSet_label),width=.62) +
  geom_hline(yintercept=bg_pct,linetype=2,linewidth=.45,colour="grey35") +
  geom_text(aes(y=Support_pct/2,label=Inside_label),
            colour="white",fontface="bold",size=2.55) +
  geom_text(aes(y=104,label=Stat_label),hjust=0,size=2.35,colour="grey20") +
  scale_fill_manual(values=col_set) +
  scale_y_continuous(
    breaks=c(0,25,50,75,100),
    limits=c(0,145),
    expand=c(0,0)
  ) +
  coord_flip(clip="off") +
  labs(
    x=NULL,
    y="Genes with direct E2-induced DRIP support (%)",
    title="Direct DRIP support",
    subtitle=paste0("Dashed line: DRIP background = ",sprintf("%.1f",bg_pct),"%"),
    tag="A"
  ) +
  theme_pub +
  theme(legend.position="none")

ggsave(
  file.path(base_dir,"GSE81851_A_DRIPfraction_publication.pdf"),
  pA,width=9.6,height=6.8,units="cm",device="pdf"
)

make_dist <- function(flag,label){
  d <- gene_rank[get(flag)%in%TRUE,.(SYMBOL,Mean_LFC_2h,Mean_LFC_24h)]
  d[,GeneSet:=label]
  d
}

dist <- rbind(
  make_dist("RL_Sig65","RL-Sig65"),
  make_dist("C5_stable_up","C5 stable-up"),
  make_dist("Candidate_24","24 candidates")
)

dist <- melt(
  dist,
  id.vars=c("SYMBOL","GeneSet"),
  measure.vars=c("Mean_LFC_2h","Mean_LFC_24h"),
  variable.name="Time",
  value.name="DRIP_LFC"
)

dist[Time=="Mean_LFC_2h",Time:="2 h"]
dist[Time=="Mean_LFC_24h",Time:="24 h"]
dist <- dist[is.finite(DRIP_LFC)]

dist$GeneSet <- factor(
  dist$GeneSet,
  levels=c("RL-Sig65","C5 stable-up","24 candidates")
)
dist$Time <- factor(dist$Time,levels=c("2 h","24 h"))

bg_med <- data.table(
  Time=factor(c("2 h","24 h"),levels=c("2 h","24 h")),
  Median=c(
    median(gene_rank$Mean_LFC_2h,na.rm=TRUE),
    median(gene_rank$Mean_LFC_24h,na.rm=TRUE)
  )
)

n_lab <- dist[,.(N=uniqueN(SYMBOL)),by=.(Time,GeneSet)]

pB <- ggplot(dist,aes(x=GeneSet,y=DRIP_LFC,fill=GeneSet)) +
  geom_hline(yintercept=0,linewidth=.35,colour="grey70") +
  geom_hline(
    data=bg_med,
    aes(yintercept=Median),
    inherit.aes=FALSE,
    linetype=2,
    linewidth=.45,
    colour="grey30"
  ) +
  geom_violin(
    trim=TRUE,
    scale="width",
    width=.82,
    alpha=.78,
    linewidth=.35,
    colour="grey30"
  ) +
  geom_boxplot(
    width=.14,
    outlier.shape=NA,
    fill="white",
    linewidth=.42,
    colour="grey20"
  ) +
  geom_text(
    data=n_lab,
    aes(x=GeneSet,y=Inf,label=paste0("n=",N)),
    inherit.aes=FALSE,
    vjust=1.25,
    size=2.2,
    colour="grey25"
  ) +
  facet_wrap(~Time,nrow=1,scales="free_y") +
  scale_fill_manual(values=col_set) +
  labs(
    x=NULL,
    y="E2-induced DRIP log2 fold change",
    title="E2-induced DRIP dynamics",
    subtitle="Dashed line indicates the genome-wide DRIP background median",
    tag="B"
  ) +
  theme_pub +
  theme(
    legend.position="none",
    axis.text.x=element_text(angle=18,hjust=1)
  )

ggsave(
  file.path(base_dir,"GSE81851_B_DRIP_publication.pdf"),
  pB,width=10.4,height=7.2,units="cm",device="pdf"
)

cand <- data.table(SYMBOL=candidate24)

rank_sub <- gene_rank[,.(SYMBOL,Mean_LFC_2h,Mean_LFC_24h)]
direct_sub <- direct[,.(SYMBOL,DRIP_2h,DRIP_24h,DRIP_persistent)]

cand <- merge(cand,rank_sub,by="SYMBOL",all.x=TRUE)
cand <- merge(cand,direct_sub,by="SYMBOL",all.x=TRUE)

for(v in c("DRIP_2h","DRIP_24h","DRIP_persistent")){
  cand[is.na(get(v)),(v):=FALSE]
}

cand[,Final3:=SYMBOL%in%final3]
cand[,Evaluable:=is.finite(Mean_LFC_2h)|is.finite(Mean_LFC_24h)]

setorder(cand,-Final3,-DRIP_persistent,-DRIP_24h,-Mean_LFC_24h)

cand[,Gene_label:=ifelse(Final3,paste0(SYMBOL,"  "),SYMBOL)]
cand$Gene_label <- factor(cand$Gene_label,levels=rev(cand$Gene_label))

heat <- rbind(
  cand[,.(Gene_label,X=1,LFC=Mean_LFC_2h)],
  cand[,.(Gene_label,X=2,LFC=Mean_LFC_24h)]
)

support <- rbind(
  cand[,.(Gene_label,X=3,Supported=DRIP_2h,Evaluable)],
  cand[,.(Gene_label,X=4,Supported=DRIP_24h,Evaluable)],
  cand[,.(Gene_label,X=5,Supported=DRIP_persistent,Evaluable)]
)

support[,Status:=fifelse(!Evaluable,"Not evaluable",
                         fifelse(Supported,"Supported","Not supported"))]

support$Status <- factor(
  support$Status,
  levels=c("Supported","Not supported","Not evaluable")
)

lim <- max(abs(heat$LFC),na.rm=TRUE)
if(!is.finite(lim) || lim==0) lim <- 1

pC <- ggplot() +
  geom_tile(
    data=heat,
    aes(x=X,y=Gene_label,fill=LFC),
    width=.82,
    height=.82,
    colour="white",
    linewidth=.35
  ) +
  geom_point(
    data=support,
    aes(x=X,y=Gene_label,shape=Status),
    size=2.55,
    stroke=.65,
    colour="#187d79"
  ) +
  scale_fill_gradient2(
    low="#6fa6cf",
    mid="#f7f7f7",
    high="#b43665",
    midpoint=0,
    limits=c(-lim,lim),
    oob=scales::squish,
    name="E2-induced\nDRIP log2FC"
  ) +
  scale_shape_manual(
    values=c(
      "Supported"=16,
      "Not supported"=1,
      "Not evaluable"=4
    ),
    name="Direct DRIP"
  ) +
  scale_x_continuous(
    breaks=1:5,
    labels=c(
      "2 h\nchange",
      "24 h\nchange",
      "2 h\nDRIP",
      "24 h\nDRIP",
      "Persistent"
    ),
    limits=c(.5,5.45),
    expand=c(0,0)
  ) +
  labs(
    x=NULL,
    y=NULL,
    title="DRIP landscape of the 24 candidate genes",
    subtitle=" Final three candidates",
    tag="C"
  ) +
  theme_minimal(base_size=8) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(size=6.7,face="bold",colour="black"),
    axis.text.y=element_text(size=6.5,face="italic",colour="black"),
    axis.ticks=element_blank(),
    legend.position="bottom",
    legend.box="vertical",
    legend.title=element_text(size=7),
    legend.text=element_text(size=6.5),
    plot.title=element_text(size=9,face="bold",hjust=.5),
    plot.subtitle=element_text(size=6.8,colour="grey35",hjust=.5),
    plot.tag=element_text(size=10,face="bold"),
    plot.margin=margin(5,5,5,5)
  )

ggsave(
  file.path(base_dir,"GSE81851_C_24candidategeneDRIP_publication.pdf"),
  pC,width=8.0,height=13.7,units="cm",device="pdf"
)

left_panel <- pA / pB + plot_layout(heights=c(.82,1.18))

final_fig <- (left_panel | pC) +
  plot_layout(widths=c(1.35,.85))

ggsave(
  file.path(base_dir,"GSE81851_DRIP_publication_V2.pdf"),
  final_fig,
  width=17,
  height=14.3,
  units="cm",
  device="pdf"
)
