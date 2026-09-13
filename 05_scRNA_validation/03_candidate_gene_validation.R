# Candidate-gene validation in GSE245601

PROJECT_DIR <- "."

options(stringsAsFactors=FALSE); set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tidyr)
  library(edgeR); library(ggplot2); library(ggrepel)
  library(patchwork); library(pheatmap)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")
input_file <- file.path(base_dir, "10.malignat",
                        "GSE245601_Malignant_Epithelial_Reclustered.rds")
candidate_file <- file.path(PROJECT_DIR, "gene_sets/candidate24_genes.txt")
out_dir <- file.path(base_dir, "12.gene")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

FM_MIN_PCT <- 0.01; FM_FDR_CUT <- 0.05; FM_LFC_CUT <- 0.25
PB_FDR_CUT <- 0.05; PB_LFC_CUT <- 0.25
PB_MIN_C4 <- 10L; PB_MIN_OTHER <- 50L
DESC_MIN_C4 <- 3L; CLUSTER_MIN_CELLS <- 10L

for (f in c(input_file, candidate_file)) {
  if (!file.exists(f)) stop("not_found: ", f)
}

safe_vec <- function(x) {
  if (is.data.frame(x)) x <- x[[1]]
  if (is.list(x)) x <- unlist(x, recursive=TRUE, use.names=FALSE)
  as.character(x)
}

fmt_fdr <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

fm_table <- function(x, prefix) {
  lfc <- intersect(c("avg_log2FC", "avg_logFC"), colnames(x))[1]
  if (is.na(lfc)) stop("FindMarkersresultsmissinglogFC")
  z <- data.frame(
    gene=rownames(x), log2FC=x[[lfc]], P=x$p_val, FDR=x$p_val_adj,
    pct1=x$pct.1, pct2=x$pct.2, row.names=NULL
  )
  colnames(z)[-1] <- paste0(prefix, "_", colnames(z)[-1])
  z
}

build_state_pb <- function(counts, patient, state, patients, prefix) {
  keep <- patient %in% patients
  info <- expand_grid(Patient=patients, State=c("Other","C4")) %>%
    mutate(key=paste(Patient, State, sep="||"),
           PB_id=sprintf(paste0(prefix,"%03d"), row_number()))
  id_map <- setNames(info$PB_id, info$key)
  id <- unname(id_map[paste(patient[keep], state[keep], sep="||")])
  f <- factor(id, levels=info$PB_id)
  ind <- Matrix::sparse.model.matrix(~0+f); colnames(ind) <- levels(f)
  mat <- as.matrix(counts[,keep,drop=FALSE] %*% ind)
  list(mat=mat[,info$PB_id,drop=FALSE], info=info)
}

pb_long <- function(pb, genes) {
  y <- DGEList(counts=pb$mat); y <- calcNormFactors(y)
  x <- cpm(y, log=TRUE, prior.count=2)
  as.data.frame(t(x[genes,,drop=FALSE])) %>%
    mutate(PB_id=rownames(.)) %>%
    left_join(pb$info, by="PB_id") %>%
    pivot_longer(all_of(genes), names_to="gene", values_to="logCPM")
}

theme_pub <- function(size=11) {
  theme_classic(base_size=size) +
    theme(
      plot.title=element_text(face="bold", hjust=0),
      plot.subtitle=element_text(color="grey30"),
      strip.background=element_blank(),
      strip.text=element_text(face="bold", size=size),
      axis.title=element_text(face="bold"),
      legend.title=element_blank()
    )
}

candidate_raw <- unique(trimws(readLines(candidate_file, warn=FALSE)))
candidate_raw <- candidate_raw[nzchar(candidate_raw)]

obj <- readRDS(input_file); DefaultAssay(obj) <- "RNA"
if (!"seurat_clusters" %in% colnames(obj[[]])) stop("missingseurat_clusters")
if (length(Layers(obj[["RNA"]])) > 1) obj <- JoinLayers(obj, assay="RNA")
if (!"data" %in% Layers(obj[["RNA"]])) obj <- NormalizeData(obj, verbose=FALSE)

gene_map <- setNames(rownames(obj), toupper(rownames(obj)))
candidate_use <- unique(unname(gene_map[toupper(candidate_raw)]))
candidate_use <- candidate_use[!is.na(candidate_use)]
if (!length(candidate_use)) stop("candidategene")

cluster_vec <- safe_vec(obj$seurat_clusters)
cluster_levels <- as.character(sort(unique(as.integer(cluster_vec))))
if (!"4" %in% cluster_levels) stop("4")

patient_col <- if ("patient_id" %in% colnames(obj[[]])) {
  "patient_id"
} else if ("sample_gsm" %in% colnames(obj[[]])) {
  "sample_gsm"
} else "orig.ident"

patient_vec <- safe_vec(obj@meta.data[[patient_col]])
all_patients <- sort(unique(patient_vec))
patient_map <- data.frame(
  Patient=all_patients, DisplayID=paste0("P",seq_along(all_patients))
)

obj$seurat_clusters <- factor(cluster_vec, levels=cluster_levels)
obj$C4_state <- factor(ifelse(cluster_vec=="4","C4","Other"),
                       levels=c("Other","C4"))

write.csv(patient_map, file.path(out_dir,"00_patient_map.csv"),
          row.names=FALSE)
writeLines(setdiff(toupper(candidate_raw),toupper(candidate_use)),
           file.path(out_dir,"00_candidategene.txt"))

print(table(obj$seurat_clusters))

fm_main <- FindMarkers(
  obj, ident.1="C4", ident.2="Other", group.by="C4_state",
  assay="RNA", test.use="wilcox", min.pct=FM_MIN_PCT,
  logfc.threshold=0, only.pos=FALSE, verbose=FALSE
) %>% fm_table("FM") %>% arrange(FM_FDR)

write.csv(fm_main, file.path(out_dir,"01_C4_vs_Other_FindMarkersresults.csv"),
          row.names=FALSE)

pair_df <- lapply(setdiff(cluster_levels,"4"), function(cl) {
  x <- FindMarkers(
    obj, ident.1="4", ident.2=cl, group.by="seurat_clusters",
    assay="RNA", test.use="wilcox", min.pct=FM_MIN_PCT,
    logfc.threshold=0, only.pos=FALSE, verbose=FALSE
  ) %>% fm_table("Pair")
  x$ref_cluster <- cl; x
}) %>% bind_rows()

pair_24 <- pair_df %>% filter(gene %in% candidate_use)
pair_stats <- pair_24 %>%
  group_by(gene) %>%
  summarise(
    pair_n=n(), pair_n_up=sum(Pair_log2FC>0,na.rm=TRUE),
    pair_frac_up=pair_n_up/pair_n,
    pair_median_log2FC=median(Pair_log2FC,na.rm=TRUE),
    pair_min_log2FC=min(Pair_log2FC,na.rm=TRUE), .groups="drop"
  )

write.csv(pair_24, file.path(out_dir,"02_24gene_C4clusterFindMarkersresults.csv"),
          row.names=FALSE)

state_count <- as.data.frame(
  table(Patient=patient_vec, State=obj$C4_state),
  stringsAsFactors=FALSE
) %>%
  pivot_wider(names_from=State, values_from=Freq, values_fill=0) %>%
  mutate(
    MainIncluded=C4>=PB_MIN_C4 & Other>=PB_MIN_OTHER,
    DescriptiveIncluded=C4>=DESC_MIN_C4 & Other>=PB_MIN_OTHER
  ) %>%
  left_join(patient_map, by="Patient")

main_patients <- state_count$Patient[state_count$MainIncluded]
desc_patients <- state_count$Patient[state_count$DescriptiveIncluded]
if (length(main_patients)<4) stop("bulk4")

write.csv(state_count, file.path(out_dir,"03_C4Othercell.csv"),
          row.names=FALSE)

counts_mat <- GetAssayData(obj, assay="RNA", layer="counts")
state_vec <- as.character(obj$C4_state)
pb_main <- build_state_pb(counts_mat, patient_vec, state_vec,
                          main_patients, "PB")

meta <- pb_main$info %>%
  mutate(
    Patient=factor(Patient,levels=main_patients),
    State=factor(State,levels=c("Other","C4"))
  )

y <- DGEList(counts=pb_main$mat)
y <- y[filterByExpr(y,group=meta$State),,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~Patient+State,data=meta)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design,robust=TRUE)
qlf <- glmQLFTest(fit,coef="StateC4")

pb_de <- topTags(qlf,n=Inf)$table %>%
  data.frame() %>% tibble::rownames_to_column("gene") %>%
  rename(PB_log2FC=logFC,PB_logCPM=logCPM,PB_P=PValue,PB_FDR=FDR) %>%
  arrange(PB_FDR)

write.csv(pb_de, file.path(out_dir,"04_C4_vs_Other_bulk_edgeRresults.csv"),
          row.names=FALSE)

paired_main <- pb_long(pb_main,candidate_use)
pb_desc <- build_state_pb(counts_mat,patient_vec,state_vec,desc_patients,"PD")
paired_desc <- pb_long(pb_desc,candidate_use)

get_stats <- function(x,prefix) {
  z <- x %>%
    select(Patient,State,gene,logCPM) %>%
    pivot_wider(names_from=State,values_from=logCPM) %>%
    mutate(delta=C4-Other) %>%
    group_by(gene) %>%
    summarise(
      n=n(), n_C4_higher=sum(delta>0,na.rm=TRUE),
      frac_C4_higher=n_C4_higher/n,
      median_delta=median(delta,na.rm=TRUE), .groups="drop"
    )
  colnames(z)[-1] <- paste0(prefix,"_",colnames(z)[-1]); z
}

main_stats <- get_stats(paired_main,"main")
desc_stats <- get_stats(paired_desc,"desc")

write.csv(paired_desc, file.path(out_dir,"05_24geneexpression.csv"),
          row.names=FALSE)

summary_df <- data.frame(gene=candidate_use) %>%
  left_join(fm_main,by="gene") %>%
  left_join(pair_stats,by="gene") %>%
  left_join(pb_de,by="gene") %>%
  left_join(main_stats,by="gene") %>%
  left_join(desc_stats,by="gene") %>%
  mutate(
    FM_up=!is.na(FM_FDR) & FM_FDR<FM_FDR_CUT &
      FM_log2FC>=FM_LFC_CUT,
    PB_up=!is.na(PB_FDR) & PB_FDR<PB_FDR_CUT &
      PB_log2FC>=PB_LFC_CUT,

    direction_same=!is.na(FM_log2FC) & !is.na(PB_log2FC) &
      FM_log2FC>0 & PB_log2FC>0,
    Tier=case_when(
      FM_up & PB_up & pair_frac_up>=0.75 &
        main_frac_C4_higher>=0.70 ~ "Tier 1",
      direction_same & (FM_up|PB_up) & pair_frac_up>=0.50 &
        main_frac_C4_higher>=0.60 ~ "Tier 2",
      TRUE ~ "Tier 3"
    )
  ) %>%
  arrange(factor(Tier,levels=c("Tier 1","Tier 2","Tier 3")),
          desc(PB_log2FC),desc(FM_log2FC))

tier1_genes <- summary_df$gene[summary_df$Tier=="Tier 1"]
if (!length(tier1_genes)) stop("Tier 1gene")

write.csv(summary_df, file.path(out_dir,"06_24genesummary.csv"),
          row.names=FALSE)
for (tier in c("Tier 1","Tier 2","Tier 3")) {
  writeLines(summary_df$gene[summary_df$Tier==tier],
             file.path(out_dir,paste0("07_",gsub(" ","",tier),"_genes.txt")))
}

scatter_df <- summary_df %>%
  mutate(
    Validation=factor(Tier,levels=c("Tier 3","Tier 2","Tier 1")),
    Label=ifelse(Tier=="Tier 1",gene,NA_character_)
  )

pA <- ggplot(scatter_df,aes(FM_log2FC,PB_log2FC)) +
  geom_hline(yintercept=0,linetype=2,color="grey65") +
  geom_vline(xintercept=0,linetype=2,color="grey65") +
  geom_hline(yintercept=PB_LFC_CUT,linetype=3,color="grey78") +
  geom_vline(xintercept=FM_LFC_CUT,linetype=3,color="grey78") +
  geom_abline(slope=1,linetype=3,color="grey82") +
  geom_point(aes(fill=Validation),shape=21,size=3.7,
             color="black",stroke=0.35) +
  geom_text_repel(aes(label=Label),size=4,fontface="bold",
                  box.padding=0.5,point.padding=0.4,
                  min.segment.length=0,seed=1234) +
  scale_fill_manual(values=c(
    "Tier 1"="#B2182B","Tier 2"="#E69F00","Tier 3"="#D9D9D9"
  )) +
  coord_equal() + theme_pub(12) +
  theme(legend.position="top") +
  labs(
    title="Cross-method validation of candidate genes",
    subtitle=paste0(length(tier1_genes)," of ",nrow(summary_df),
                    " genes met all validation criteria"),
    x="Single-cell FindMarkers log2 fold-change",
    y="Patient-level pseudobulk log2 fold-change",
    fill=NULL
  )

ggsave(file.path(out_dir,"08A_24gene_consistencyplot.pdf"),
       pA,width=7.2,height=6.4)

effect_order <- summary_df %>% arrange(PB_log2FC) %>% pull(gene)
effect_df <- summary_df %>% mutate(gene=factor(gene,levels=effect_order))
effect_long <- bind_rows(
  summary_df %>% transmute(gene,Method="FindMarkers",Effect=FM_log2FC),
  summary_df %>% transmute(gene,Method="Pseudobulk edgeR",Effect=PB_log2FC)
) %>% mutate(
  gene=factor(gene,levels=effect_order),
  Method=factor(Method,levels=c("FindMarkers","Pseudobulk edgeR"))
)

gene_labels <- setNames(
  ifelse(effect_order %in% tier1_genes,paste0(effect_order,"  *"),effect_order),
  effect_order
)

pB <- ggplot(effect_df,aes(y=gene)) +
  geom_vline(xintercept=0,linetype=2,color="grey65") +
  geom_vline(xintercept=0.25,linetype=3,color="grey78") +
  geom_segment(aes(x=FM_log2FC,xend=PB_log2FC,yend=gene),
               color="grey75",linewidth=0.55) +
  geom_point(data=effect_long,aes(x=Effect,fill=Method),
             shape=21,size=3,color="black",stroke=0.3) +
  scale_fill_manual(values=c(
    "FindMarkers"="#4C78A8","Pseudobulk edgeR"="#E45756"
  )) +
  scale_y_discrete(labels=gene_labels) +
  theme_pub(10.5) +
  theme(
    legend.position="top",axis.title.y=element_blank(),
    panel.grid.major.y=element_line(color="grey94",linewidth=0.3)
  ) +
  labs(
    title="Effect-size agreement across 24 candidate genes",
    subtitle="* Tier 1 externally validated genes",
    x="C4-like versus non-C4 log2 fold-change",
    fill=NULL
  )

ggsave(file.path(out_dir,"08B_24gene_plot.pdf"),
       pB,width=7.8,height=9)

delta_df <- paired_desc %>%
  select(Patient,State,gene,logCPM) %>%
  pivot_wider(names_from=State,values_from=logCPM) %>%
  mutate(delta=C4-Other) %>%
  left_join(patient_map,by="Patient") %>%
  select(gene,DisplayID,delta)

delta_order <- summary_df %>%
  arrange(factor(Tier,levels=c("Tier 1","Tier 2","Tier 3")),
          desc(desc_median_delta)) %>% pull(gene)

delta_mat <- delta_df %>%
  pivot_wider(names_from=DisplayID,values_from=delta) %>% as.data.frame()
rownames(delta_mat) <- delta_mat$gene; delta_mat$gene <- NULL
delta_mat <- as.matrix(delta_mat[delta_order,,drop=FALSE])

heat_lim <- as.numeric(quantile(abs(delta_mat),0.95,na.rm=TRUE))
if (!is.finite(heat_lim) || heat_lim==0) heat_lim <- 1
delta_plot <- pmax(pmin(delta_mat,heat_lim),-heat_lim)

ann_row <- data.frame(
  Tier=factor(summary_df$Tier[match(rownames(delta_mat),summary_df$gene)],
              levels=c("Tier 1","Tier 2","Tier 3"))
)
rownames(ann_row) <- rownames(delta_mat)
run_tier <- rle(as.character(ann_row$Tier))
gaps_row <- head(cumsum(run_tier$lengths),-1)

pdf(file.path(out_dir,"08C_24gene_plot.pdf"),
    width=8.5,height=9.5)
pheatmap(
  delta_plot,cluster_rows=FALSE,cluster_cols=FALSE,
  color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  breaks=seq(-heat_lim,heat_lim,length.out=101),
  border_color="white",fontsize_row=9,fontsize_col=10,
  annotation_row=ann_row,
  annotation_colors=list(Tier=c(
    "Tier 1"="#B2182B","Tier 2"="#E69F00","Tier 3"="#BDBDBD"
  )),
  gaps_row=gaps_row,
  main="Patient-level expression shift in C4-like cells\nC4 minus non-C4 pseudobulk logCPM"
)
dev.off()

paired_tier1 <- paired_main %>%
  filter(gene %in% tier1_genes) %>%
  left_join(patient_map,by="Patient") %>%
  mutate(
    gene=factor(gene,levels=tier1_genes),
    State=factor(State,levels=c("Other","C4"),
                 labels=c("non-C4","C4-like"))
  )

facet_info <- summary_df %>%
  filter(gene %in% tier1_genes) %>%
  mutate(facet=paste0(
    gene,"\nlog2FC=",sprintf("%.2f",PB_log2FC),
    ", FDR ",fmt_fdr(PB_FDR)
  ))
facet_labels <- setNames(facet_info$facet,facet_info$gene)

pD <- ggplot(paired_tier1,aes(State,logCPM,group=DisplayID)) +
  geom_line(color="grey65",linewidth=0.55,alpha=0.85) +
  geom_point(aes(fill=State),shape=21,size=3,
             color="black",stroke=0.35) +
  stat_summary(aes(group=State),fun=median,geom="crossbar",
               width=0.45,linewidth=0.7,color="black") +
  facet_wrap(~gene,scales="free_y",nrow=1,
             labeller=as_labeller(facet_labels)) +
  scale_fill_manual(values=c("non-C4"="#D9D9D9","C4-like"="#B2182B")) +
  theme_pub(12) +
  theme(legend.position="none",axis.text.x=element_text(face="bold")) +
  labs(
    title="Patient-level validation of Tier 1 genes",
    subtitle=paste0("Formal paired pseudobulk analysis; n = ",
                    length(main_patients)," patients"),
    x=NULL,y="Pseudobulk logCPM"
  )

ggsave(file.path(out_dir,"08D_Tier1gene_plot.pdf"),
       pD,width=10.5,height=4.8)

cluster_count <- as.data.frame(
  table(Patient=patient_vec,Cluster=cluster_vec),
  stringsAsFactors=FALSE
) %>%
  filter(Freq>=CLUSTER_MIN_CELLS) %>%
  mutate(key=paste(Patient,Cluster,sep="||"),
         PB_id=sprintf("PC%03d",row_number()))

cluster_map <- setNames(cluster_count$PB_id,cluster_count$key)
cell_id <- unname(cluster_map[paste(patient_vec,cluster_vec,sep="||")])
keep_cell <- !is.na(cell_id)
f <- factor(cell_id[keep_cell],levels=cluster_count$PB_id)
ind <- Matrix::sparse.model.matrix(~0+f); colnames(ind) <- levels(f)

cluster_mat <- as.matrix(counts_mat[,keep_cell,drop=FALSE] %*% ind)
cluster_mat <- cluster_mat[,cluster_count$PB_id,drop=FALSE]
yc <- DGEList(counts=cluster_mat); yc <- calcNormFactors(yc)
cluster_logcpm <- cpm(yc,log=TRUE,prior.count=2)

cluster_long <- as.data.frame(
  t(cluster_logcpm[tier1_genes,,drop=FALSE])
) %>%
  mutate(PB_id=rownames(.)) %>%
  left_join(
    cluster_count %>% select(PB_id,Patient,Cluster,CellNumber=Freq),
    by="PB_id"
  ) %>%
  left_join(patient_map,by="Patient") %>%
  pivot_longer(all_of(tier1_genes),names_to="gene",values_to="logCPM") %>%
  mutate(
    Cluster=factor(Cluster,levels=cluster_levels),
    gene=factor(gene,levels=tier1_genes),
    Highlight=ifelse(Cluster=="4","C4-like","Other clusters")
  )

write.csv(cluster_long,
          file.path(out_dir,"09_Tier1gene_clusterbulkexpression.csv"),
          row.names=FALSE)

pE <- ggplot(cluster_long,aes(Cluster,logCPM)) +
  geom_boxplot(aes(fill=Highlight),width=0.62,outlier.shape=NA,
               alpha=0.88,linewidth=0.45) +
  geom_jitter(aes(fill=Highlight),shape=21,width=0.11,size=2.15,
              color="black",stroke=0.25,alpha=0.9) +
  facet_wrap(~gene,scales="free_y",nrow=1) +
  scale_fill_manual(values=c("C4-like"="#B2182B",
                             "Other clusters"="#D9D9D9")) +
  theme_pub(12) +
  theme(legend.position="top") +
  labs(
    title="Tier 1 genes across malignant epithelial clusters",
    subtitle=paste0(
      "Each point represents one patient x cluster pseudobulk; minimum ",
      CLUSTER_MIN_CELLS," cells"
    ),
    x="Validation cluster",y="Pseudobulk logCPM",fill=NULL
  )

ggsave(file.path(out_dir,"08E_Tier1gene_0-4clusterplot.pdf"),
       pE,width=10.5,height=4.8)

main_fig <- (pA+pB+plot_layout(widths=c(1,1.05))) /
  (pD+pE+plot_layout(widths=c(1,1))) +
  plot_annotation(
    tag_levels="A",
    title="External validation of C5-related candidate genes",
    theme=theme(plot.title=element_text(face="bold",size=16))
  )

ggsave(file.path(out_dir,"08F_plot_candidategene.pdf"),
       main_fig,width=15,height=11)

saveRDS(
  list(
    candidate_genes=candidate_use,tier1_genes=tier1_genes,
    FindMarkers=fm_main,pairwise_FindMarkers=pair_24,
    pseudobulk_DE=pb_de,paired_main=paired_main,
    paired_descriptive=paired_desc,cluster_expression=cluster_long,
    summary=summary_df,main_patients=main_patients,
    descriptive_patients=desc_patients
  ),
  file.path(out_dir,"12A_24gene_validation_publication_bundle.rds")
)

options(stringsAsFactors = FALSE)
set.seed(1234)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

base_dir <- file.path(PROJECT_DIR, "results/GSE245601")

old_dir <- file.path(base_dir, "12.gene")
out_dir <- file.path(base_dir, "12.gene", "2")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

bundle_file <- file.path(
  old_dir,
  "12A_24gene_validation_publication_bundle.rds"
)

patient_map_file <- file.path(
  old_dir,
  "00_patient_map.csv"
)

for (f in c(bundle_file, patient_map_file)) {
  if (!file.exists(f)) {
    stop("not_found: ", f)
  }
}

col_rose   <- "#b43665"
col_blue   <- "#6fa6cf"
col_purple <- "#8560af"
col_yellow <- "#efb421"
col_orange <- "#e58027"
col_teal   <- "#187d79"

col_grey1 <- "#E4E4E4"
col_grey2 <- "#B8B8B8"
col_grey3 <- "#6A6A6A"
col_black <- "#202020"

col_region <- "#F5DDE5"

fmt_fdr <- function(x) {

  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x < 0.001,
      "<0.001",
      sprintf("%.3f", x)
    )
  )
}

theme_pub <- function(base_size = 11) {

  theme_classic(base_size = base_size) +

    theme(

      plot.title = element_blank(),
      plot.subtitle = element_blank(),

      axis.title = element_text(
        face = "bold",
        color = col_black
      ),

      axis.text = element_text(
        color = col_black
      ),

      axis.text.x = element_text(
        face = "bold"
      ),

      strip.background = element_blank(),

      strip.text = element_text(
        face = "bold",
        color = col_black,
        size = base_size
      ),

      legend.title = element_blank(),

      legend.text = element_text(
        color = col_black
      ),

      axis.line = element_line(
        linewidth = 0.55,
        color = col_black
      ),

      axis.ticks = element_line(
        linewidth = 0.5,
        color = col_black
      ),

      plot.margin = margin(
        7, 9, 7, 9
      )
    )
}

res <- readRDS(bundle_file)

patient_map <- read.csv(
  patient_map_file,
  stringsAsFactors = FALSE
)

summary_df <- res$summary
validated_genes <- res$tier1_genes
paired_main <- res$paired_main
cluster_long <- res$cluster_expression
main_patients <- res$main_patients

if (length(validated_genes) < 1) {
  stop("not_foundfinalgene")
}

scatter_df <- summary_df %>%

  mutate(

    Validation = ifelse(
      gene %in% validated_genes,
      "Passed all validation criteria",
      "Other candidates"
    ),

    Validation = factor(
      Validation,
      levels = c(
        "Other candidates",
        "Passed all validation criteria"
      )
    ),

    Label = ifelse(
      gene %in% validated_genes,
      gene,
      NA_character_
    )
  )

x_min <- min(
  scatter_df$FM_log2FC,
  na.rm = TRUE
)

x_max <- max(
  scatter_df$FM_log2FC,
  na.rm = TRUE
)

y_min <- min(
  scatter_df$PB_log2FC,
  na.rm = TRUE
)

y_max <- max(
  scatter_df$PB_log2FC,
  na.rm = TRUE
)

pad_x <- max(
  0.25,
  diff(range(scatter_df$FM_log2FC, na.rm = TRUE)) * 0.10
)

pad_y <- max(
  0.25,
  diff(range(scatter_df$PB_log2FC, na.rm = TRUE)) * 0.10
)

x_left <- floor(
  (x_min - pad_x) * 10
) / 10

x_right <- ceiling(
  (x_max + pad_x) * 10
) / 10

y_low <- floor(
  (y_min - pad_y) * 10
) / 10

y_high <- ceiling(
  (y_max + pad_y) * 10
) / 10

criteria_text <- paste0(
  "Genes highlighted in rose passed all criteria:\n",
  "FDR < 0.05 and log2FC \u2265 0.25 in both analyses\n",
  "Positive in \u226575% of C4-vs-cluster contrasts\n",
  "C4 higher in \u226570% of evaluable patients"
)

pA <- ggplot(
  scatter_df,
  aes(FM_log2FC, PB_log2FC)
) +

  annotate(
    "rect",

    xmin = 0.25,
    xmax = x_right,

    ymin = 0.25,
    ymax = y_high,

    fill = col_region,
    alpha = 0.55,
    color = NA
  ) +

  geom_vline(
    xintercept = 0.25,
    linetype = 2,
    linewidth = 0.55,
    color = col_grey3
  ) +

  geom_hline(
    yintercept = 0.25,
    linetype = 2,
    linewidth = 0.55,
    color = col_grey3
  ) +

  geom_point(
    aes(fill = Validation),

    shape = 21,
    size = 4.1,

    color = "black",
    stroke = 0.40
  ) +

  geom_text_repel(

    aes(label = Label),

    size = 3.9,
    fontface = "bold",

    box.padding = 0.45,
    point.padding = 0.40,

    min.segment.length = 0,

    segment.color = col_grey3,
    segment.size = 0.35,

    seed = 1234
  ) +

  annotate(

    "label",

    x = x_right -
      0.03 * (x_right - x_left),

    y = y_high -
      0.07 * (y_high - y_low),

    label = paste0(
      "Effect-size concordance region\n",
      "Both log2FC \u2265 0.25"
    ),

    hjust = 1,
    vjust = 1,

    size = 3.4,
    fontface = "bold",

    label.size = 0.20,

    fill = "white",
    color = col_rose
  ) +

  annotate(

    "label",

    x = 0.25,

    y = y_high -
      0.38 * (y_high - y_low),

    label = paste0(
      "Single-cell cutoff\n",
      "log2FC = 0.25"
    ),

    hjust = -0.05,
    vjust = 1,

    size = 3.0,

    label.size = 0.18,

    fill = "white",
    color = col_black
  ) +

  annotate(

    "label",

    x = x_right -
      0.03 * (x_right - x_left),

    y = 0.25,

    label = paste0(
      "Pseudobulk cutoff\n",
      "log2FC = 0.25"
    ),

    hjust = 1,
    vjust = -0.15,

    size = 3.0,

    label.size = 0.18,

    fill = "white",
    color = col_black
  ) +

  annotate(

    "label",

    x = x_left +
      0.03 * (x_right - x_left),

    y = y_high -
      0.05 * (y_high - y_low),

    label = criteria_text,

    hjust = 0,
    vjust = 1,

    size = 2.75,

    lineheight = 1.08,

    label.size = 0.20,

    fill = "white",

    color = col_black
  ) +

  scale_fill_manual(

    values = c(

      "Other candidates" =
        col_blue,

      "Passed all validation criteria" =
        col_rose
    )
  ) +

  coord_cartesian(

    xlim = c(
      x_left,
      x_right
    ),

    ylim = c(
      y_low,
      y_high
    ),

    clip = "off"
  ) +

  labs(

    x = paste0(
      "Single-cell differential expression, ",
      "log2FC"
    ),

    y = paste0(
      "Patient-paired pseudobulk, ",
      "log2FC"
    )
  ) +

  theme_pub(11.5) +

  theme(

    legend.position = c(
      0.74,
      0.13
    ),

    legend.background =
      element_blank(),

    legend.key.size =
      unit(0.75, "lines"),

    legend.text =
      element_text(size = 9.5)
  )

effect_order <- summary_df %>%

  arrange(PB_log2FC) %>%

  pull(gene)

effect_df <- summary_df %>%

  mutate(

    gene = factor(
      gene,
      levels = effect_order
    )
  )

effect_long <- bind_rows(

  summary_df %>%

    transmute(

      gene,

      Method =
        "Single-cell DE",

      Effect =
        FM_log2FC
    ),

  summary_df %>%

    transmute(

      gene,

      Method =
        "Patient-paired pseudobulk",

      Effect =
        PB_log2FC
    )

) %>%

  mutate(

    gene = factor(
      gene,
      levels = effect_order
    ),

    Method = factor(

      Method,

      levels = c(
        "Single-cell DE",
        "Patient-paired pseudobulk"
      )
    ),

    Validated =
      gene %in% validated_genes
  )

b_x_max <- max(
  effect_long$Effect,
  na.rm = TRUE
)

b_x_min <- min(
  effect_long$Effect,
  na.rm = TRUE
)

gene_labels <- setNames(
  effect_order,
  effect_order
)

pB <- ggplot(
  effect_df,
  aes(y = gene)
) +

  geom_vline(

    xintercept = 0.25,

    linetype = 2,

    linewidth = 0.60,

    color = col_grey3
  ) +

  geom_segment(

    aes(
      x = FM_log2FC,
      xend = PB_log2FC,
      yend = gene
    ),

    color = col_grey2,

    linewidth = 0.70
  ) +

  geom_point(

    data = effect_long %>%
      filter(!Validated),

    aes(
      x = Effect,
      fill = Method
    ),

    shape = 21,

    size = 4.0,

    color = "black",

    stroke = 0.38
  ) +

  geom_point(

    data = effect_long %>%
      filter(Validated),

    aes(
      x = Effect,
      fill = Method
    ),

    shape = 21,

    size = 4.8,

    color = col_rose,

    stroke = 1.05
  ) +

  annotate(

    "label",

    x = 0.25,

    y = length(effect_order) - 0.25,

    label = paste0(
      "Effect-size cutoff\n",
      "log2FC = 0.25"
    ),

    hjust = -0.06,

    vjust = 1,

    size = 3.2,

    fontface = "bold",

    label.size = 0.20,

    fill = "white",

    color = col_black
  ) +

  scale_fill_manual(

    values = c(

      "Single-cell DE" =
        col_purple,

      "Patient-paired pseudobulk" =
        col_teal
    )
  ) +

  scale_y_discrete(
    labels = gene_labels
  ) +

  scale_x_continuous(

    expand = expansion(
      mult = c(0.04, 0.05)
    )
  ) +

  labs(

    x = paste0(
      "Log2 fold-change ",
      "(C4 versus other clusters)"
    ),

    y = NULL
  ) +

  theme_pub(10.8) +

  theme(

    panel.grid.major.y =
      element_line(
        color = "#F0F0F0",
        linewidth = 0.40
      ),

    legend.position = "top",

    legend.text =
      element_text(size = 10),

    axis.text.y =
      element_text(
        face = "bold",
        size = 9.7
      )
  )

paired_plot <- paired_main %>%

  filter(
    gene %in% validated_genes
  ) %>%

  left_join(
    patient_map,
    by = "Patient"
  ) %>%

  mutate(

    gene = factor(
      gene,
      levels = validated_genes
    ),

    State = factor(

      State,

      levels = c(
        "Other",
        "C4"
      ),

      labels = c(
        "Other clusters",
        "C4"
      )
    ),

    DisplayID = factor(

      DisplayID,

      levels =
        unique(patient_map$DisplayID)
    )
  )

facet_info <- summary_df %>%

  filter(
    gene %in% validated_genes
  ) %>%

  mutate(

    facet_label = paste0(

      gene,

      "\nlog2FC = ",
      sprintf("%.2f", PB_log2FC),

      ", FDR ",
      fmt_fdr(PB_FDR)
    )
  )

facet_labels <- setNames(
  facet_info$facet_label,
  facet_info$gene
)

pC <- ggplot(

  paired_plot,

  aes(
    State,
    logCPM,
    group = DisplayID
  )

) +

  geom_line(

    color = col_grey2,

    linewidth = 0.65,

    alpha = 0.90
  ) +

  geom_point(

    aes(fill = State),

    shape = 21,

    size = 3.4,

    color = "black",

    stroke = 0.38
  ) +

  stat_summary(

    aes(group = State),

    fun = median,

    geom = "crossbar",

    width = 0.43,

    linewidth = 0.80,

    color = col_black
  ) +

  facet_wrap(

    ~gene,

    scales = "free_y",

    nrow = 1,

    labeller =
      as_labeller(facet_labels)
  ) +

  scale_fill_manual(

    values = c(

      "Other clusters" =
        col_blue,

      "C4" =
        col_rose
    )
  ) +

  labs(

    x = NULL,

    y =
      "Pseudobulk logCPM"
  ) +

  theme_pub(11.2) +

  theme(

    legend.position = "none",

    axis.text.x =
      element_text(
        face = "bold"
      ),

    strip.text =
      element_text(
        size = 10.8,
        face = "bold"
      )
  )

cluster_plot <- cluster_long %>%

  filter(
    gene %in% validated_genes
  ) %>%

  mutate(

    gene = factor(
      gene,
      levels = validated_genes
    ),

    Cluster =
      as.character(Cluster),

    Cluster_lab = factor(

      paste0(
        "C",
        Cluster
      ),

      levels = paste0(
        "C",
        sort(
          unique(
            as.integer(Cluster)
          )
        )
      )
    )
  )

cluster_cols <- c(

  "C0" = col_blue,
  "C1" = col_purple,
  "C2" = col_yellow,
  "C3" = col_orange,
  "C4" = col_rose
)

pD <- ggplot(

  cluster_plot,

  aes(
    Cluster_lab,
    logCPM
  )

) +

  geom_boxplot(

    aes(fill = Cluster_lab),

    width = 0.60,

    outlier.shape = NA,

    alpha = 0.90,

    linewidth = 0.50,

    color = "black"
  ) +

  geom_jitter(

    aes(fill = Cluster_lab),

    shape = 21,

    width = 0.10,

    size = 2.3,

    color = "black",

    stroke = 0.28,

    alpha = 0.92
  ) +

  facet_wrap(

    ~gene,

    scales = "free_y",

    nrow = 1
  ) +

  scale_fill_manual(
    values = cluster_cols
  ) +

  labs(

    x =
      "Validation cluster",

    y =
      "Pseudobulk logCPM"
  ) +

  theme_pub(11.2) +

  theme(

    legend.position = "none",

    strip.text =
      element_text(
        size = 10.9,
        face = "bold"
      )
  )

ggsave(

  file.path(
    out_dir,
    "12D_A_candidategene.pdf"
  ),

  pA,

  width = 9.2,
  height = 4.9,

  device = cairo_pdf
)

ggsave(

  file.path(
    out_dir,
    "12D_B_24gene.pdf"
  ),

  pB,

  width = 7.0,
  height = 10.9,

  device = cairo_pdf
)

ggsave(

  file.path(
    out_dir,
    "12D_C_finalgene.pdf"
  ),

  pC,

  width = 9.2,
  height = 4.2,

  device = cairo_pdf
)

ggsave(

  file.path(
    out_dir,
    "12D_D_finalgenecluster.pdf"
  ),

  pD,

  width = 9.2,
  height = 4.2,

  device = cairo_pdf
)

design <- c(

  area(
    t = 1,
    l = 2,
    b = 1,
    r = 2
  ),

  area(
    t = 1,
    l = 1,
    b = 3,
    r = 1
  ),

  area(
    t = 2,
    l = 2,
    b = 2,
    r = 2
  ),

  area(
    t = 3,
    l = 2,
    b = 3,
    r = 2
  )
)

main_fig <- pA + pB + pC + pD +

  plot_layout(

    design = design,

    widths = c(
      1.02,
      1.34
    ),

    heights = c(
      1.04,
      1.00,
      1.00
    )
  ) +

  plot_annotation(

    tag_levels = "A",

    theme = theme(

      plot.tag = element_text(

        face = "bold",

        size = 16,

        color = col_black
      )
    )
  )

pdf_file <- file.path(
  out_dir,
  "12D_plot_candidategene_analysis version.pdf"
)

png_file <- file.path(
  out_dir,
  "12D_plot_candidategene_analysis version.png"
)

ggsave(

  filename = pdf_file,

  plot = main_fig,

  width = 16.0,

  height = 11.5,

  device = cairo_pdf
)

ggsave(

  filename = png_file,

  plot = main_fig,

  width = 16.0,

  height = 11.5,

  dpi = 600,

  bg = "white"
)

caption <- c(

  "Panel A",
  paste0(
    "The shaded region indicates concordant effect sizes ",
    "(log2FC >= 0.25) in both single-cell differential ",
    "expression and patient-paired pseudobulk analyses."
  ),
  paste0(
    "Genes highlighted in rose passed all validation criteria: ",
    "FDR < 0.05 and log2FC >= 0.25 in both analyses, ",
    "positive direction in >=75% of C4-versus-cluster contrasts, ",
    "and higher C4 expression in >=70% of evaluable patients."
  ),

  "",

  "Panel B",
  paste0(
    "Paired effect sizes for all candidate genes. ",
    "The dashed vertical line indicates the predefined ",
    "effect-size cutoff of log2FC = 0.25."
  ),

  "",

  "Panel C",
  paste0(
    "Patient-paired pseudobulk expression of genes passing ",
    "all validation criteria in C4 versus other malignant ",
    "epithelial clusters."
  ),

  "",

  "Panel D",
  paste0(
    "Patient-by-cluster pseudobulk expression of validated ",
    "genes across C0-C4."
  )
)

writeLines(

  caption,

  file.path(
    out_dir,
    "12D_plot.txt"
  )
)
