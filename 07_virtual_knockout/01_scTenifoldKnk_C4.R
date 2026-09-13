# scTenifoldKnk virtual knockout in GSE245601 C4 malignant epithelial cells

#### HMGB2 / H2AZ2 / UBALD2

rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(1234)

if(!requireNamespace("scTenifoldKnk",quietly=TRUE))
  install.packages("scTenifoldKnk")
if(!requireNamespace("scTenifoldNet",quietly=TRUE))
  install.packages("scTenifoldNet")
if(!requireNamespace("patchwork",quietly=TRUE))
  install.packages("patchwork")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(scTenifoldKnk)
  library(scTenifoldNet)
  library(ggplot2)
  library(patchwork)
})

project_dir <- "."

input_file <- file.path(
  project_dir, "results", "GSE245601", "10.malignat",
  "GSE245601_Malignant_Epithelial_Reclustered.rds"
)

out_dir <- file.path(project_dir, "results", "scTenifoldKnk", "C4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rl65_file <- file.path(project_dir, "gene_sets", "RL_Sig65.txt")
c5_file <- file.path(project_dir, "gene_sets", "C5_stable_up_genes.txt")
candidate24_file <- file.path(project_dir, "gene_sets", "candidate24_genes.txt")

need <- c(input_file,rl65_file,c5_file,candidate24_file)
if(any(!file.exists(need))){
  stop("Missing required files:\n", paste(need[!file.exists(need)], collapse="\n"))
}

KO_GENES <- c("HMGB2","H2AZ2","UBALD2")

SEED <- 1234L
HVG_N <- 2000L
MIN_GENE_CELLS <- 25L
MIN_KO_CELLS <- 10L

N_NET <- 20L
N_SUBCELL_MAX <- 200L
N_CORES <- max(1L,min(4L,parallel::detectCores(logical=FALSE)))

canon <- function(x){
  x <- toupper(trimws(as.character(x)))
  x[x=="H2AFV"] <- "H2AZ2"
  x
}

read_gene_set <- function(file){
  x <- canon(readLines(file,warn=FALSE))
  unique(x[!is.na(x) & nzchar(x)])
}

fmt_fdr <- function(x){
  out <- rep("NA",length(x))
  ok <- !is.na(x)

  out[ok & x<1e-99] <- "<1e-99"
  ii <- ok & x>=1e-99 & x<0.001
  out[ii] <- formatC(x[ii],format="e",digits=1)
  ii <- ok & x>=0.001
  out[ii] <- sprintf("%.3f",x[ii])

  out
}

obj <- readRDS(input_file)
DefaultAssay(obj) <- "RNA"

if(length(Layers(obj[["RNA"]]))>1){
  obj <- JoinLayers(obj,assay="RNA")
}

if(!"C_cluster" %in% colnames(obj[[]])){
  if(!"seurat_clusters" %in% colnames(obj[[]]))
    stop("C_cluster or seurat_clusters is required.")

  obj$C_cluster <- paste0("C",obj$seurat_clusters)
}

if(!"C4" %in% as.character(unique(obj$C_cluster))){
  stop("C4 was not found. Available clusters: ",
       paste(sort(unique(as.character(obj$C_cluster))),collapse=", "))
}

c4 <- subset(obj,subset=C_cluster=="C4")
DefaultAssay(c4) <- "RNA"

if(ncol(c4)<150)
  stop("C4 contains fewer than 150 cells; network analysis was not run.")

batch_col <- if("sample_gsm" %in% colnames(c4[[]])){
  "sample_gsm"
} else {
  "orig.ident"
}

sample_tab <- sort(table(c4[[batch_col,drop=TRUE]]),decreasing=TRUE)

print(sample_tab)

c4 <- NormalizeData(c4,verbose=FALSE)
c4 <- FindVariableFeatures(
  c4,
  selection.method="vst",
  nfeatures=HVG_N,
  verbose=FALSE
)

hvg <- canon(VariableFeatures(c4))

counts <- GetAssayData(
  c4,
  assay="RNA",
  layer="counts"
)

new_names <- canon(rownames(counts))

if(anyDuplicated(new_names)){
  dup <- unique(new_names[duplicated(new_names)])
  stop(
    "Duplicated gene symbols after canonicalization: ",
    paste(dup,collapse=", "),
    "\nCheck H2AFV/H2AZ2 in the input object."
  )
}

rownames(counts) <- new_names

rl65 <- read_gene_set(rl65_file)
c5_sig <- read_gene_set(c5_file)
candidate24 <- read_gene_set(candidate24_file)

gene_sets <- list(
  `RL-Sig65` = rl65,
  `C5 stable-up` = c5_sig,
  `24 candidates` = candidate24
)

n_exp <- Matrix::rowSums(counts>0)
names(n_exp) <- rownames(counts)

eligible <- names(n_exp)[
  n_exp>=MIN_GENE_CELLS &
    !grepl("^MT-",names(n_exp),ignore.case=TRUE)
]

ko_check <- data.table(
  Gene=KO_GENES,
  Expressed_cells=sapply(
    KO_GENES,
    function(g) if(g %in% names(n_exp)) n_exp[g] else 0
  )
)

print(ko_check)

if(any(ko_check$Expressed_cells<MIN_KO_CELLS)){
  bad <- ko_check[Expressed_cells<MIN_KO_CELLS,Gene]
  stop(
    "KO genes with insufficient C4 expression: ",
    paste(bad,collapse=", ")
  )
}

mandatory <- unique(c(
  KO_GENES,
  intersect(rl65,eligible),
  intersect(candidate24,eligible)
))

network_genes <- unique(c(
  intersect(hvg,eligible),
  mandatory
))

network_genes <- intersect(network_genes,rownames(counts))

counts_net <- counts[network_genes,,drop=FALSE]

expr <- scTenifoldNet::cpmNormalization(counts_net)
expr <- as.matrix(expr)

gene_sd <- apply(expr,1,sd)
keep <- is.finite(gene_sd) & gene_sd>0
expr <- expr[keep,,drop=FALSE]

if(any(!KO_GENES %in% rownames(expr))){
  stop(
    "KO genes missing from the final network: ",
    paste(setdiff(KO_GENES,rownames(expr)),collapse=", ")
  )
}

N_SUBCELL <- min(
  N_SUBCELL_MAX,
  floor(ncol(expr)*0.80)
)

if(N_SUBCELL>=ncol(expr))
  N_SUBCELL <- ncol(expr)-1L

if(N_SUBCELL<100)
  stop("Fewer than 100 cells are available per network replicate.")

coverage <- rbindlist(lapply(names(gene_sets),function(nm){
  gs <- gene_sets[[nm]]

  data.table(
    GeneSet=nm,
    Input_N=length(gs),
    Network_N=length(intersect(gs,rownames(expr))),
    Coverage=length(intersect(gs,rownames(expr)))/length(gs)
  )
}))

print(coverage)

if(any(coverage$Network_N<5))
  stop("At least one core gene set has fewer than five genes in the network.")

qc_table <- rbindlist(list(
  data.table(
    Section="Dataset",
    Item=c(
      "All_malignant_cells","C4_cells","C4_samples",
      "Network_genes","Networks","Cells_per_network"
    ),
    Value=as.character(c(
      ncol(obj),ncol(c4),length(sample_tab),
      nrow(expr),N_NET,N_SUBCELL
    ))
  ),

  data.table(
    Section="Sample",
    Item=names(sample_tab),
    Value=as.character(as.numeric(sample_tab))
  ),

  data.table(
    Section="KO_expression",
    Item=ko_check$Gene,
    Value=as.character(ko_check$Expressed_cells)
  ),

  data.table(
    Section="GeneSet_coverage",
    Item=coverage$GeneSet,
    Value=paste0(
      coverage$Network_N,"/",
      coverage$Input_N," (",
      sprintf("%.1f",100*coverage$Coverage),"%)"
    )
  ),

  data.table(
    Section="Software",
    Item=c("scTenifoldKnk","Seurat"),
    Value=c(
      as.character(packageVersion("scTenifoldKnk")),
      as.character(packageVersion("Seurat"))
    )
  )
))

fwrite(
  qc_table,
  file.path(out_dir,"02_C4_input_network_QC.csv")
)

fits <- setNames(vector("list",length(KO_GENES)),KO_GENES)
dr_list <- vector("list",length(KO_GENES))

for(i in seq_along(KO_GENES)){

  ko_gene <- KO_GENES[i]

  set.seed(SEED)

  fit <- scTenifoldKnk::scTenifoldKnk(
    countMatrix=expr,
    gKO=ko_gene,
    qc=FALSE,
    nc_nNet=N_NET,
    nc_nCells=N_SUBCELL,
    nc_nComp=3,
    nc_scaleScores=TRUE,
    nc_symmetric=FALSE,
    nc_lambda=0,
    nc_q=0.90,
    td_K=3,
    td_maxIter=1000,
    td_maxError=1e-05,
    td_nDecimal=3,
    ma_nDim=2,
    nCores=N_CORES
  )

  if(!all(c("tensorNetworks","manifoldAlignment",
            "diffRegulation") %in% names(fit))){
    stop(ko_gene, ": unexpected scTenifoldKnk output structure.")
  }

  d <- as.data.table(fit$diffRegulation)

  required_cols <- c(
    "gene","distance","Z","FC","p.value","p.adj"
  )

  if(!all(required_cols %in% names(d))){
    stop(
      ko_gene, ": diffRegulation is missing columns: ",
      paste(setdiff(required_cols,names(d)),collapse=", ")
    )
  }

  d[,gene:=canon(gene)]

  if(anyDuplicated(d$gene))
    stop(ko_gene, ": duplicated genes were found in diffRegulation.")

  d[,KO:=ko_gene]
  d[,Self_KO:=gene==ko_gene]
  d[,FDR05:=!is.na(p.adj) & p.adj<0.05]
  d[,FDR10:=!is.na(p.adj) & p.adj<0.10]

  setcolorder(
    d,
    c(
      "KO","gene","Self_KO",
      "distance","Z","FC",
      "p.value","p.adj",
      "FDR05","FDR10"
    )
  )

  fits[[ko_gene]] <- fit
  dr_list[[i]] <- d

}

all_dr <- rbindlist(dr_list,use.names=TRUE,fill=TRUE)

expected_n <- length(KO_GENES)*nrow(expr)

if(nrow(all_dr)!=expected_n){
  stop(
    "Unexpected merged result size: ", nrow(all_dr),
    "; expected: ", expected_n
  )
}

saveRDS(
  fits,
  file.path(out_dir,"01_scTenifoldKnk_three_gene_raw_results.rds")
)

fwrite(
  all_dr,
  file.path(out_dir,"03_three_gene_differential_regulation.csv")
)

max_sparse_diff <- function(A,B){
  D <- A-B
  xx <- D@x
  if(length(xx)==0) return(0)
  max(abs(xx),na.rm=TRUE)
}

wt_check <- data.table(
  Comparison=c(
    "HMGB2_vs_H2AZ2",
    "HMGB2_vs_UBALD2",
    "H2AZ2_vs_UBALD2"
  ),
  Max_abs_difference=c(
    max_sparse_diff(
      fits$HMGB2$tensorNetworks$WT,
      fits$H2AZ2$tensorNetworks$WT
    ),
    max_sparse_diff(
      fits$HMGB2$tensorNetworks$WT,
      fits$UBALD2$tensorNetworks$WT
    ),
    max_sparse_diff(
      fits$H2AZ2$tensorNetworks$WT,
      fits$UBALD2$tensorNetworks$WT
    )
  )
)

print(wt_check)

if(any(wt_check$Max_abs_difference>1e-6)){
  warning(
    "WT networks differ numerically across KO runs; inspect stability analyses."
  )
}

program_test <- function(ko_gene,set_name,set_genes){

  d <- all_dr[
    KO==ko_gene &
      Self_KO==FALSE &
      is.finite(Z)
  ]

  gs <- intersect(
    unique(canon(set_genes)),
    d$gene
  )

  if(length(gs)<5){
    return(data.table(
      KO=ko_gene,
      GeneSet=set_name,
      N_in_network=length(gs),
      Median_Z=NA_real_,
      Background_median_Z=NA_real_,
      Delta_median_Z=NA_real_,
      Wilcox_P=NA_real_,
      Top10_N=NA_integer_,
      Top10_overlap=NA_integer_,
      Top10_OR=NA_real_,
      Fisher_P=NA_real_
    ))
  }

  z_set <- d[gene%in%gs,Z]
  z_bg <- d[!gene%in%gs,Z]

  wt <- wilcox.test(
    z_set,z_bg,
    alternative="greater",
    exact=FALSE
  )

  top_cut <- quantile(
    d$Z,
    probs=0.90,
    na.rm=TRUE,
    names=FALSE
  )

  d[,Top10:=Z>=top_cut]

  a <- d[gene%in%gs & Top10==TRUE,.N]
  b <- d[gene%in%gs & Top10==FALSE,.N]
  c <- d[!gene%in%gs & Top10==TRUE,.N]
  dd <- d[!gene%in%gs & Top10==FALSE,.N]

  ft <- fisher.test(
    matrix(c(a,b,c,dd),nrow=2,byrow=TRUE),
    alternative="greater"
  )

  data.table(
    KO=ko_gene,
    GeneSet=set_name,
    N_in_network=length(gs),
    Median_Z=median(z_set,na.rm=TRUE),
    Background_median_Z=median(z_bg,na.rm=TRUE),
    Delta_median_Z=median(z_set,na.rm=TRUE)-median(z_bg,na.rm=TRUE),
    Wilcox_P=wt$p.value,
    Top10_N=sum(d$Top10),
    Top10_overlap=a,
    Top10_OR=unname(ft$estimate),
    Fisher_P=ft$p.value
  )
}

program_stats <- rbindlist(
  lapply(KO_GENES,function(ko_gene){
    rbindlist(
      lapply(names(gene_sets),function(nm){
        program_test(
          ko_gene,
          nm,
          gene_sets[[nm]]
        )
      })
    )
  })
)

program_stats[,Wilcox_FDR:=p.adjust(Wilcox_P,method="BH")]
program_stats[,Fisher_FDR:=p.adjust(Fisher_P,method="BH")]

fwrite(
  program_stats,
  file.path(out_dir,"04_program_perturbation_statistics.csv")
)

print(program_stats)

WT <- fits[[1]]$tensorNetworks$WT

if(!all(KO_GENES %in% rownames(WT)))
  stop("Candidate genes are missing from the WT network.")

edge_top <- rbindlist(
  lapply(KO_GENES,function(ko_gene){

    ww <- as.numeric(WT[ko_gene,])
    names(ww) <- colnames(WT)

    dd <- data.table(
      KO=ko_gene,
      Target=canon(names(ww)),
      Weight=ww
    )

    dd <- dd[
      Target!=ko_gene &
        is.finite(Weight) &
        Weight!=0
    ]

    dd[,AbsWeight:=abs(Weight)]
    setorder(dd,-AbsWeight)

    head(dd,10)
  })
)

neighbor_genes <- unique(edge_top$Target)

edge_all <- rbindlist(
  lapply(KO_GENES,function(ko_gene){

    ww <- as.numeric(WT[ko_gene,])
    names(ww) <- canon(colnames(WT))

    data.table(
      KO=ko_gene,
      Target=neighbor_genes,
      Weight=ww[neighbor_genes]
    )
  })
)

edge_all[is.na(Weight),Weight:=0]
edge_all[,AbsWeight:=abs(Weight)]
edge_all[,RL_Sig65:=Target%in%rl65]
edge_all[,C5_associated:=Target%in%c5_sig]
edge_all[,Candidate24:=Target%in%candidate24]

fwrite(
  edge_all,
  file.path(out_dir,"05_WT_direct_network_neighborhood.csv")
)
