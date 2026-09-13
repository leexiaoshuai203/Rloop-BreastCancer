# scTenifoldKnk virtual knockout in GSE306201 C5 malignant epithelial cells

#### HMGB2 / H2AZ2 / UBALD2
#### Discovery malignant epithelial C5

rm(list=ls()); options(stringsAsFactors=FALSE); set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
  library(scTenifoldKnk); library(scTenifoldNet)
})

project_dir <- "."

input_file <- file.path(
  project_dir, "results", "GSE306201", "Malignant",
  "Malignant_Epithelial_Reclustered.rds"
)

out_dir <- file.path(project_dir, "results", "scTenifoldKnk", "C5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rl65_file <- file.path(project_dir, "gene_sets", "RL_Sig65.txt")
candidate_file <- file.path(project_dir, "gene_sets", "candidate24_genes.txt")

need <- c(input_file,rl65_file,candidate_file)
if(any(!file.exists(need))) stop("Missing required files:\n", paste(need[!file.exists(need)], collapse="\n"))

KO_GENES <- c("HMGB2","H2AZ2","UBALD2")
SEED <- 1234L; HVG_N <- 2000L; N_NET <- 20L
N_CORES <- max(1L,min(4L,parallel::detectCores(logical=FALSE)))

canon <- function(x){
  x <- toupper(trimws(as.character(x)))
  x[x=="H2AFV"] <- "H2AZ2"
  x
}
read_genes <- function(f){
  x <- canon(readLines(f,warn=FALSE))
  unique(x[!is.na(x) & nzchar(x)])
}

obj <- readRDS(input_file)
DefaultAssay(obj) <- "RNA"

if(length(Layers(obj[["RNA"]]))>1) obj <- JoinLayers(obj,assay="RNA")
if(!"seurat_clusters" %in% colnames(obj[[]])) stop("seurat_clusters is required.")
if(!"5" %in% as.character(unique(obj$seurat_clusters)))
  stop("Cluster 5 was not found. Available clusters: ",paste(sort(unique(obj$seurat_clusters)),collapse=", "))

c5 <- subset(obj,subset=seurat_clusters=="5")
DefaultAssay(c5) <- "RNA"

if(ncol(c5)<100) warning("C5 contains fewer than 100 cells; virtual-knockout stability may be limited.")
if(ncol(c5)<60) stop("C5 contains fewer than 60 cells; network analysis was not run.")

batch_col <- if("orig.ident" %in% colnames(c5[[]])) "orig.ident" else NULL
if(!is.null(batch_col)){
  sample_tab <- sort(table(c5[[batch_col,drop=TRUE]]),decreasing=TRUE)
print(sample_tab)
} else sample_tab <- NULL

MIN_GENE_CELLS <- min(25L,max(10L,ceiling(ncol(c5)*0.10)))

c5 <- NormalizeData(c5,verbose=FALSE)
c5 <- FindVariableFeatures(c5,selection.method="vst",nfeatures=HVG_N,verbose=FALSE)

counts <- GetAssayData(c5,assay="RNA",layer="counts")
hvg <- canon(VariableFeatures(c5))

new_names <- canon(rownames(counts))
if(anyDuplicated(new_names)){
  dup <- unique(new_names[duplicated(new_names)])
  stop("Duplicated gene symbols after canonicalization: ",paste(dup,collapse=", "))
}
rownames(counts) <- new_names

rl65 <- read_genes(rl65_file)
candidate24 <- read_genes(candidate_file)

n_exp <- Matrix::rowSums(counts>0); names(n_exp) <- rownames(counts)
eligible <- names(n_exp)[n_exp>=MIN_GENE_CELLS &
                           !grepl("^MT-",names(n_exp),ignore.case=TRUE)]

ko_check <- data.table(
  Gene=KO_GENES,
  Expressed_cells=sapply(KO_GENES,function(g) if(g %in% names(n_exp)) n_exp[g] else 0)
)
print(ko_check)

if(any(ko_check$Expressed_cells<10)){
  bad <- ko_check[Expressed_cells<10,Gene]
  stop("KO genes expressed in fewer than 10 C5 cells: ",paste(bad,collapse=", "))
}

mandatory <- unique(c(KO_GENES,intersect(rl65,eligible),intersect(candidate24,eligible)))
network_genes <- unique(c(intersect(hvg,eligible),mandatory))
network_genes <- intersect(network_genes,rownames(counts))

counts_net <- counts[network_genes,,drop=FALSE]
expr <- scTenifoldNet::cpmNormalization(counts_net)
expr <- as.matrix(expr)

gene_sd <- apply(expr,1,sd)
expr <- expr[is.finite(gene_sd) & gene_sd>0,,drop=FALSE]

if(any(!KO_GENES %in% rownames(expr)))
  stop("KO genes missing from the final network: ", paste(setdiff(KO_GENES, rownames(expr)), collapse=", "))

N_SUBCELL <- min(200L,floor(ncol(expr)*0.80))
if(N_SUBCELL>=ncol(expr)) N_SUBCELL <- ncol(expr)-1L
if(N_SUBCELL<50) stop("Fewer than 50 cells are available per network replicate.")
if(N_SUBCELL<100) warning("Fewer than 100 cells are available per network replicate; stability may be limited.")

coverage <- data.table(
  GeneSet=c("RL-Sig65","24 candidates"),
  Input_N=c(length(rl65),length(candidate24)),
  Network_N=c(length(intersect(rl65,rownames(expr))),
              length(intersect(candidate24,rownames(expr))))
)
coverage[,Coverage:=Network_N/Input_N]

qc <- rbindlist(list(
  data.table(Section="Dataset",
             Item=c("All_malignant_cells","C5_cells","C5_samples",
                    "Network_genes","Networks","Cells_per_network"),
             Value=as.character(c(ncol(obj),ncol(c5),
                                  ifelse(is.null(sample_tab),NA,length(sample_tab)),
                                  nrow(expr),N_NET,N_SUBCELL))),
  data.table(Section="KO_expression",Item=ko_check$Gene,
             Value=as.character(ko_check$Expressed_cells)),
  data.table(Section="GeneSet_coverage",Item=coverage$GeneSet,
             Value=paste0(coverage$Network_N,"/",coverage$Input_N,
                          " (",sprintf("%.1f",100*coverage$Coverage),"%)"))
))
fwrite(qc,file.path(out_dir,"02_C5_input_network_QC.csv"))

fits <- setNames(vector("list",length(KO_GENES)),KO_GENES)
dr_list <- vector("list",length(KO_GENES))

for(i in seq_along(KO_GENES)){
  ko_gene <- KO_GENES[i]

  set.seed(SEED)

  fit <- scTenifoldKnk::scTenifoldKnk(
    countMatrix=expr,gKO=ko_gene,qc=FALSE,
    nc_nNet=N_NET,nc_nCells=N_SUBCELL,nc_nComp=3,
    nc_scaleScores=TRUE,nc_symmetric=FALSE,nc_lambda=0,nc_q=0.90,
    td_K=3,td_maxIter=1000,td_maxError=1e-05,td_nDecimal=3,
    ma_nDim=2,nCores=N_CORES
  )

  if(!all(c("tensorNetworks","manifoldAlignment","diffRegulation") %in% names(fit)))
    stop(ko_gene, ": unexpected scTenifoldKnk output structure.")

  d <- as.data.table(fit$diffRegulation)
  need_cols <- c("gene","distance","Z","FC","p.value","p.adj")
  if(!all(need_cols %in% names(d)))
    stop(ko_gene, ": diffRegulation is missing columns: ", paste(setdiff(need_cols,names(d)), collapse=", "))

  d[,gene:=canon(gene)]
  d[,KO:=ko_gene]
  d[,Self_KO:=gene==ko_gene]
  d[,FDR05:=!is.na(p.adj) & p.adj<0.05]
  d[,FDR10:=!is.na(p.adj) & p.adj<0.10]
  setcolorder(d,c("KO","gene","Self_KO","distance","Z","FC",
                  "p.value","p.adj","FDR05","FDR10"))

  fits[[ko_gene]] <- fit
  dr_list[[i]] <- d

}

all_dr <- rbindlist(dr_list,use.names=TRUE,fill=TRUE)
expected_n <- length(KO_GENES)*nrow(expr)
if(nrow(all_dr)!=expected_n)
  stop("Unexpected merged result size: ", nrow(all_dr), "; expected: ", expected_n)

saveRDS(fits,file.path(out_dir,"01_scTenifoldKnk_three_gene_raw_results.rds"))
fwrite(all_dr,file.path(out_dir,"03_three_gene_differential_regulation.csv"))

max_diff <- function(A,B){
  D <- A-B
  if(length(D@x)==0) return(0)
  max(abs(D@x),na.rm=TRUE)
}

wt_check <- data.table(
  Comparison=c("HMGB2_vs_H2AZ2","HMGB2_vs_UBALD2","H2AZ2_vs_UBALD2"),
  Max_abs_difference=c(
    max_diff(fits$HMGB2$tensorNetworks$WT,fits$H2AZ2$tensorNetworks$WT),
    max_diff(fits$HMGB2$tensorNetworks$WT,fits$UBALD2$tensorNetworks$WT),
    max_diff(fits$H2AZ2$tensorNetworks$WT,fits$UBALD2$tensorNetworks$WT)
  )
)
print(wt_check)
fwrite(wt_check,file.path(out_dir,"04_WT_network_consistency_check.csv"))

ko_summary <- all_dr[Self_KO==FALSE,.(
  Genes_tested=.N,
  FDR05=sum(FDR05,na.rm=TRUE),
  FDR10=sum(FDR10,na.rm=TRUE),
  Median_Z=median(Z,na.rm=TRUE),
  Mean_Z=mean(Z,na.rm=TRUE),
  Max_Z=max(Z,na.rm=TRUE)
),by=KO]

print(ko_summary)
fwrite(ko_summary,file.path(out_dir,"05_three_gene_virtual_knockout_summary.csv"))
