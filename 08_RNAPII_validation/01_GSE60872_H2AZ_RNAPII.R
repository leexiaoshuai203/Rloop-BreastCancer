# GSE60872 H2A.Z-RNAPII analysis

PROJECT_DIR <- "."

options(stringsAsFactors=FALSE)
set.seed(1234)

if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")

cran <- c("data.table","ggplot2","patchwork","scales")
bio  <- c("AnnotationDbi","org.Hs.eg.db")

for(p in cran){
  if(!requireNamespace(p,quietly=TRUE))
    install.packages(p,dependencies=TRUE)
}
for(p in bio){
  if(!requireNamespace(p,quietly=TRUE))
    BiocManager::install(p,ask=FALSE,update=FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

project <- PROJECT_DIR

in_dir <- file.path(PROJECT_DIR, "data", "raw", "GSE60872")

out_dir <- file.path(PROJECT_DIR, "results", "GSE60872")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

rl65_file <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")

c5_file <- file.path(PROJECT_DIR, "gene_sets/C5_stable_up_genes.txt")

cand24_file <- file.path(PROJECT_DIR, "gene_sets/candidate24_genes.txt")

drip_file <- file.path(PROJECT_DIR, "results", "GSE81851", "enrichment", "03_gene_level_E2_DRIP_effects.csv")

need <- c(rl65_file,c5_file,cand24_file,drip_file)
if(any(!file.exists(need))){
  stop("missing: \n",paste(need[!file.exists(need)],collapse="\n"))
}

csvs <- list.files(
  in_dir,
  pattern="\\.csv$",
  recursive=TRUE,
  full.names=TRUE,
  ignore.case=TRUE
)

pick_file <- function(key){
  z <- csvs[grepl(key,basename(csvs),fixed=TRUE)]
  if(length(z)!=1){
    stop(
      ": ",key,
      "\ncurrent: \n",paste(z,collapse="\n")
    )
  }
  z
}

ctrl_file <- pick_file("PolII-control-density")
kd_file   <- pick_file("PolII-siH2AZ-density")
input_file<- pick_file("MCF7-input-density")

canon <- function(x){
  x <- toupper(trimws(as.character(x)))
  x[x=="H2AFV"] <- "H2AZ2"
  x
}

read_genes <- function(f){
  x <- canon(readLines(f,warn=FALSE))
  unique(x[!is.na(x) & nzchar(x)])
}

read_den <- function(f,suffix){

  d <- fread(f)

  if(ncol(d)!=3){
    stop("CSV3: ",basename(f))
  }

  setnames(d,1,"REFSEQ")

  if(!all(c("TSS","Gene.Body") %in% names(d))){
    stop("missingTSSGene.Body: ",basename(f))
  }

  setnames(
    d,
    c("TSS","Gene.Body"),
    c(
      paste0("TSS_",suffix),
      paste0("GB_",suffix)
    )
  )

  d[,REFSEQ:=sub("\\..*$","",as.character(REFSEQ))]
  d
}

safe_log2ratio <- function(a,b){
  z <- rep(NA_real_,length(a))
  ok <- is.finite(a) & is.finite(b) & a>0 & b>0
  z[ok] <- log2(a[ok]/b[ok])
  z
}

savepdf <- function(name,p,w,h){
  ggsave(
    file.path(out_dir,name),p,
    width=w,height=h,units="cm",
    device=cairo_pdf
  )
}

ctrl <- read_den(ctrl_file,"ctrl")
kd   <- read_den(kd_file,"kd")
inp  <- read_den(input_file,"input")

x <- Reduce(
  function(a,b) merge(a,b,by="REFSEQ",all=FALSE),
  list(ctrl,kd,inp)
)

sym <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys=x$REFSEQ,
  column="SYMBOL",
  keytype="REFSEQ",
  multiVals="first"
)

x[,SYMBOL:=canon(unname(sym))]

mapping_stat <- data.table(
  N_RefSeq=nrow(x),
  N_mapped=sum(!is.na(x$SYMBOL) & nzchar(x$SYMBOL)),
  Mapping_fraction=mean(!is.na(x$SYMBOL) & nzchar(x$SYMBOL))
)

fwrite(
  mapping_stat,
  file.path(out_dir,"00_H2AZ_RefSeqstatistics.csv")
)

x[,`:=`(
  TSS_ctrl_bg=TSS_ctrl-TSS_input,
  TSS_kd_bg=TSS_kd-TSS_input,
  GB_ctrl_bg=GB_ctrl-GB_input,
  GB_kd_bg=GB_kd-GB_input
)]

x[,Delta_TSS_bg:=safe_log2ratio(TSS_kd_bg,TSS_ctrl_bg)]
x[,Delta_GB_bg :=safe_log2ratio(GB_kd_bg,GB_ctrl_bg)]

x[,PI_ctrl_bg:=ifelse(
  TSS_ctrl_bg>0 & GB_ctrl_bg>0,
  TSS_ctrl_bg/GB_ctrl_bg,
  NA_real_
)]

x[,PI_kd_bg:=ifelse(
  TSS_kd_bg>0 & GB_kd_bg>0,
  TSS_kd_bg/GB_kd_bg,
  NA_real_
)]

x[,Delta_PI_bg:=safe_log2ratio(PI_kd_bg,PI_ctrl_bg)]

x[,Delta_TSS_norm:=safe_log2ratio(TSS_kd,TSS_ctrl)]
x[,Delta_GB_norm :=safe_log2ratio(GB_kd,GB_ctrl)]

x[,PI_ctrl_norm:=ifelse(
  TSS_ctrl>0 & GB_ctrl>0,
  TSS_ctrl/GB_ctrl,
  NA_real_
)]

x[,PI_kd_norm:=ifelse(
  TSS_kd>0 & GB_kd>0,
  TSS_kd/GB_kd,
  NA_real_
)]

x[,Delta_PI_norm:=safe_log2ratio(
  PI_kd_norm,
  PI_ctrl_norm
)]

tmp <- x[
  !is.na(SYMBOL) &
    nzchar(SYMBOL)
]

tmp[,rank_TSS:=ifelse(
  is.finite(TSS_ctrl_bg),
  TSS_ctrl_bg,
  -Inf
)]

setorder(
  tmp,
  SYMBOL,
  -rank_TSS,
  -TSS_ctrl
)

gene <- tmp[!duplicated(SYMBOL)]
gene[,rank_TSS:=NULL]

gene[,PI_eligible :=
       is.finite(TSS_ctrl_bg) &
       TSS_ctrl_bg>=0.001 &
       is.finite(PI_ctrl_bg) &
       is.finite(PI_kd_bg)]

gene[PI_eligible==FALSE,Delta_PI_bg:=NA_real_]

global_test <- function(d,ctrl_col,kd_col,label){

  d <- d[
    is.finite(get(ctrl_col)) &
      is.finite(get(kd_col)) &
      get(ctrl_col)>0 &
      get(kd_col)>0
  ]

  a <- d[[ctrl_col]]
  b <- d[[kd_col]]

  p_unpaired <- wilcox.test(
    b,a,
    paired=FALSE,
    exact=FALSE
  )$p.value

  p_paired <- wilcox.test(
    b,a,
    paired=TRUE,
    exact=FALSE
  )$p.value

  data.table(
    Method=label,
    N=nrow(d),
    Median_control=median(a),
    Median_KD=median(b),
    Median_log2_change=median(log2(b/a)),
    P_MannWhitney=p_unpaired,
    P_paired=p_paired
  )
}

g1 <- global_test(
  gene[PI_eligible==TRUE],
  "PI_ctrl_bg","PI_kd_bg",
  "Input-corrected PI"
)

g2 <- global_test(
  gene,
  "PI_ctrl_norm","PI_kd_norm",
  "Normalized-density ratio sensitivity"
)

global <- rbind(g1,g2)
fwrite(
  global,
  file.path(out_dir,"01_H2AZ_PI.csv")
)

print(global)

if(g1$Median_KD > g1$Median_control){
}else{
}

rl65   <- read_genes(rl65_file)
c5     <- read_genes(c5_file)
cand24 <- read_genes(cand24_file)

sets <- list(
  "RL-Sig65"=rl65,
  "C5 stable-up"=c5,
  "24 candidates"=cand24
)

metric_map <- c(
  "TSS RNAPII"="Delta_TSS_bg",
  "Gene-body RNAPII"="Delta_GB_bg",
  "Pausing index"="Delta_PI_bg"
)

test_set <- function(set_name,genes,metric_name,col){

  d <- gene[
    is.finite(get(col)),
    .(SYMBOL,Score=get(col))
  ]

  g <- intersect(genes,d$SYMBOL)

  a <- d[SYMBOL %chin% g,Score]
  b <- d[!SYMBOL %chin% g,Score]

  if(length(a)<3 || length(b)<3){
    p <- NA_real_
  }else{
    p <- wilcox.test(
      a,b,
      alternative="two.sided",
      exact=FALSE
    )$p.value
  }

  data.table(
    GeneSet=set_name,
    Metric=metric_name,
    N_set=length(a),
    N_background=length(b),
    Median_set=median(a,na.rm=TRUE),
    Median_background=median(b,na.rm=TRUE),
    Delta_median=median(a,na.rm=TRUE)-median(b,na.rm=TRUE),
    P=p
  )
}

set_stat <- rbindlist(lapply(names(sets),function(s){
  rbindlist(lapply(names(metric_map),function(m){
    test_set(
      s,
      sets[[s]],
      m,
      unname(metric_map[m])
    )
  }))
}))

set_stat[,FDR:=p.adjust(P,"BH")]

fwrite(
  set_stat,
  file.path(out_dir,"02_H2AZ_gene_setPolIIstatistics.csv")
)

sens_map <- c(
  "TSS RNAPII"="Delta_TSS_norm",
  "Gene-body RNAPII"="Delta_GB_norm",
  "Density-ratio PI proxy"="Delta_PI_norm"
)

sens_stat <- rbindlist(lapply(names(sets),function(s){
  rbindlist(lapply(names(sens_map),function(m){
    test_set(
      s,
      sets[[s]],
      m,
      unname(sens_map[m])
    )
  }))
}))

sens_stat[,FDR:=p.adjust(P,"BH")]

fwrite(
  sens_stat,
  file.path(out_dir,"03_H2AZ_normalized_densityanalysis.csv")
)

dr <- fread(drip_file)

req_dr <- c(
  "SYMBOL",
  "Mean_LFC_2h",
  "Mean_LFC_24h",
  "Strict_E2_DRIP"
)

if(!all(req_dr %in% names(dr))){
  stop("DRIPmissing: ",
       paste(setdiff(req_dr,names(dr)),collapse=", "))
}

dr[,SYMBOL:=canon(SYMBOL)]
dr <- dr[!duplicated(SYMBOL)]

dr[,Strict:=toupper(
  as.character(Strict_E2_DRIP)
) %in% c("TRUE","T","1")]

mrg <- merge(
  gene,
  dr[,.(SYMBOL,Mean_LFC_2h,Mean_LFC_24h,Strict)],
  by="SYMBOL",
  all=FALSE
)

test_drip_status <- function(metric_name,col){

  d <- mrg[
    is.finite(get(col)),
    .(Strict,Score=get(col))
  ]

  a <- d[Strict==TRUE,Score]
  b <- d[Strict==FALSE,Score]

  p <- if(length(a)>=3 && length(b)>=3){
    wilcox.test(
      a,b,
      alternative="two.sided",
      exact=FALSE
    )$p.value
  }else{
    NA_real_
  }

  data.table(
    Metric=metric_name,
    N_Strict=length(a),
    N_Background=length(b),
    Median_Strict=median(a,na.rm=TRUE),
    Median_Background=median(b,na.rm=TRUE),
    Delta_median=median(a,na.rm=TRUE)-median(b,na.rm=TRUE),
    P=p
  )
}

drip_status <- rbindlist(lapply(names(metric_map),function(m){
  test_drip_status(
    m,
    unname(metric_map[m])
  )
}))

drip_status[,FDR:=p.adjust(P,"BH")]

fwrite(
  drip_status,
  file.path(out_dir,"04_H2AZ_DRIPPolII.csv")
)

group_sets <- c(
  list("All genes"=unique(mrg$SYMBOL)),
  sets
)

cor_list <- list()
kk <- 1

for(grp in names(group_sets)){

  gg <- group_sets[[grp]]

  for(tm in c("2 h","24 h")){

    drip_col <- if(tm=="2 h"){
      "Mean_LFC_2h"
    }else{
      "Mean_LFC_24h"
    }

    for(mm in names(metric_map)){

      pol_col <- unname(metric_map[mm])

      d <- mrg[
        SYMBOL %chin% gg &
          is.finite(get(drip_col)) &
          is.finite(get(pol_col)),
        .(
          DRIP=get(drip_col),
          PolII=get(pol_col)
        )
      ]

      if(nrow(d)>=10){
        ct <- suppressWarnings(
          cor.test(
            d$DRIP,
            d$PolII,
            method="spearman",
            exact=FALSE
          )
        )

        rho <- unname(ct$estimate)
        p <- ct$p.value
      }else{
        rho <- p <- NA_real_
      }

      cor_list[[kk]] <- data.table(
        Group=grp,
        Time=tm,
        Metric=mm,
        N=nrow(d),
        Rho=rho,
        P=p
      )

      kk <- kk+1
    }
  }
}

cor_stat <- rbindlist(cor_list)
cor_stat[,FDR:=p.adjust(P,"BH")]

fwrite(
  cor_stat,
  file.path(out_dir,"05_H2AZ_DRIP.csv")
)

keep_cols <- c(
  "SYMBOL","REFSEQ",
  "TSS_ctrl","TSS_kd","TSS_input",
  "GB_ctrl","GB_kd","GB_input",
  "TSS_ctrl_bg","TSS_kd_bg",
  "GB_ctrl_bg","GB_kd_bg",
  "Delta_TSS_bg","Delta_GB_bg",
  "PI_ctrl_bg","PI_kd_bg","Delta_PI_bg",
  "Delta_TSS_norm","Delta_GB_norm","Delta_PI_norm",
  "PI_eligible"
)

c24 <- merge(
  data.table(SYMBOL=cand24),
  gene[,..keep_cols],
  by="SYMBOL",
  all.x=TRUE
)

final3 <- merge(
  data.table(
    SYMBOL=c("HMGB2","H2AZ2","UBALD2")
  ),
  gene[,..keep_cols],
  by="SYMBOL",
  all.x=TRUE
)

fwrite(
  c24,
  file.path(out_dir,"06_H2AZ_24candidategenePolIIresults.csv")
)

fwrite(
  final3,
  file.path(out_dir,"07_H2AZ_finalcandidateresults.csv")
)

pi_long <- rbind(
  gene[
    PI_eligible==TRUE &
      is.finite(PI_ctrl_bg) &
      PI_ctrl_bg>0,
    .(
      Method="Input-corrected PI",
      Condition="Control",
      PI=PI_ctrl_bg
    )
  ],
  gene[
    PI_eligible==TRUE &
      is.finite(PI_kd_bg) &
      PI_kd_bg>0,
    .(
      Method="Input-corrected PI",
      Condition="H2A.Z KD",
      PI=PI_kd_bg
    )
  ],
  gene[
    is.finite(PI_ctrl_norm) &
      PI_ctrl_norm>0,
    .(
      Method="Normalized-density ratio",
      Condition="Control",
      PI=PI_ctrl_norm
    )
  ],
  gene[
    is.finite(PI_kd_norm) &
      PI_kd_norm>0,
    .(
      Method="Normalized-density ratio",
      Condition="H2A.Z KD",
      PI=PI_kd_norm
    )
  ]
)

pi_long[,Condition:=factor(
  Condition,
  levels=c("Control","H2A.Z KD")
)]

p1 <- ggplot(
  pi_long,
  aes(Condition,log2(PI),fill=Condition)
)+
  geom_violin(
    trim=TRUE,
    scale="width",
    alpha=.40,
    linewidth=.35
  )+
  geom_boxplot(
    width=.18,
    outlier.shape=NA,
    linewidth=.45
  )+
  facet_wrap(~Method,scales="free_y")+
  scale_fill_manual(
    values=c(
      "Control"="#6fa6cf",
      "H2A.Z KD"="#b43665"
    )
  )+
  labs(
    x=NULL,
    y="log2 pausing index",
    title="Global RNAPII pausing after H2A.Z perturbation"
  )+
  theme_classic(base_size=9)+
  theme(
    plot.title=element_text(face="bold",hjust=.5),
    axis.text=element_text(color="black"),
    strip.background=element_blank(),
    strip.text=element_text(face="bold"),
    legend.position="none"
  )

savepdf("01_H2AZ_PI.pdf",p1,14,8.5)

set_stat[,GeneSet:=factor(
  GeneSet,
  levels=c(
    "24 candidates",
    "C5 stable-up",
    "RL-Sig65"
  )
)]

set_stat[,Metric:=factor(
  Metric,
  levels=c(
    "TSS RNAPII",
    "Gene-body RNAPII",
    "Pausing index"
  )
)]

set_stat[,Sig:=ifelse(
  is.na(FDR),"",
  ifelse(FDR<0.001,"***",
         ifelse(FDR<0.01,"**",
                ifelse(FDR<0.05,"*","ns")))
)]

p2 <- ggplot(
  set_stat,
  aes(
    GeneSet,
    Delta_median,
    color=Metric,
    group=Metric
  )
)+
  geom_hline(yintercept=0,linetype=2,color="grey60")+
  geom_point(
    position=position_dodge(width=.52),
    size=2.8
  )+
  geom_text(
    aes(label=Sig),
    position=position_dodge(width=.52),
    vjust=-1,
    size=2.8,
    show.legend=FALSE
  )+
  scale_color_manual(
    values=c(
      "TSS RNAPII"="#b43665",
      "Gene-body RNAPII"="#187d79",
      "Pausing index"="#efb421"
    )
  )+
  coord_flip()+
  labs(
    x=NULL,
    y="Δ median log2 change vs background",
    title="H2A.Z perturbation across study-defined gene sets",
    color=NULL
  )+
  theme_classic(base_size=9)+
  theme(
    plot.title=element_text(face="bold",hjust=.5),
    axis.text=element_text(color="black"),
    legend.position="bottom"
  )

savepdf("02_H2AZ_gene_setPolII.pdf",p2,14,8.5)

d3 <- mrg[
  is.finite(Delta_PI_bg)
]

d3[,DRIP_status:=ifelse(
  Strict,
  "Strict E2-DRIP",
  "Background"
)]

d3[,DRIP_status:=factor(
  DRIP_status,
  levels=c("Background","Strict E2-DRIP")
)]

pi_drip_row <- drip_status[
  Metric=="Pausing index"
]

subtxt <- if(nrow(pi_drip_row)==1){
  paste0(
    "Δmedian=",
    sprintf("%.3f",pi_drip_row$Delta_median),
    "; FDR=",
    signif(pi_drip_row$FDR,3)
  )
}else{
  NULL
}

p3 <- ggplot(
  d3,
  aes(DRIP_status,Delta_PI_bg,fill=DRIP_status)
)+
  geom_hline(yintercept=0,linetype=2,color="grey60")+
  geom_violin(
    trim=TRUE,
    scale="width",
    alpha=.38,
    linewidth=.35
  )+
  geom_boxplot(
    width=.18,
    outlier.shape=NA,
    linewidth=.45
  )+
  scale_fill_manual(
    values=c(
      "Background"="grey75",
      "Strict E2-DRIP"="#b43665"
    )
  )+
  labs(
    x=NULL,
    y="Δ log2 pausing index",
    title="H2A.Z perturbation in E2-DRIP-associated genes",
    subtitle=subtxt
  )+
  theme_classic(base_size=9)+
  theme(
    plot.title=element_text(face="bold",hjust=.5),
    plot.subtitle=element_text(hjust=.5),
    axis.text=element_text(color="black"),
    legend.position="none"
  )

savepdf("03A_H2AZ_DRIPPI.pdf",p3,10,8.5)

scatter <- mrg[
  is.finite(Mean_LFC_24h) &
    is.finite(Delta_PI_bg)
]

cr <- cor_stat[
  Group=="All genes" &
    Time=="24 h" &
    Metric=="Pausing index"
]

sub2 <- if(nrow(cr)==1){
  paste0(
    "Spearman ρ=",
    sprintf("%.3f",cr$Rho),
    "; FDR=",
    signif(cr$FDR,3)
  )
}else{
  NULL
}

p3b <- ggplot(
  scatter,
  aes(Mean_LFC_24h,Delta_PI_bg)
)+
  geom_hline(yintercept=0,linetype=2,color="grey75")+
  geom_vline(xintercept=0,linetype=2,color="grey75")+
  geom_point(
    size=.75,
    alpha=.22,
    color="#187d79"
  )+
  geom_smooth(
    method="lm",
    formula=y~x,
    se=FALSE,
    linewidth=.7,
    color="#b43665"
  )+
  labs(
    x="Mean E2-induced DRIP log2FC (24 h)",
    y="Δ log2 pausing index",
    title="E2-DRIP signal and H2A.Z-dependent RNAPII pausing",
    subtitle=sub2
  )+
  theme_classic(base_size=9)+
  theme(
    plot.title=element_text(face="bold",hjust=.5),
    plot.subtitle=element_text(hjust=.5),
    axis.text=element_text(color="black")
  )

savepdf("03B_H2AZ_DRIP.pdf",p3b,10,8.5)

hm <- melt(
  c24[,.(SYMBOL,Delta_TSS_bg,Delta_GB_bg,Delta_PI_bg)],
  id.vars="SYMBOL",
  variable.name="Metric",
  value.name="Value"
)

hm[,Metric:=factor(
  Metric,
  levels=c(
    "Delta_TSS_bg",
    "Delta_GB_bg",
    "Delta_PI_bg"
  ),
  labels=c(
    "TSS RNAPII",
    "Gene-body RNAPII",
    "Pausing index"
  )
)]

ord <- c24[
  order(-Delta_PI_bg,na.last=TRUE),
  SYMBOL
]

hm[,SYMBOL:=factor(SYMBOL,levels=rev(ord))]

p4 <- ggplot(
  hm,
  aes(Metric,SYMBOL,fill=Value)
)+
  geom_tile(color="white",linewidth=.35)+
  scale_fill_gradient2(
    low="#6fa6cf",
    mid="white",
    high="#b43665",
    midpoint=0,
    na.value="grey90"
  )+
  labs(
    x=NULL,
    y=NULL,
    fill="log2 change",
    title="RNAPII response of the 24 candidate genes"
  )+
  theme_classic(base_size=8.5)+
  theme(
    plot.title=element_text(face="bold",hjust=.5),
    axis.text=element_text(color="black"),
    legend.position="right"
  )

savepdf("04_H2AZ_24candidategeneplot.pdf",p4,10.5,13)

p_all <- p1 / p2 / p3 +
  plot_layout(heights=c(1,1,1))

savepdf(
  "05_H2AZ_plot.pdf",
  p_all,
  16.6,
  20
)
