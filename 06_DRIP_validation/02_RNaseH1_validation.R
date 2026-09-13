# GSE241307 RNase H1 sensitivity analysis

PROJECT_DIR <- "."

options(stringsAsFactors=FALSE); set.seed(1234)

if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
bio <- c("rtracklayer","GenomicRanges","GenomicFeatures","AnnotationDbi",
         "org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg38.knownGene")
for(p in bio) if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,ask=FALSE,update=FALSE)
for(p in c("data.table","ggplot2")) if(!requireNamespace(p,quietly=TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(rtracklayer)
  library(GenomicRanges); library(GenomicFeatures); library(AnnotationDbi)
  library(org.Hs.eg.db); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

base_dir <- file.path(PROJECT_DIR, "data", "raw", "GSE241307")
out_dir <- file.path(PROJECT_DIR, "results", "GSE241307")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

project <- PROJECT_DIR
rl65_file <- file.path(PROJECT_DIR, "gene_sets/RL_Sig65.txt")
c5_file <- file.path(PROJECT_DIR, "gene_sets/C5_stable_up_genes.txt")
cand_file <- file.path(PROJECT_DIR, "gene_sets/candidate24_genes.txt")

find_bw <- function(gsm){
  f <- list.files(base_dir,paste0("^",gsm,".*\\.bigWig$"),full.names=TRUE,ignore.case=TRUE)
  f <- f[file.exists(f) & !dir.exists(f)]
  if(length(f)!=1) stop(gsm," bigWignumber1, check.")
  f
}

bw_ev <- find_bw("GSM7720902")
bw_oe <- find_bw("GSM7720903")
bw_ctrl <- find_bw("GSM7720906")

need <- c(bw_ev,bw_oe,bw_ctrl,rl65_file,c5_file,cand_file)
if(any(!file.exists(need))) stop("missing: \n",paste(need[!file.exists(need)],collapse="\n"))

clean_gene <- function(x){
  x <- toupper(trimws(as.character(x))); x[x=="H2AFV"] <- "H2AZ2"
  unique(x[!is.na(x) & nzchar(x)])
}

rl65 <- clean_gene(readLines(rl65_file,warn=FALSE))
c5 <- clean_gene(readLines(c5_file,warn=FALSE))
cand24 <- clean_gene(readLines(cand_file,warn=FALSE))
final3 <- c("HMGB2","H2AZ2","UBALD2")
sets <- list(RL_Sig65=rl65,C5_stable_up=c5,Candidate_24=cand24)

gb <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                             single.strand.genes.only=TRUE)

keep_chr <- as.character(seqnames(gb)) %in% paste0("chr",c(1:22,"X","Y"))
gb <- gb[keep_chr]; entrez <- names(gb)

symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=entrez,column="SYMBOL",
                                keytype="ENTREZID",multiVals="first")
symbol <- toupper(as.character(symbol)); symbol[symbol=="H2AFV"] <- "H2AZ2"

ok <- !is.na(symbol) & nzchar(symbol)
gb <- gb[ok]; entrez <- entrez[ok]; symbol <- symbol[ok]

locus <- gb
plus <- as.character(strand(locus))=="+"
minus <- as.character(strand(locus))=="-"
start(locus)[plus] <- pmax(1,start(locus)[plus]-1000)
end(locus)[minus] <- end(locus)[minus]+1000

prom <- promoters(gb,upstream=1000,downstream=1000)

bw_chr <- seqlevels(BigWigFile(bw_ev))
ok <- as.character(seqnames(locus)) %in% bw_chr
gb <- gb[ok]; locus <- locus[ok]; prom <- prom[ok]
symbol <- symbol[ok]; entrez <- entrez[ok]

bw_mean <- function(file,regions){

  out <- rep(NA_real_,length(regions))
  bw <- BigWigFile(file)
  valid <- as.character(seqnames(regions)) %in% seqlevels(bw)
  out[valid] <- 0

  for(chr in unique(as.character(seqnames(regions[valid])))){

    ids <- which(valid & as.character(seqnames(regions))==chr)
    q <- regions[ids]; strand(q) <- "*"

    sig <- rtracklayer::import(bw,which=reduce(q))
    if(length(sig)==0) next
    strand(sig) <- "*"

    h <- findOverlaps(q,sig,ignore.strand=TRUE)
    if(length(h)==0) next

    qi <- queryHits(h); si <- subjectHits(h)
    ov_width <- pmin(end(q)[qi],end(sig)[si]) -
      pmax(start(q)[qi],start(sig)[si]) + 1

    score <- as.numeric(mcols(sig)$score[si])
    contribution <- score * ov_width

    sums <- rowsum(contribution,qi,reorder=FALSE)
    hit_id <- as.integer(rownames(sums))
    value <- numeric(length(q))
    value[hit_id] <- sums[,1] / width(q)[hit_id]

    out[ids] <- value
    rm(sig,h); gc(FALSE)
  }

  if(length(out)!=length(regions)) stop("BigWig.")
  out
}

s_ev <- bw_mean(bw_ev,locus)
s_oe <- bw_mean(bw_oe,locus)
s_ctrl <- bw_mean(bw_ctrl,locus)

raw <- data.table(SYMBOL=symbol,ENTREZID=entrez,EV=s_ev,
                  RNaseH1_OE=s_oe,RNase_treated=s_ctrl)

gene <- raw[is.finite(EV) & is.finite(RNaseH1_OE) & is.finite(RNase_treated),
            .(N_loci=.N,EV=mean(EV),RNaseH1_OE=mean(RNaseH1_OE),
              RNase_treated=mean(RNase_treated)),by=SYMBOL]

gene[,Delta_OE:=RNaseH1_OE-EV]
gene[,Delta_RNaseCtrl:=RNase_treated-EV]
gene[,RNaseH1_sensitive:=Delta_OE<0]
gene[,`:=`(RL_Sig65=SYMBOL%in%rl65,C5_stable_up=SYMBOL%in%c5,
           Candidate_24=SYMBOL%in%cand24,Final_3=SYMBOL%in%final3)]

setorder(gene,Delta_OE)
fwrite(gene,file.path(out_dir,"01_gene_level_RNaseH_signal_changes.csv"))

test_set <- function(name,g){

  g <- intersect(clean_gene(g),gene$SYMBOL)
  d <- gene[SYMBOL%in%g]; bg <- gene[!SYMBOL%in%g]

  p1 <- wilcox.test(d$EV,d$RNaseH1_OE,paired=TRUE,
                    alternative="greater",exact=FALSE)$p.value
  p2 <- wilcox.test(d$Delta_OE,bg$Delta_OE,
                    alternative="less",exact=FALSE)$p.value
  p3 <- wilcox.test(d$EV,d$RNase_treated,paired=TRUE,
                    alternative="greater",exact=FALSE)$p.value

  data.frame(GeneSet=name,N_evaluable=nrow(d),
             Median_EV=median(d$EV),Median_OE=median(d$RNaseH1_OE),
             Median_delta=median(d$Delta_OE),
             Fraction_decreased=mean(d$Delta_OE<0),
             Paired_P=p1,Background_delta=median(bg$Delta_OE),
             Selective_P=p2,Median_RNaseCtrl=median(d$RNase_treated),
             RNaseCtrl_P=p3)
}

stats <- rbindlist(lapply(names(sets),function(n) test_set(n,sets[[n]])))
stats[,Paired_FDR:=p.adjust(Paired_P,"BH")]
stats[,Selective_FDR:=p.adjust(Selective_P,"BH")]
stats[,RNaseCtrl_FDR:=p.adjust(RNaseCtrl_P,"BH")]

fwrite(stats,file.path(out_dir,"02_gene_set_RNaseH_sensitivity_statistics.csv"))

pick_gene <- function(g){
  id <- which(symbol==g)
  if(length(id)==0) return(NA_integer_)
  id[which.max(width(gb[id]))]
}

ids <- sapply(final3,pick_gene)
valid <- !is.na(ids)
ids2 <- as.integer(ids[valid])

res3 <- data.table(
  SYMBOL=final3[valid],
  Promoter_EV=bw_mean(bw_ev,prom[ids2]),
  Promoter_OE=bw_mean(bw_oe,prom[ids2]),
  Promoter_RNaseCtrl=bw_mean(bw_ctrl,prom[ids2]),
  GeneBody_EV=bw_mean(bw_ev,gb[ids2]),
  GeneBody_OE=bw_mean(bw_oe,gb[ids2]),
  GeneBody_RNaseCtrl=bw_mean(bw_ctrl,gb[ids2])
)

res3[,Promoter_delta:=Promoter_OE-Promoter_EV]
res3[,GeneBody_delta:=GeneBody_OE-GeneBody_EV]
res3[,Promoter_sensitive:=Promoter_delta<0]
res3[,GeneBody_sensitive:=GeneBody_delta<0]

candidate <- merge(data.table(SYMBOL=final3),res3,by="SYMBOL",all.x=TRUE,sort=FALSE)
candidate <- candidate[match(final3,SYMBOL)]

if(nrow(candidate)!=3) stop("candidategeneresults3, analysis.")

fwrite(candidate,file.path(out_dir,"03_candidate_gene_RNaseH_results.csv"))

plot_df <- rbindlist(lapply(names(sets),function(n){
  d <- gene[SYMBOL%in%sets[[n]],.(SYMBOL,Delta_OE,Delta_RNaseCtrl)]
  d[,GeneSet:=n]; d
}))

plot_long <- melt(plot_df,id.vars=c("SYMBOL","GeneSet"),
                  measure.vars=c("Delta_OE","Delta_RNaseCtrl"),
                  variable.name="Comparison",value.name="Delta")

plot_long[Comparison=="Delta_OE",Comparison:="RNase H1 OE - EV"]
plot_long[Comparison=="Delta_RNaseCtrl",Comparison:="RNase-treated - EV"]
plot_long$GeneSet <- factor(plot_long$GeneSet,
                            levels=c("RL_Sig65","C5_stable_up","Candidate_24"))

p <- ggplot(plot_long,aes(GeneSet,Delta)) +
  geom_hline(yintercept=0,linetype=2,linewidth=.35) +
  geom_violin(trim=TRUE,linewidth=.35) +
  geom_boxplot(width=.18,outlier.shape=NA,linewidth=.4) +
  facet_wrap(~Comparison,ncol=1,scales="free_y") +
  labs(x=NULL,y="Change in DRIP signal") +
  theme_classic(base_size=8) +
  theme(strip.background=element_blank(),strip.text=element_text(face="bold"),
        axis.text.x=element_text(angle=20,hjust=1))

ggsave(file.path(out_dir,"04_RNaseHplot.pdf"),
       p,width=8.5,height=10,units="cm")

hm <- pick_gene("HMGB2")

if(!is.na(hm)){

  region <- gb[hm]
  start(region) <- pmax(1,start(region)-5000)
  end(region) <- end(region)+5000

  bins <- unlist(tile(region,n=250))
  pos <- (start(bins)+end(bins))/2
  hm_strand <- as.character(strand(gb[hm]))
  tss <- if(hm_strand=="+") start(gb[hm]) else end(gb[hm])
  pos_kb <- if(hm_strand=="+") (pos-tss)/1000 else (tss-pos)/1000

  prof <- rbind(
    data.table(Position_kb=pos_kb,Signal=bw_mean(bw_ev,bins),Group="EV"),
    data.table(Position_kb=pos_kb,Signal=bw_mean(bw_oe,bins),Group="RNase H1 OE"),
    data.table(Position_kb=pos_kb,Signal=bw_mean(bw_ctrl,bins),Group="RNase-treated")
  )

  prof$Group <- factor(prof$Group,levels=c("EV","RNase H1 OE","RNase-treated"))

  p2 <- ggplot(prof,aes(Position_kb,Signal,linetype=Group)) +
    annotate("rect",xmin=-1,xmax=1,ymin=-Inf,ymax=Inf,alpha=.08) +
    geom_line(linewidth=.65) + geom_vline(xintercept=0,linetype=3,linewidth=.35) +
    labs(x="Position relative to HMGB2 TSS (kb)",y="DRIP signal",
         linetype=NULL,title="HMGB2 locus") +
    theme_classic(base_size=8) +
    theme(plot.title=element_text(hjust=.5,face="bold"),legend.position="top")

  ggsave(file.path(out_dir,"05_HMGB2DRIPplot.pdf"),
         p2,width=9,height=6.5,units="cm")
}

print(stats[,.(GeneSet,N_evaluable,Median_EV,Median_OE,Median_delta,
               Fraction_decreased,Paired_FDR,Selective_FDR,RNaseCtrl_FDR)])

print(candidate)

if(!requireNamespace("patchwork",quietly=TRUE)) install.packages("patchwork")
suppressPackageStartupMessages(library(patchwork))

col_set <- c("RL-Sig65"="#b43665","C5 stable-up"="#6fa6cf","24 candidates"="#8560af")
col_group <- c("EV"="#187d79","RNase H1 OE"="#efb421","RNase-treated"="#e58027")

fmt_fdr <- function(x){
  out <- rep("NA",length(x)); ok <- !is.na(x)
  out[ok & x < 1e-99] <- "<1e-99"
  out[ok & x >= 1e-99 & x < 0.001] <- formatC(x[ok & x >= 1e-99 & x < 0.001],format="e",digits=1)
  out[ok & x >= 0.001] <- sprintf("%.3f",x[ok & x >= 0.001])
  out
}

theme_pub <- theme_classic(base_size=8) +
  theme(
    axis.text=element_text(size=7,colour="black"),
    axis.title=element_text(size=8,colour="black"),
    axis.line=element_line(linewidth=.4,colour="black"),
    axis.ticks=element_line(linewidth=.35,colour="black"),
    legend.title=element_blank(),
    legend.text=element_text(size=7),
    strip.background=element_blank(),
    strip.text=element_text(size=8,face="bold"),
    plot.title=element_text(size=9,face="bold",hjust=.5),
    plot.margin=margin(5,6,5,6)
  )

plotA <- copy(stats)
plotA[,GeneSet_label:=c("RL-Sig65","C5 stable-up","24 candidates")]
plotA[,Decrease_label:=paste0(sprintf("%.1f",Fraction_decreased*100),"% decreased")]
plotA[,FDR_label:=paste0("Selective FDR ",fmt_fdr(Selective_FDR))]
plotA[,Label:=paste0(Decrease_label,"\n",FDR_label)]

plotA$GeneSet_label <- factor(plotA$GeneSet_label,
                              levels=c("24 candidates","C5 stable-up","RL-Sig65"))

x_all <- c(plotA$Median_delta,plotA$Background_delta,0)
x_range <- diff(range(x_all,na.rm=TRUE))
if(x_range==0) x_range <- 1
label_x <- max(x_all,na.rm=TRUE) + 0.10*x_range

pA <- ggplot(plotA,aes(y=GeneSet_label)) +
  geom_vline(xintercept=0,linetype=2,linewidth=.4,colour="grey55") +
  geom_segment(aes(x=Background_delta,xend=Median_delta,yend=GeneSet_label),
               linewidth=.8,colour="grey70") +
  geom_point(aes(x=Background_delta),shape=23,size=2.6,stroke=.7,fill="white",colour="grey35") +
  geom_point(aes(x=Median_delta,colour=GeneSet_label),size=3.2) +
  geom_text(aes(x=label_x,label=Label),hjust=0,size=2.45,lineheight=.95) +
  scale_colour_manual(values=col_set) +
  scale_x_continuous(expand=expansion(mult=c(.08,.42))) +
  labs(
    x="Change in DRIP signal (RNase H1 OE - EV)",
    y=NULL,
    title="RNase H1-sensitive transcriptional programs"
  ) +
  coord_cartesian(clip="off") +
  theme_pub +
  theme(legend.position="none")

ggsave(file.path(out_dir,"04_RNaseHgene_set_publication.pdf"),
       pA,width=9.2,height=6.8,units="cm",device="pdf")

cand_long <- rbind(
  data.table(SYMBOL=candidate$SYMBOL,Region="Promoter",Group="EV",Signal=candidate$Promoter_EV),
  data.table(SYMBOL=candidate$SYMBOL,Region="Promoter",Group="RNase H1 OE",Signal=candidate$Promoter_OE),
  data.table(SYMBOL=candidate$SYMBOL,Region="Promoter",Group="RNase-treated",Signal=candidate$Promoter_RNaseCtrl),
  data.table(SYMBOL=candidate$SYMBOL,Region="Gene body",Group="EV",Signal=candidate$GeneBody_EV),
  data.table(SYMBOL=candidate$SYMBOL,Region="Gene body",Group="RNase H1 OE",Signal=candidate$GeneBody_OE),
  data.table(SYMBOL=candidate$SYMBOL,Region="Gene body",Group="RNase-treated",Signal=candidate$GeneBody_RNaseCtrl)
)

cand_long$SYMBOL <- factor(cand_long$SYMBOL,levels=c("UBALD2","H2AZ2","HMGB2"))
cand_long$Region <- factor(cand_long$Region,levels=c("Promoter","Gene body"))
cand_long$Group <- factor(cand_long$Group,levels=c("EV","RNase H1 OE","RNase-treated"))

cand_seg <- rbind(
  data.table(SYMBOL=candidate$SYMBOL,Region="Promoter",
             EV=candidate$Promoter_EV,OE=candidate$Promoter_OE),
  data.table(SYMBOL=candidate$SYMBOL,Region="Gene body",
             EV=candidate$GeneBody_EV,OE=candidate$GeneBody_OE)
)

cand_seg$SYMBOL <- factor(cand_seg$SYMBOL,levels=c("UBALD2","H2AZ2","HMGB2"))
cand_seg$Region <- factor(cand_seg$Region,levels=c("Promoter","Gene body"))

pB <- ggplot(cand_long,aes(x=Signal,y=SYMBOL)) +
  geom_segment(data=cand_seg,aes(x=EV,xend=OE,y=SYMBOL,yend=SYMBOL),
               inherit.aes=FALSE,linewidth=.7,colour="grey65") +
  geom_point(aes(colour=Group,shape=Group),size=2.6,stroke=.7) +
  facet_wrap(~Region,ncol=1,scales="free_x") +
  scale_colour_manual(values=col_group) +
  scale_shape_manual(values=c("EV"=16,"RNase H1 OE"=17,"RNase-treated"=15)) +
  labs(x="Average DRIP signal",y=NULL,title="Candidate gene loci") +
  theme_pub +
  theme(
    legend.position="top",
    legend.direction="horizontal",
    panel.spacing=unit(.5,"lines")
  )

ggsave(file.path(out_dir,"05_candidate_gene_RNaseH_results_publication.pdf"),
       pB,width=8.2,height=8.5,units="cm",device="pdf")

hm <- pick_gene("HMGB2")

if(is.na(hm)) stop("not_foundHMGB2 hg38gene.")

region <- gb[hm]
start(region) <- pmax(1,start(region)-5000)
end(region) <- end(region)+5000

bins <- unlist(tile(region,n=300))
pos <- (start(bins)+end(bins))/2

hm_strand <- as.character(strand(gb[hm]))
tss <- if(hm_strand=="+") start(gb[hm]) else end(gb[hm])
pos_kb <- if(hm_strand=="+") (pos-tss)/1000 else (tss-pos)/1000

prof <- rbind(
  data.table(Position_kb=pos_kb,Signal=bw_mean(bw_ev,bins),Group="EV"),
  data.table(Position_kb=pos_kb,Signal=bw_mean(bw_oe,bins),Group="RNase H1 OE"),
  data.table(Position_kb=pos_kb,Signal=bw_mean(bw_ctrl,bins),Group="RNase-treated")
)

prof$Group <- factor(prof$Group,levels=c("EV","RNase H1 OE","RNase-treated"))

pC <- ggplot(prof,aes(Position_kb,Signal,colour=Group,linetype=Group)) +
  annotate("rect",xmin=-1,xmax=1,ymin=-Inf,ymax=Inf,fill="#b43665",alpha=.07) +
  geom_vline(xintercept=0,linetype=3,linewidth=.4,colour="grey30") +
  geom_line(linewidth=.8) +
  scale_colour_manual(values=col_group) +
  scale_linetype_manual(values=c("EV"="solid","RNase H1 OE"="longdash","RNase-treated"="dotted")) +
  annotate("text",x=0,y=Inf,label="TSS",vjust=1.3,size=2.4,colour="grey25") +
  labs(
    x="Position relative to HMGB2 TSS (kb)",
    y="DRIP signal",
    title="HMGB2 locus"
  ) +
  theme_pub +
  theme(
    legend.position="top",
    legend.direction="horizontal"
  )

ggsave(file.path(out_dir,"06_HMGB2DRIPplot_publication.pdf"),
       pC,width=9.5,height=6.8,units="cm",device="pdf")

final_fig <- (pA | pB) / pC +
  plot_layout(widths=c(1.08,.92),heights=c(1,0.82)) +
  plot_annotation(tag_levels="A") &
  theme(
    plot.tag=element_text(size=10,face="bold",colour="black"),
    plot.tag.position=c(.01,.99)
  )

ggsave(file.path(out_dir,"07_RNaseH_publication.pdf"),
       final_fig,width=17,height=13.5,units="cm",device="pdf")

options(stringsAsFactors=FALSE); set.seed(1234)

if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
for(p in c("data.table","ggplot2","scales"))
  if(!requireNamespace(p,quietly=TRUE)) install.packages(p)
for(p in c("rtracklayer","GenomicRanges","GenomicFeatures","AnnotationDbi",
           "org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg38.knownGene"))
  if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,ask=FALSE,update=FALSE)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
  library(rtracklayer); library(GenomicRanges); library(GenomicFeatures)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

pdf_dev <- if(capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

project <- PROJECT_DIR
base_dir <- file.path(PROJECT_DIR, "data", "raw", "GSE241307")
out_dir  <- file.path(PROJECT_DIR, "results", "GSE241307")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

find_bw <- function(gsm){
  f <- list.files(base_dir,paste0("^",gsm,".*\\.bigWig$"),
                  full.names=TRUE,recursive=TRUE,ignore.case=TRUE)
  f <- f[file.exists(f) & !dir.exists(f)]
  if(length(f)==0) stop("not_found ",gsm," bigWig.")
  if(length(f)>1) message(gsm," , : ",f[1])
  f[1]
}

bw_ev   <- find_bw("GSM7720902")
bw_oe   <- find_bw("GSM7720903")
bw_ctrl <- find_bw("GSM7720906")

genes <- c("HMGB2","H2AZ2","UBALD2")
UP_KB <- 5; DOWN_KB <- 8; N_BIN <- 260L

col_group <- c(
  "EV"="#187d79",
  "RNase H1 OE"="#efb421",
  "RNase-treated control"="#b43665"
)
lty_group <- c(
  "EV"="solid",
  "RNase H1 OE"="longdash",
  "RNase-treated control"="dotted"
)

gb <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene,single.strand.genes.only=TRUE)
gb <- gb[as.character(seqnames(gb)) %in% paste0("chr",c(1:22,"X","Y"))]
entrez <- names(gb)

symbol <- mapIds(org.Hs.eg.db,keys=entrez,column="SYMBOL",
                 keytype="ENTREZID",multiVals="first")
symbol <- toupper(as.character(symbol))
symbol[symbol=="H2AFV"] <- "H2AZ2"

ok <- !is.na(symbol) & nzchar(symbol)
gb <- gb[ok]; symbol <- symbol[ok]

bw_chr <- seqlevels(BigWigFile(bw_ev))
ok <- as.character(seqnames(gb)) %in% bw_chr
gb <- gb[ok]; symbol <- symbol[ok]

pick_gene <- function(g){
  id <- which(symbol==g)
  if(length(id)==0) return(NA_integer_)
  id[which.max(width(gb[id]))]
}

ids <- sapply(genes,pick_gene)
if(any(is.na(ids))) stop("hg38annotationmissing: ",paste(genes[is.na(ids)],collapse=", "))

bw_mean <- function(file,regions){
  out <- rep(NA_real_,length(regions))
  bw <- BigWigFile(file)
  valid <- as.character(seqnames(regions)) %in% seqlevels(bw)
  out[valid] <- 0

  for(chr in unique(as.character(seqnames(regions[valid])))){
    ii <- which(valid & as.character(seqnames(regions))==chr)
    q <- regions[ii]; strand(q) <- "*"

    sig <- rtracklayer::import(bw,which=reduce(q))
    if(length(sig)==0) next
    strand(sig) <- "*"

    h <- findOverlaps(q,sig,ignore.strand=TRUE)
    if(length(h)==0) next

    qi <- queryHits(h); si <- subjectHits(h)
    ov <- pmin(end(q)[qi],end(sig)[si]) - pmax(start(q)[qi],start(sig)[si]) + 1
    contribution <- as.numeric(mcols(sig)$score[si]) * ov

    z <- rowsum(contribution,qi,reorder=FALSE)
    hit <- as.integer(rownames(z))
    value <- numeric(length(q))
    value[hit] <- z[,1] / width(q)[hit]
    out[ii] <- value
  }
  out
}

make_profile <- function(g){
  gid <- pick_gene(g)
  gr <- gb[gid]
  st <- as.character(strand(gr))
  tss <- if(st=="+") start(gr) else end(gr)

  if(st=="+"){
    s <- max(1,tss-UP_KB*1000); e <- tss+DOWN_KB*1000
  } else {
    s <- max(1,tss-DOWN_KB*1000); e <- tss+UP_KB*1000
  }

  win <- GRanges(seqnames=seqnames(gr),ranges=IRanges(s,e),strand="*")
  bins <- unlist(tile(win,n=N_BIN))
  mid <- (start(bins)+end(bins))/2

  pos <- if(st=="+") (mid-tss)/1000 else (tss-mid)/1000

  rbind(
    data.table(Gene=g,Position_kb=pos,Signal=bw_mean(bw_ev,bins),Group="EV"),
    data.table(Gene=g,Position_kb=pos,Signal=bw_mean(bw_oe,bins),Group="RNase H1 OE"),
    data.table(Gene=g,Position_kb=pos,Signal=bw_mean(bw_ctrl,bins),
               Group="RNase-treated control")
  )
}

prof <- rbindlist(lapply(genes,make_profile))
prof <- prof[is.finite(Position_kb) & is.finite(Signal)]
prof[,Gene:=factor(Gene,levels=genes)]
prof[,Group:=factor(Group,levels=names(col_group))]

fwrite(prof,file.path(out_dir,"Figure6E_candidate_gene_TSS_DRIP_profiles.csv"))

p <- ggplot(prof,aes(Position_kb,Signal,colour=Group,linetype=Group)) +
  annotate("rect",xmin=-1,xmax=1,ymin=-Inf,ymax=Inf,
           fill="#b43665",alpha=.055) +
  geom_vline(xintercept=0,linetype=3,linewidth=.4,colour="grey30") +
  geom_line(linewidth=.72,na.rm=TRUE) +
  facet_wrap(~Gene,nrow=1,scales="free_y") +
  scale_colour_manual(values=col_group,name=NULL) +
  scale_linetype_manual(values=lty_group,name=NULL) +
  scale_x_continuous(limits=c(-UP_KB,DOWN_KB),breaks=c(-5,0,5),expand=c(0,0)) +
  labs(x="Position relative to TSS (kb)",y="DRIP signal") +
  theme_classic(base_size=8) +
  theme(
    axis.text=element_text(size=7,colour="black"),
    axis.title=element_text(size=8,colour="black"),
    axis.line=element_line(linewidth=.45,colour="black"),
    axis.ticks=element_line(linewidth=.4,colour="black"),
    strip.background=element_blank(),
    strip.text=element_text(size=8.5,face="bold.italic",colour="black"),
    legend.position="top",
    legend.direction="horizontal",
    legend.text=element_text(size=7),
    legend.key.width=grid::unit(.9,"cm"),
    panel.spacing=grid::unit(.75,"lines"),
    panel.grid.major.y=element_line(colour="grey93",linewidth=.25),
    panel.grid.minor=element_blank(),
    plot.margin=margin(4,5,4,5)
  )

ggsave(file.path(out_dir,"Figure6E_candidategene_TSSlocus_.pdf"),
       p,width=16,height=5.8,units="cm",device=pdf_dev)

for(g in genes){
  pg <- ggplot(prof[Gene==g],
               aes(Position_kb,Signal,colour=Group,linetype=Group)) +
    annotate("rect",xmin=-1,xmax=1,ymin=-Inf,ymax=Inf,
             fill="#b43665",alpha=.055) +
    geom_vline(xintercept=0,linetype=3,linewidth=.4,colour="grey30") +
    geom_line(linewidth=.78,na.rm=TRUE) +
    scale_colour_manual(values=col_group,name=NULL) +
    scale_linetype_manual(values=lty_group,name=NULL) +
    scale_x_continuous(limits=c(-UP_KB,DOWN_KB),breaks=c(-5,0,5),expand=c(0,0)) +
    labs(x=paste0("Position relative to ",g," TSS (kb)"),
         y="DRIP signal",title=paste0(g," locus")) +
    theme_classic(base_size=8) +
    theme(
      axis.text=element_text(size=7,colour="black"),
      axis.title=element_text(size=8,colour="black"),
      plot.title=element_text(size=9,face="bold",hjust=.5),
      legend.position="top",
      legend.direction="horizontal",
      legend.text=element_text(size=7),
      legend.key.width=grid::unit(.9,"cm"),
      panel.grid.major.y=element_line(colour="grey93",linewidth=.25),
      panel.grid.minor=element_blank()
    )

  ggsave(file.path(out_dir,paste0("Figure6E_",g,"_TSSlocus.pdf")),
         pg,width=8.2,height=5.8,units="cm",device=pdf_dev)
}
