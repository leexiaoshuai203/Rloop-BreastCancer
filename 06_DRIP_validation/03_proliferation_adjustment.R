# Pan-proliferation competition and DRIP specificity analysis

#### broad proliferation PC1 + gene-level adjustment + matched null
#### + within-high-proliferation validation
rm(list=ls()); gc(); options(stringsAsFactors=FALSE); set.seed(20260816)

project <- "."
OUT <- file.path(project, "results", "DRIP_validation", "pan_proliferation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
INPUT_DIR <- file.path(project, "data", "processed", "NAC")
RL_FILE <- file.path(project, "gene_sets", "RL_Sig65.txt")
C5_FILE <- file.path(project, "gene_sets", "C5_stable_up_genes.txt")
GMT_FILE <- file.path(project, "data", "gene_sets", "h.all.v2026.1.Hs.symbols.gmt")
DRIP_FILE <- file.path(
  project, "results", "GSE81851", "enrichment",
  "03_gene_level_E2_DRIP_effects.csv"
)
need <- c(RL_FILE,C5_FILE,GMT_FILE,DRIP_FILE)
if (any(!file.exists(need))) stop("Missing required files:\n", paste(need[!file.exists(need)], collapse="\n"))

cran <- c("data.table","ggplot2","patchwork","metafor","brglm2","openxlsx","scales","msigdbr")
bioc <- c("GSVA","GSEABase","TxDb.Hsapiens.UCSC.hg19.knownGene","org.Hs.eg.db","GenomicFeatures",
          "GenomicRanges","AnnotationDbi","BSgenome.Hsapiens.UCSC.hg19","BSgenome","Biostrings")
if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
for(p in cran) if(!requireNamespace(p,quietly=TRUE)) install.packages(p,dependencies=TRUE)
for(p in bioc) if(!requireNamespace(p,quietly=TRUE)) BiocManager::install(p,ask=FALSE,update=FALSE)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(metafor); library(brglm2)
  library(openxlsx); library(scales); library(GSVA); library(GSEABase)
})

PAL <- c("RL-Sig65"="#b43665","C5 stable-up"="#187d79","Pan-proliferation"="#7b6ca8",
         "MKI67"="#777777","PCNA"="#6fa6cf")
theme_pub <- function(s=9) theme_classic(base_size=s,base_family="Arial")+theme(
  axis.text=element_text(color="black"),axis.title=element_text(face="bold",color="black"),
  strip.background=element_blank(),strip.text=element_text(face="bold"),
  panel.grid.major.y=element_line(color="#EEEEEE",linewidth=.25),panel.grid.minor=element_blank())
savepdf <- function(n,p,w,h) ggsave(file.path(OUT,n),p,width=w,height=h,units="cm",device=cairo_pdf,bg="white")
canon <- function(x){x<-toupper(trimws(as.character(x))); x[x=="H2AFV"]<-"H2AZ2"; x}
clean <- function(x){x<-canon(x); unique(x[!is.na(x)&nzchar(x)])}
zscore <- function(x){x<-as.numeric(x); s<-sd(x,na.rm=TRUE); if(!is.finite(s)||s==0)return(rep(0,length(x))); (x-mean(x,na.rm=TRUE))/s}
get_gse <- function(x){h<-regmatches(x,regexpr("GSE[0-9]+",x,ignore.case=TRUE)); h<-toupper(h); ifelse(nchar(h)>0,h,NA_character_)}
get_resp <- function(s,f=""){
  r<-rep(NA_character_,length(s)); r[grepl("([_.-])pCR$",s,ignore.case=TRUE)]<-"pCR"; r[grepl("([_.-])RD$",s,ignore.case=TRUE)]<-"RD"
  if(all(is.na(r))){if(grepl("pCR",f,ignore.case=TRUE))r[]<-"pCR"; if(grepl("([_.-])RD",f,ignore.case=TRUE))r[]<-"RD"}; r
}
collapse_dup <- function(m){
  if(!anyDuplicated(rownames(m)))return(m)
  sm<-rowsum(m,group=rownames(m),reorder=FALSE,na.rm=TRUE); n<-as.numeric(table(factor(rownames(m),levels=rownames(sm)))); sm/n
}
read_expr <- function(path){
  d<-fread(path,data.table=FALSE,check.names=FALSE,showProgress=FALSE); g<-canon(d[[1]])
  m<-as.matrix(d[,-1,drop=FALSE]); suppressWarnings(mode(m)<-"numeric"); rownames(m)<-g
  m<-m[!is.na(rownames(m))&rownames(m)!="",,drop=FALSE]; m<-m[rowMeans(is.na(m))<=.5,,drop=FALSE]
  for(i in seq_len(nrow(m))){q<-is.na(m[i,]); if(any(q)){v<-median(m[i,],na.rm=TRUE); if(!is.finite(v))v<-0; m[i,q]<-v}}
  m<-collapse_dup(m); cc<-get_gse(basename(path))
  if(is.na(cc)){u<-unique(na.omit(get_gse(colnames(m)))); if(length(u)!=1)stop("Unable to identify a GEO accession for: ", path); cc<-u}
  sn<-colnames(m); hit<-grepl(cc,sn,ignore.case=TRUE); sn[!hit]<-paste0(cc,"__",sn[!hit]); sn<-make.unique(sn); colnames(m)<-sn
  r<-get_resp(sn,basename(path)); if(any(is.na(r)))stop("Unable to identify pCR/RD labels in: ", basename(path))
  list(cohort=cc,expression=m,response=setNames(r,sn))
}
combine_cohort <- function(objs){
  out<-list()
  for(cc in unique(vapply(objs,`[[`,character(1),"cohort"))){
    x<-objs[vapply(objs,function(z)z$cohort==cc,logical(1))]
    if(length(x)==1){out[[cc]]<-x[[1]]; next}
    g<-Reduce(intersect,lapply(x,function(z)rownames(z$expression)))
    mm<-do.call(cbind,lapply(x,function(z)z$expression[g,,drop=FALSE])); rr<-unlist(lapply(x,`[[`,"response"),use.names=TRUE)
    if(anyDuplicated(colnames(mm))){nn<-make.unique(colnames(mm)); names(rr)<-nn; colnames(mm)<-nn}
    out[[cc]]<-list(cohort=cc,expression=mm,response=rr[colnames(mm)])
  }; out
}
read_gmt <- function(f){
  g<-GSEABase::getGmt(f,geneIdType=GSEABase::SymbolIdentifier()); x<-lapply(g,GSEABase::geneIds)
  names(x)<-vapply(g,GSEABase::setName,character(1)); x
}
ssgsea <- function(m,sets){
  ss<-lapply(sets,function(x)intersect(x,rownames(m)))
  if(any(lengths(ss)<5)) stop("Gene sets with fewer than five matched genes: ", paste(names(ss)[lengths(ss)<5], collapse=", "))
  if(exists("ssgseaParam",where=asNamespace("GSVA"),inherits=FALSE)){
    p<-GSVA::ssgseaParam(exprData=m,geneSets=ss,alpha=.25,normalize=FALSE); as.matrix(GSVA::gsva(p,verbose=FALSE))
  }else as.matrix(GSVA::gsva(expr=m,gset.idx.list=ss,method="ssgsea",ssgsea.norm=FALSE,verbose=FALSE))
}
safe_fit <- function(d,form){
  vars<-all.vars(form); d<-d[complete.cases(d[,..vars])]
  if(nrow(d)<30||uniqueN(d$Strict)<2)return(NULL)
  tryCatch(glm(form,data=d,family=binomial(),method=brglm2::brglmFit),error=function(e)NULL)
}
term_row <- function(fit,term){
  if(is.null(fit))return(NULL); q<-summary(fit)$coefficients; if(!term%in%rownames(q))return(NULL)
  b<-q[term,"Estimate"]; se<-q[term,"Std. Error"]
  data.table(beta=b,SE=se,OR=exp(b),CI_low=exp(b-1.96*se),CI_high=exp(b+1.96*se),P=q[term,"Pr(>|z|)"])
}
meta_fit <- function(scores,var,label){
  rr<-rbindlist(lapply(unique(scores$Cohort),function(cc){
    d<-scores[Cohort==cc,.(y,x=get(var))]; if(nrow(d)<10||uniqueN(d$y)<2)return(NULL)
    f<-tryCatch(glm(y~x,data=d,family=binomial(),method=brglm2::brglmFit),error=function(e)NULL)
    if(is.null(f))return(NULL); q<-summary(f)$coefficients["x",]; data.table(Cohort=cc,beta=q["Estimate"],SE=q["Std. Error"])
  }),fill=TRUE)
  if(nrow(rr)<2)return(data.table(Effect=label,OR=NA,CI_low=NA,CI_high=NA,P=NA,I2=NA))
  m<-metafor::rma.uni(yi=rr$beta,sei=rr$SE,method="REML")
  data.table(Effect=label,OR=exp(as.numeric(m$b)),CI_low=exp(m$ci.lb),CI_high=exp(m$ci.ub),P=as.numeric(m$pval),I2=as.numeric(m$I2))
}

rl65<-clean(readLines(RL_FILE,warn=FALSE)); c5<-clean(readLines(C5_FILE,warn=FALSE))
hall<-read_gmt(GMT_FILE); names(hall)<-gsub("^HALLMARK_","",names(hall))
cin70<-clean(c("TPX2","PRC1","FOXM1","CDK1","TGIF2","MCM2","H2AFZ","TOP2A","PCNA","UBE2C","MELK","TRIP13",
               "NCAPD2","MCM7","RNASEH2A","RAD51AP1","KIF20A","CDC45","MAD2L1","ESPL1","CCNB2","FEN1","TTK","CCT5",
               "RFC4","ATAD2","CKAP5","NUP205","CDC20","CKS2","RRM2","ELAVL1","CCNB1","RRM1","AURKB","MSH6","EZH2",
               "CTPS1","DKC1","OIP5","CDCA8","PTTG1","CEP55","H2AFX","CMAS","NCAPH","MCM10","LSM4","NCAPG2","ASF1B",
               "ZWINT","PBK","CDCA3","ECT2","CDC6","UNG","MTCH2","RAD21","ACTL6A","GPI","SRSF2","HDGF","NXT1","NEK2",
               "DHCR7","AURKA","NDUFAB1","MIIP","KIF4A"))

cols<-as.data.table(msigdbr::msigdbr_collections())
ccol<-if("gs_collection"%in%names(cols))"gs_collection" else "gs_cat"
scol<-if("gs_subcollection"%in%names(cols))"gs_subcollection" else "gs_subcat"
find_sub<-function(collection,patterns){
  x<-unique(cols[get(ccol)==collection,get(scol)]); x<-x[!is.na(x)]
  for(pt in patterns){h<-x[grepl(pt,x,ignore.case=TRUE)]; if(length(h))return(h[1])}
  stop("Unable to identify the requested MSigDB subcollection: ", collection, " / ", paste(patterns, collapse="/"))
}
msig_call<-function(collection,subcollection){
  fm<-names(formals(msigdbr::msigdbr)); a<-list(species="Homo sapiens")
  if("collection"%in%fm)a$collection<-collection else a$category<-collection
  if("subcollection"%in%fm)a$subcollection<-subcollection else a$subcategory<-subcollection
  as.data.table(do.call(msigdbr::msigdbr,a))
}
extract_set<-function(tab,exact){z<-tab[toupper(gs_name)==toupper(exact)]; if(!nrow(z))stop("MSigDB gene set not found: ", exact); clean(z$gene_symbol)}
m_kegg<-msig_call("C2",find_sub("C2",c("KEGG_LEGACY","KEGG")))
m_react<-msig_call("C2",find_sub("C2","REACTOME"))
m_gobp<-msig_call("C5",find_sub("C5",c("GO:BP","BP")))

pan_sets<-list(
  CIN70=cin70,H_E2F=clean(hall[["E2F_TARGETS"]]),H_G2M=clean(hall[["G2M_CHECKPOINT"]]),
  H_MYC_V1=clean(hall[["MYC_TARGETS_V1"]]),H_MYC_V2=clean(hall[["MYC_TARGETS_V2"]]),
  H_MITOTIC_SPINDLE=clean(hall[["MITOTIC_SPINDLE"]]),
  KEGG_CELL_CYCLE=extract_set(m_kegg,"KEGG_CELL_CYCLE"),
  REACTOME_CELL_CYCLE=extract_set(m_react,"REACTOME_CELL_CYCLE"),
  GOBP_MITOTIC_CELL_CYCLE=extract_set(m_gobp,"GOBP_MITOTIC_CELL_CYCLE"),
  GOBP_CHROMOSOME_SEGREGATION=extract_set(m_gobp,"GOBP_CHROMOSOME_SEGREGATION"),
  GOBP_DNA_REPLICATION=extract_set(m_gobp,"GOBP_DNA_REPLICATION"))
if(any(lengths(pan_sets)<20))stop("Unexpectedly small proliferation gene sets: ", paste(names(pan_sets)[lengths(pan_sets)<20], collapse=", "))
targets<-list(RL_Sig65=rl65,C5_stable_up=c5); all_sets<-c(targets,pan_sets)
set_tab<-data.table(GeneSet=names(pan_sets),N=lengths(pan_sets))

fs<-list.files(INPUT_DIR,pattern="\\.txt$",full.names=TRUE); coh<-combine_cohort(lapply(fs,read_expr))
cohort_order<-names(coh); cohort_order<-cohort_order[order(as.numeric(gsub("GSE","",cohort_order)))]
score_list<-list(); expr_rank_list<-list(); coverage_list<-list()

for(cc in cohort_order){
  obj<-coh[[cc]]; m<-obj$expression; if(!all(c("MKI67","PCNA")%in%rownames(m)))stop(cc, " is missing MKI67 or PCNA.")
  sc<-ssgsea(m,all_sets); sn<-intersect(names(obj$response),colnames(sc))
  pz<-t(apply(sc[names(pan_sets),sn,drop=FALSE],1,zscore)); pd<-as.data.table(t(pz)); setnames(pd,names(pan_sets))
  pd[,`:=`(Sample=sn,Cohort=cc,Response=obj$response[sn],y=ifelse(obj$response[sn]=="pCR",1L,0L),
           RL_Sig65=zscore(sc["RL_Sig65",sn]),C5_stable_up=zscore(sc["C5_stable_up",sn]),
           MKI67=zscore(m["MKI67",sn]),PCNA=zscore(m["PCNA",sn]))]; score_list[[cc]]<-pd
  av<-rowMeans(m,na.rm=TRUE)
  expr_rank_list[[cc]]<-data.table(Gene=rownames(m),Rank=rank(av,ties.method="average",na.last="keep")/sum(is.finite(av)))
  coverage_list[[cc]]<-rbindlist(lapply(names(pan_sets),function(nm)
    data.table(Cohort=cc,GeneSet=nm,Total=length(pan_sets[[nm]]),Matched=sum(pan_sets[[nm]]%in%rownames(m)))))
}
scores<-rbindlist(score_list,fill=TRUE); coverage<-rbindlist(coverage_list); pan_names<-names(pan_sets)
pc<-prcomp(as.matrix(scores[,..pan_names]),center=FALSE,scale.=FALSE); pc1<-pc$x[,1]; load<-pc$rotation[,1]
if(cor(pc1,rowMeans(scores[,..pan_names]),use="complete.obs")<0){pc1<--pc1; load<--load}
scores[,PanProlifPC1:=zscore(pc1)]; pc_var<-summary(pc)$importance[2,1]
load_tab<-data.table(Program=names(load),Loading=as.numeric(load)); load_tab[,Program:=factor(Program,levels=Program[order(Loading)])]
scores[,RL_resPan:=zscore(resid(lm(RL_Sig65~PanProlifPC1+MKI67))),by=Cohort]
scores[,C5_resPan:=zscore(resid(lm(C5_stable_up~PanProlifPC1+MKI67))),by=Cohort]

dr<-fread(DRIP_FILE); dr[,SYMBOL:=canon(SYMBOL)]; dr<-dr[!duplicated(SYMBOL)]
dr[,Strict:=toupper(as.character(Strict_E2_DRIP))%in%c("TRUE","T","1")]; univ<-dr$SYMBOL

gene_rho<-rbindlist(lapply(cohort_order,function(cc){
  m<-coh[[cc]]$expression; ss<-scores[Cohort==cc,.(Sample,PanProlifPC1)]
  pv<-setNames(ss$PanProlifPC1,ss$Sample); sn<-intersect(colnames(m),names(pv)); g<-intersect(univ,rownames(m))
  rr<-apply(m[g,sn,drop=FALSE],1,function(x)suppressWarnings(cor(x,pv[sn],method="spearman",use="complete.obs")))
  data.table(Gene=g,Cohort=cc,Rho=as.numeric(rr))
}),fill=TRUE)
gene_pan<-gene_rho[is.finite(Rho),.(PanProlifRho=median(Rho),Ncohort=.N),by=Gene][Ncohort>=3]

er<-Map(function(x,nm){y<-copy(x);setnames(y,"Rank",nm);y},expr_rank_list,names(expr_rank_list))
expr_ref<-Reduce(function(a,b)merge(a,b,by="Gene",all=TRUE),er); rc<-setdiff(names(expr_ref),"Gene")
expr_ref[,`:=`(ExprRank=rowMeans(.SD,na.rm=TRUE),NexprCohort=rowSums(!is.na(.SD))),.SDcols=rc]
expr_ref<-expr_ref[NexprCohort>=3,.(Gene,ExprRank,NexprCohort)]

burden<-data.table(Gene=univ)
for(nm in names(pan_sets))burden[,(nm):=as.integer(Gene%in%pan_sets[[nm]])]
burden[,ProlifBurden:=rowSums(.SD),.SDcols=names(pan_sets)]; burden<-burden[,.(Gene,ProlifBurden)]

##======================== 6. gene length + promoter GC ==========================##
txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
bs<-BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
gr<-GenomicFeatures::genes(txdb); gr<-gr[as.character(GenomicRanges::seqnames(gr))%in%names(bs)]
pr<-GenomicRanges::trim(GenomicRanges::promoters(gr,upstream=1000,downstream=1000))
gc<-rowSums(Biostrings::letterFrequency(BSgenome::getSeq(bs,pr),c("G","C"),as.prob=TRUE))
sy<-canon(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=names(gr),keytype="ENTREZID",column="SYMBOL",multiVals="first"))
anno<-data.table(Gene=sy,LogLength=log10(GenomicRanges::width(gr)),PromoterGC=gc)
anno<-anno[!is.na(Gene)&Gene!="",.(LogLength=median(LogLength),PromoterGC=median(PromoterGC)),by=Gene]

md<-data.table(Gene=univ,Strict=dr$Strict)
md<-Reduce(function(a,b)merge(a,b,by="Gene",all.x=TRUE),list(md,expr_ref,gene_pan,burden,anno))
md[,`:=`(C5=as.integer(Gene%in%c5),RL65=as.integer(Gene%in%rl65))]
md<-md[complete.cases(ExprRank,PanProlifRho,ProlifBurden,LogLength,PromoterGC)]

##======================== 7. gene-level logistic ==========================##
fit_target<-function(target,label,universe,stage,rhs){
  d<-copy(universe); d[,Target:=get(target)]
  q<-term_row(safe_fit(d,as.formula(paste0("Strict~Target+",rhs))),"Target"); if(is.null(q))return(data.table())
  q[,`:=`(Target=label,Stage=stage,N=nrow(d),N_target=sum(d$Target==1),
          DRIP_target=mean(d[Target==1]$Strict),DRIP_background=mean(d[Target==0]$Strict))]; q
}
stages<-list(
  "Expression + length + GC"="ExprRank+LogLength+PromoterGC",
  "+ continuous pan-proliferation"="ExprRank+LogLength+PromoterGC+PanProlifRho",
  "Full: + proliferation burden"="ExprRank+LogLength+PromoterGC+PanProlifRho+ProlifBurden",
  "Nonlinear pan-proliferation"="ExprRank+LogLength+PromoterGC+splines::ns(PanProlifRho,df=3)+ProlifBurden")
logit<-rbindlist(lapply(names(stages),function(st)
  rbind(fit_target("C5","C5 stable-up",md,st,stages[[st]]),
        fit_target("RL65","RL-Sig65",md,st,stages[[st]]))),fill=TRUE)
logit[,FDR:=p.adjust(P,"BH")]

##======================== 8. matched permutation ==========================##
match_perm<-function(target,label,d,covs,universe_label,B=10000,K=50){
  z<-copy(d[complete.cases(d[,..covs])]); z[,Target:=get(target)]
  for(v in covs)z[,(paste0("Z_",v)):=zscore(get(v))]
  tar<-z[Target==1,Gene]; pool0<-z[Target==0,Gene]
  near<-lapply(tar,function(g){
    q<-z[Gene==g]; pool<-z[Gene%in%pool0]; ds<-rep(0,nrow(pool))
    for(v in covs)ds<-ds+(pool[[paste0("Z_",v)]]-q[[paste0("Z_",v)]])^2
    pool$Gene[order(ds)][seq_len(min(K,nrow(pool)))]
  })
  keep<-lengths(near)>=10; tar<-tar[keep]; near<-near[keep]; obs<-mean(z[match(tar,Gene)]$Strict)
  set.seed(20260816+sum(utf8ToInt(label))+nrow(z))
  null<-replicate(B,{ctrl<-vapply(near,function(x)sample(x,1),character(1)); mean(z[match(ctrl,Gene)]$Strict)})
  ctrl1<-vapply(near,function(x)x[1],character(1))
  bal<-rbindlist(lapply(covs,function(v){
    a<-z[match(tar,Gene)][[v]]; b<-z[match(ctrl1,Gene)][[v]]; ps<-sqrt((var(a)+var(b))/2)
    data.table(Target=label,Universe=universe_label,Covariate=v,SMD=ifelse(is.finite(ps)&&ps>0,(mean(a)-mean(b))/ps,0))
  }))
  res<-data.table(Target=label,Universe=universe_label,N=length(tar),Observed=obs,Null_mean=mean(null),
                  Null_low=quantile(null,.025),Null_high=quantile(null,.975),Enrichment=obs/mean(null),
                  P=(1+sum(null>=obs))/(B+1),Max_abs_SMD=max(abs(bal$SMD)))
  list(res=res,balance=bal)
}
full_cov<-c("ExprRank","LogLength","PromoterGC","PanProlifRho","ProlifBurden")
m1<-match_perm("C5","C5 stable-up",md,full_cov,"All evaluable genes")
m2<-match_perm("RL65","RL-Sig65",md,full_cov,"All evaluable genes")

q75<-quantile(md$PanProlifRho,.75,na.rm=TRUE); high<-md[PanProlifRho>=q75]
high_logit<-rbind(
  fit_target("C5","C5 stable-up",high,"Top quartile PanProlifRho","ExprRank+LogLength+PromoterGC+PanProlifRho+ProlifBurden"),
  fit_target("RL65","RL-Sig65",high,"Top quartile PanProlifRho","ExprRank+LogLength+PromoterGC+PanProlifRho+ProlifBurden"),fill=TRUE)
high_logit[,FDR:=p.adjust(P,"BH")]
h1<-match_perm("C5","C5 stable-up",high,full_cov,"Top quartile PanProlifRho")
h2<-match_perm("RL65","RL-Sig65",high,full_cov,"Top quartile PanProlifRho")
match_res<-rbind(m1$res,m2$res,h1$res,h2$res); match_res[,FDR:=p.adjust(P,"BH")]
balance<-rbind(m1$balance,m2$balance,h1$balance,h2$balance)

nac<-rbind(meta_fit(scores,"RL_Sig65","RL-Sig65"),meta_fit(scores,"C5_stable_up","C5 stable-up"),
           meta_fit(scores,"PanProlifPC1","Pan-proliferation"),meta_fit(scores,"MKI67","MKI67"),
           meta_fit(scores,"PCNA","PCNA"),meta_fit(scores,"RL_resPan","RL residual"),
           meta_fit(scores,"C5_resPan","C5 residual"),fill=TRUE); nac[,FDR:=p.adjust(P,"BH")]

p1<-ggplot(load_tab,aes(Program,Loading))+geom_col(width=.7)+coord_flip()+
  labs(x=NULL,y="PC1 loading",title=paste0("Pan-proliferation PC1 · variance explained ",percent(pc_var,accuracy=.1)))+theme_pub(8.5)
savepdf("01_PanProliferation_PC1_loadings.pdf",p1,10,8)

nn<-nac[Effect%in%c("RL-Sig65","C5 stable-up","Pan-proliferation","MKI67","PCNA")&is.finite(OR)]
nn[,Effect:=factor(Effect,levels=rev(c("RL-Sig65","C5 stable-up","Pan-proliferation","MKI67","PCNA")))]
p2<-ggplot(nn,aes(Effect,OR,ymin=CI_low,ymax=CI_high,color=Effect))+
  geom_hline(yintercept=1,linetype=2,color="grey55")+geom_errorbar(width=.12,linewidth=.6)+geom_point(size=3)+
  scale_color_manual(values=PAL)+scale_y_log10()+coord_flip()+labs(x=NULL,y="Pooled OR for pCR per 1-SD increase")+
  theme_pub(9)+theme(legend.position="none")
savepdf("02_NAC_PanProliferation_reference.pdf",p2,10,7)

lp<-logit[Stage!="Nonlinear pan-proliferation"]; lp[,Stage:=factor(Stage,levels=names(stages)[1:3])]
p3<-ggplot(lp,aes(Stage,OR,ymin=CI_low,ymax=CI_high,color=Target,group=Target))+
  geom_hline(yintercept=1,linetype=2,color="grey55")+geom_line(linewidth=.45)+geom_errorbar(width=.08,linewidth=.6)+
  geom_point(size=2.8)+scale_y_log10()+scale_color_manual(values=c("RL-Sig65"="#b43665","C5 stable-up"="#187d79"))+
  labs(x=NULL,y="Adjusted OR for direct E2-induced DRIP",color=NULL)+theme_pub(8.5)+
  theme(axis.text.x=element_text(angle=25,hjust=1),legend.position="top")
savepdf("03_DRIP_pan_proliferation_multivariable_adjustment.pdf",p3,12,7)

match_res[,Universe:=factor(Universe,levels=c("All evaluable genes","Top quartile PanProlifRho"))]
p4<-ggplot(match_res,aes(Target,Observed,color=Target))+geom_errorbar(aes(ymin=Null_low,ymax=Null_high),
                                                                      width=.10,color="grey60",linewidth=.65)+geom_point(aes(y=Null_mean),shape=4,size=2.5,color="grey30")+
  geom_point(size=3)+facet_wrap(~Universe,nrow=1)+scale_color_manual(values=c("RL-Sig65"="#b43665","C5 stable-up"="#187d79"))+
  labs(x=NULL,y="Observed direct-DRIP fraction\n(null mean and 95% interval shown)",color=NULL)+theme_pub(8.5)+theme(legend.position="top")
savepdf("04_DRIP_pan_proliferation_matched_permutation.pdf",p4,13,6.5)

hl<-rbind(logit[Stage=="Full: + proliferation burden",.(Target,OR,CI_low,CI_high,P,FDR,Universe="All evaluable genes")],
          high_logit[,.(Target,OR,CI_low,CI_high,P,FDR,Universe="Top quartile PanProlifRho")])
hl[,Universe:=factor(Universe,levels=c("All evaluable genes","Top quartile PanProlifRho"))]
p5<-ggplot(hl,aes(Universe,OR,ymin=CI_low,ymax=CI_high,color=Target))+
  geom_hline(yintercept=1,linetype=2,color="grey55")+geom_errorbar(width=.10,position=position_dodge(.22),linewidth=.6)+
  geom_point(position=position_dodge(.22),size=3)+scale_y_log10()+
  scale_color_manual(values=c("RL-Sig65"="#b43665","C5 stable-up"="#187d79"))+
  labs(x=NULL,y="Adjusted OR for direct E2-induced DRIP",color=NULL)+theme_pub(8.5)+
  theme(axis.text.x=element_text(angle=18,hjust=1),legend.position="top")
savepdf("05_DRIP_high_proliferation_sensitivity.pdf",p5,10,6.5)

main<-(p1|p2)/(p3|p4)/(p5|plot_spacer())+plot_layout(heights=c(.95,1,1))
savepdf("06_PanProliferation_summary.pdf",main,16.6,20)

wb<-createWorkbook(); add<-function(n,x){addWorksheet(wb,n);writeData(wb,n,x)}
add("PanProlif_gene_sets",set_tab); add("Coverage",coverage); add("PC1_loadings",load_tab); add("NAC_scores",scores)
add("NAC_meta",nac); add("Gene_pan_rho",gene_rho); add("Gene_metadata",md); add("DRIP_logistic",logit)
add("DRIP_highProlif_logistic",high_logit); add("DRIP_matched",match_res); add("Match_balance",balance)
saveWorkbook(wb,file.path(OUT,"00_PanProliferation_results.xlsx"),overwrite=TRUE)

c5_full<-logit[Target=="C5 stable-up"&Stage=="Full: + proliferation burden"]
rl_full<-logit[Target=="RL-Sig65"&Stage=="Full: + proliferation burden"]
c5_high<-high_logit[Target=="C5 stable-up"]; rl_high<-high_logit[Target=="RL-Sig65"]
c5_match<-match_res[Target=="C5 stable-up"&Universe=="All evaluable genes"]
rl_match<-match_res[Target=="RL-Sig65"&Universe=="All evaluable genes"]
c5_hmatch<-match_res[Target=="C5 stable-up"&Universe=="Top quartile PanProlifRho"]
pass1<-nrow(c5_full)>0&&c5_full$OR>1&&c5_full$FDR<.05
pass2<-nrow(c5_match)>0&&c5_match$FDR<.05
pass3<-nrow(c5_high)>0&&c5_high$OR>1&&c5_high$FDR<.05
strong<-pass1&&pass2&&pass3
fmt<-function(x,d=3)ifelse(length(x)&&is.finite(x),format(signif(x,d),scientific=TRUE),"NA")

note<-c(
  "FINAL PAN-PROLIFERATION CEILING ANALYSIS","",
  paste0("Pan-proliferation programs: ",length(pan_sets),"; PC1 variance explained = ",round(100*pc_var,1),"%"),
  paste0("High-proliferation universe fixed at top quartile PanProlifRho; cutoff = ",round(q75,3)),"",
  paste0("C5 full logistic: OR=",round(c5_full$OR,3)," (",round(c5_full$CI_low,3),"-",round(c5_full$CI_high,3),
         "), FDR=",fmt(c5_full$FDR)," -> ",ifelse(pass1,"SUPPORTED","NOT SUPPORTED")),
  paste0("C5 full matched: observed=",round(c5_match$Observed,3),", null=",round(c5_match$Null_mean,3),
         ", FDR=",fmt(c5_match$FDR),", max|SMD|=",round(c5_match$Max_abs_SMD,3)," -> ",ifelse(pass2,"SUPPORTED","NOT SUPPORTED")),
  paste0("C5 high-proliferation logistic: OR=",round(c5_high$OR,3)," (",round(c5_high$CI_low,3),"-",
         round(c5_high$CI_high,3),"), FDR=",fmt(c5_high$FDR)," -> ",ifelse(pass3,"SUPPORTED","NOT SUPPORTED")),
  paste0("C5 high-proliferation matched sensitivity: observed=",round(c5_hmatch$Observed,3),", null=",
         round(c5_hmatch$Null_mean,3),", FDR=",fmt(c5_hmatch$FDR)),"",
  paste0("RL-Sig65 full logistic: OR=",round(rl_full$OR,3),", FDR=",fmt(rl_full$FDR)),
  paste0("RL-Sig65 full matched: observed=",round(rl_match$Observed,3),", null=",round(rl_match$Null_mean,3),
         ", FDR=",fmt(rl_match$FDR)),
  paste0("RL-Sig65 high-proliferation logistic: OR=",round(rl_high$OR,3),", FDR=",fmt(rl_high$FDR)),"",
  if(strong)
    "PRIMARY CONCLUSION: C5 excess direct E2-induced DRIP enrichment is not fully accounted for by broad canonical proliferation/mitotic transcription and remains detectable within a highly proliferation-associated gene background."
  else
    "PRIMARY CONCLUSION: The pre-specified three-layer criterion was not fully met. Do not claim C5 R-loop enrichment beyond broad proliferation.",
  "",
  "LANGUAGE BOUNDARY:",
  "Do NOT write 'independent of proliferation'.",
  "If all three primary layers pass, write that C5 enrichment was not fully accounted for by broad canonical proliferation/mitotic programs and persisted within a highly proliferation-associated gene background.",
  "No further proliferation signatures or threshold searches should be added after this analysis.")
writeLines(note,file.path(OUT,"00_PanProliferation_conclusion.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"99_sessionInfo.txt"))
