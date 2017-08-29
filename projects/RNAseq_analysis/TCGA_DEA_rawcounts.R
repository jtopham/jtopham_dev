#!/opt/R-3.2.0/bin/Rscript

## The purpose of this script is to:
##     run DESeq2 D.E.A. between KMT2D LOF/WT groups for each CA type with >=10 LOF mutations
##     note that samples with a non-LOF mutation in KMT2D are excluded
##     this script specifically runs DEA using RAW COUNTS for TCGA samples (http://bioinformatics.oxfordjournals.org/content/31/22/3666.long)
##     it is meant to be parallelized on marrahost01 using cmd line args

suppressWarnings(suppressMessages(library(DESeq2)))
#library("DESeq2", lib.loc="~/R/x86_64-unknown-linux-gnu-library/3.2")
#library("S4Vectors", lib.loc="~/R/x86_64-unknown-linux-gnu-library/3.2")

options(warn=-1)    # suppress warnings

args <- commandArgs(TRUE)                # use cmd line args
this.subset <- as.numeric(args[1])       # integer; 1 to 11

##------------------------------------------------------------------------------------------------------------
##         DATA
##------------------------------------------------------------------------------------------------------------

load('/projects/jtopham_prj/Thesis_prj/KMT2D_EZH2_paper/data/TCGA_raw_counts.RData')

# read in genotype data for all patients
load('/projects/jtopham_prj/Thesis_prj/data/R_objects/geno_clin_f.RData')

##------------------------------------------------------------------------------------------------------------
##         SCRIPT
##------------------------------------------------------------------------------------------------------------
colnames(TCGA.raw.counts)[1]<-"gene.id"

# note that 352 patients have at least one LOF KMT2D mutation
kmt2d.pts<-unique(geno.clin$cut.id[geno.clin$mapped.id=="KMT2D"&geno.clin$VARIANT_CLASSIFICATION%in%c("Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation")])
kmt2d.nonlof.pts<-unique(geno.clin$cut.id[geno.clin$mapped.id=="KMT2D"])
kmt2d.nonlof.pts<-kmt2d.nonlof.pts[!kmt2d.nonlof.pts%in%kmt2d.pts]

# determine which samples are tumor vs normal
pt.manifest<-data.frame(cbind(sample.id=colnames(TCGA.raw.counts)[-1],state=as.character(data.frame(do.call(rbind,strsplit(colnames(TCGA.raw.counts)[-1],"-")))[,4]),
                              cut.id=paste(data.frame(do.call(rbind,strsplit(colnames(TCGA.raw.counts)[-1],"-")))[,1],
                                           data.frame(do.call(rbind,strsplit(colnames(TCGA.raw.counts)[-1],"-")))[,2],
                                           data.frame(do.call(rbind,strsplit(colnames(TCGA.raw.counts)[-1],"-")))[,3],
                                           data.frame(do.call(rbind,strsplit(colnames(TCGA.raw.counts)[-1],"-")))[,4],sep='-')),stringsAsFactors = F)

# filter for only samples with genotype data available (6,850 samples)
TCGA.raw.counts<-TCGA.raw.counts[,c(1,(which(pt.manifest$cut.id%in%geno.clin$cut.id)+1))]

# transpose
tmp<-TCGA.raw.counts
TCGA.raw.counts<-data.frame(t(TCGA.raw.counts[,-1]),stringsAsFactors = F)
colnames(TCGA.raw.counts)<-tmp$gene.id
TCGA.raw.counts<-cbind(sample.id=pt.manifest$cut.id[pt.manifest$cut.id%in%geno.clin$cut.id],TCGA.raw.counts);rownames(TCGA.raw.counts)<-NULL
TCGA.raw.counts$sample.id<-as.character(TCGA.raw.counts$sample.id)

# annotate with CA type
TCGA.raw.counts<-cbind(type=geno.clin$cancer_type[match(TCGA.raw.counts$sample.id,geno.clin$cut.id)],TCGA.raw.counts)
TCGA.raw.counts$type<-as.character(TCGA.raw.counts$type)

# 11 cancer types have >=10 LOF mutations
ca.types<-names(which(table(TCGA.raw.counts$type[TCGA.raw.counts$sample.id%in%kmt2d.pts])>=10))

# 3,732 samples remaining, 285 of which have >=1 KMT2D LOF mutation
TCGA.raw.counts<-TCGA.raw.counts[TCGA.raw.counts$type%in%ca.types,]

# remove the 23 duplicated samples. 284 samples with MLL2 LOF now
TCGA.raw.counts<-TCGA.raw.counts[!duplicated(TCGA.raw.counts$sample.id),]

# remove samples with a non-LOF mutation in KMT2D. 3,415 samples now
TCGA.raw.counts<-TCGA.raw.counts[!TCGA.raw.counts$sample.id%in%kmt2d.nonlof.pts,]

# split DF into list of single-CA DFs
# and then remove cancer_type column
df.list = split(TCGA.raw.counts, f = TCGA.raw.counts$type)
names(df.list)<-unlist(lapply(df.list,FUN=function(x){return(unique(x[,1]))}))
df.list<-lapply(df.list,FUN=function(x){x<-x[,-1]})

# transpose each indiv. DF, set gene names as rownames
df.list<-lapply(df.list,FUN=function(x){
  tmp<-x$sample.id
  x<-data.frame(t(x[,-1]),stringsAsFactors = F)
  colnames(x)<-tmp
  for(i in 1:ncol(x)){x[,i]<-as.integer(x[,i])}
  return(x)
 }
)

#------------------------------------------------------------------------
# run DEA on each CA type
#------------------------------------------------------------------------

# construct experimental design object for each mini DF
colData<-lapply(df.list,FUN=function(x){
  y<-data.frame(cbind(condition=colnames(x)%in%kmt2d.pts,type=c(rep("single-read",ncol(x)))),stringsAsFactors = F)
  y$condition[y$condition==T]<-"KMT2D.LOF"
  y$condition[y$condition==F]<-"KMT2D.WT"
  y$condition<-factor(y$condition,levels=c("KMT2D.WT","KMT2D.LOF"))  # levels set in order of: untreated, treated - as per vignette
  rownames(y)<-colnames(x)
  return(y)
 }
)

# construct DESeq objects
DSq.Data<-list(0)
for (i in 1:length(df.list)){
  DSq.Data[[i]]<-DESeqDataSetFromMatrix(countData=df.list[[i]],colData=colData[[i]],design= ~condition)
}

# clear some big objects to save memory
rm(TCGA.raw.counts,geno.clin)

# perform DESeq2 on subset of data
DSq.Object<-DESeq(DSq.Data[[this.subset]])

# not extracting norm counts
#normalized.counts<-lapply(DSq.Object,FUN=function(x){return(as.data.frame(counts(x, normalized=TRUE)))})

# generate results tables
RNAseq.results<-data.frame(results(DSq.Object),stringsAsFactors = F)
RNAseq.results<-RNAseq.results[order(RNAseq.results$padj),]
RNAseq.results<-cbind(RNAseq.results,neg.log10.padj=-log10(RNAseq.results$padj))  # add column for -log10 tranformed adjusted p values
RNAseq.results<-cbind(RNAseq.results,abs.FC=abs(RNAseq.results$log2FoldChange))   # add column for absolute val fold change
RNAseq.results<-cbind(RNAseq.results,log2baseMean=log2(RNAseq.results$baseMean))  # add column for log2 baseMean
RNAseq.results<-cbind(RNAseq.results,pval.threshold=RNAseq.results$padj<0.05)     # add column for selective coloring
RNAseq.results<-cbind(type=ca.types[this.subset],RNAseq.results)                  # add column for cancer type

this.filename<-paste0('/projects/jtopham_prj/Thesis_prj/KMT2D_EZH2_paper/results/TCGA_DEA_res_rawcounts/',ca.types[this.subset],'_DEAres.RData')
save(RNAseq.results,file=this.filename)

print(paste0('Done ', ca.types[this.subset], ' analysis!'))                       # print to stdout
