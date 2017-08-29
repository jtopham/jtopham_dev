## The purpose of this script is to perform all necessary pre-processing on the ROSMAP data files
## This includes removal of the first 10 principal components, followed by computation
## of residuals (via linear regression) based on four recorded confounding factors.
##
## Note that genotype data is excluded from PCA and residual computation
##
## Also note that this script is not meant to be executed in one-go. Rather, it needs to be 
## run carefully. Fortunately, only the genotype data should need to be altered (as
## new SNPs are incorporated from the imputed Broad data) 
##
## Finally, note that the script contains some code for analysis on the method of 
## PCA / residual computation, which can obviously be skipped (it can take quite a
## while to run). This code is located near the middle of the script


library(data.table)
library(ggplot2)
library(reshape2)

# Set up functions

## This function simply orders a matrix based on rownames

my.sorter <- function(matrix){
  return(matrix[order(rownames(matrix)),])
}

## This function removes columns of a dataframe that have more than 25% NA values
## Input is a dataframe

my.NAremover <- function(input.matrix){
  output.matrix <- input.matrix[, colSums(is.na(input.matrix)) <= round(0.25*nrow(input.matrix),digits=0)]
  return(output.matrix)
}


## This function removes the first x number of principal components from a dataframe
## Note that input must be a matrix with only data, without NAs or -Inf/Inf values


my.PCA <- function(input.matrix, num.pcs){
  A <- input.matrix %*% t(input.matrix)                         # calc covariance matrix
  E <- eigen(A,TRUE)                                            # get eigenvalues and vectors for A
  
  ggplot(data.frame(cbind(seq(1,10),E$values[1:10])),           # scree plot of first 10 PCs
         aes(x=X1, y=X2)) + geom_bar(stat="identity",
                                     color="black",
                                     fill="green",alpha=0.5)   
  
  FV <- as.matrix(E$vectors[,(num.pcs+1):ncol(E$vectors)])      # create feature vector; first x PCs removed
  
  tmp.matrix<- t(FV) %*% input.matrix                           # transform original matrix
  output.matrix<- FV %*% tmp.matrix
  rownames(output.matrix)<-rownames(input.matrix)
  return(output.matrix)
}
  
## This function returns residuals for given matrix using glmnet package
## Requires as input a matrix of values (rownames are sampleIDs) and
## a dataframe of confounding factors for which to use in model
## Outputs a matrix of same dimension. Requires glmnet package loaded

find.residuals <- function(matrix,con.factors){
  return(apply(matrix,2,FUN=function(x){
    model<-lm(dat ~ age_death + msex + EV1 + EV2 + EV3 + studyn, data=data.frame((cbind(con.factors,dat=x)),stringsAsFactors = F))
    return(model$residuals)
  }))
}

## This function serves to return a statistical value for
## a linear fit between individual predictors (ie. genes, probes)
## and a group of covariates (ie. age, gender, ancestry). It
## will return either a correlation coefficient for the fit
## for that predictor, or the pvalue for the entire model
## (calculated from f statistics)

get.stats <- function(matrix, con.factors,type){
  if(!(type=="rsq"||type=="pval")){return('Error: \'type\' parameter must be \'rsq\' or \'pval\'')}
  return(
    unname(apply(matrix,2,FUN=function(x){
    model<-lm(dat ~ age_death + msex + EV1 + EV2 + EV3 + studyn, data=data.frame((cbind(con.factors,dat=x)),stringsAsFactors = F))
    if(type=="rsq"){return(round(summary(model)$r.squared,digits=4))}
    if(type=="pval"){return(unname(pf(summary(model)$fstatistic[1],summary(model)$fstatistic[2],summary(model)$fstatistic[3],lower.tail=F)))}
})))
}


# Read-in dataframes, patient ID mapping files and confounding factor data

met.df<-fread("/Users/james/Desktop/jtopham_prj/data/rosmap/AMP-AD_ROSMAP_Rush-Broad_IlluminaHumanMethylation450_740_imputed.tsv", sep='\t',header=T,stringsAsFactors=F)
acet.df<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/chipSeqResiduals.csv", sep='\t',header=T,stringsAsFactors=F,check.names=F)
geno.df<-fread("/Users/james/Desktop/jtopham_prj/data/rosmap/rosmap_genotype_recode.tsv", sep='\t',header=T,stringsAsFactors=F)
exp.df<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/AMP_AD_ROSMAP_Broad-Rush_IlluminaHiSeq_RSEM_gene_FPKM_normalized_plate1-6.tsv", sep='\t',header=T,stringsAsFactors=F)

pt.mapper<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/AMP-AD_ROSMAP_Rush-Broad_IDKey.csv", sep=',',header=T,stringsAsFactors=F)
pt.mapper.exp<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/exp_ptmap.txt", sep='\t',header=F,stringsAsFactors=F)      # seperate mapping list for expression dataset needed
pt.mapper.met<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/p3.txt", sep='\t',header=T,stringsAsFactors=F)             # another mapping file needed for methylation dataset
pt.list<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/finalptlist.txt", sep='\t', header=F,stringsAsFactors=F)

c.factors<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/technical_factors_expression_ROSMAP.txt", sep='\t',header=T,stringsAsFactors=F)

# Housekeeping

geno.df<-my.NAremover(data.frame(geno.df))                                        # remove SNPs with many NA values
geno.data<-as.matrix(geno.df[,7:ncol(geno.df)])
rownames(geno.data)<-geno.df[,2]

exp.df<-exp.df[apply(exp.df[,3:ncol(exp.df)],1,FUN=function(x){sum(x>1)>=100}),]  # filter for genes that have FPKM > 1 in at least 100 samples

exp.df[,3:ncol(exp.df)]<-exp.df[,3:ncol(exp.df)]+1                                # add psuedocount to keep values == 0 zero when log transforming
exp.df<-cbind(exp.df[,1:2],apply(exp.df[,3:ncol(exp.df)],2,log2))                 # log transform using base of 2

colnames(exp.df)<-gsub('X','',colnames(exp.df))                                   # remove X in colnames that R adds in 
exp.df[,2]<-gsub('[.][0-9]*','',exp.df[,2])                           # remove decimal point after each ENSG gene ID

acet.df<-na.omit(acet.df)                                                         # some peaks have only NA values, remove them here

pt.mapper.exp<-cbind(gsub(":.*",'',pt.mapper.exp[,1]),gsub(".*:",'',pt.mapper.exp[,1]))
pt.mapper.met[,2]<-gsub('MAP[0]*','',pt.mapper.met[,2])
pt.mapper.met[,2]<-gsub('ROS[0]*','',pt.mapper.met[,2])

pt.mapper<-rbind(pt.mapper[,1:2],cbind(projid=pt.mapper[,1],gwas_id=pt.mapper[,3]),
                 cbind(projid=pt.mapper.exp[,1],gwas_id=pt.mapper.exp[,2]),
                 cbind(projid=pt.mapper.met[,2],gwas_id=pt.mapper.met[,1]))


# Filter dataframes for patients who were determined to have all 4 data types available


# Filter expression dataset
tmp.map<-pt.mapper[match(colnames(exp.df[,3:ncol(exp.df)]),pt.mapper[,2]),]             # find which exp samples are in mapping file
exp.df<-exp.df[,colnames(exp.df[,3:ncol(exp.df)])%in%c(tmp.map[,2],colnames(exp.df)[1:2])]   # filter for those in mapping file (33 missing)
colnames(exp.df)[3:ncol(exp.df)]<-na.omit(tmp.map[,1])                                  # map sampleIDs to project IDs
exp.df<-exp.df[,colnames(exp.df)%in%c(pt.list[,1],colnames(exp.df)[1:2])]               # filter for patients of interest

# Filter acetylation dataset
acet.df<-acet.df[,colnames(acet.df)%in%pt.list[,1]]

# Filter genotype dataset
rownames(geno.data)<-gsub('ROS','',rownames(geno.data))
rownames(geno.data)<-gsub('MAP','',rownames(geno.data))
geno.data<-geno.data[rownames(geno.data)%in%pt.list[,1],]

colnames(geno.data)<-gsub('_[A-Z]$','',colnames(geno.data))                           # convert SNP ID format to match mapping file (used in pipeline)
colnames(geno.data)<-gsub('[.]','-',colnames(geno.data))

# Filter methylation dataset
met.df<-data.frame(met.df)
tmp.map1<-gsub('[.]','-',colnames(met.df))
colnames(met.df)<-tmp.map1                                      # set colnames to parsed names
tmp.map2<-pt.mapper[match(tmp.map1,pt.mapper[,2]),]             # find which met samples are in mapping file (0 are missing)

colnames(met.df)[2:ncol(met.df)]<-na.omit(tmp.map2[,1])                                # map sampleIDs to project IDs
met.df<-met.df[,colnames(met.df)%in%c(pt.list[,1],colnames(met.df)[1])]                # filter for patients of interest


# Perform PCA to remove first 10 principal components
# Takes about 2 min each dataset. It will crash if any
# NA or -Inf/Inf values in data

# Expression
exp.data<-t(as.matrix(exp.df[,3:ncol(exp.df)]))
colnames(exp.data)<-exp.df[,2]

exp.PCA<-my.PCA(exp.data,10)

# Methylation
met.data<-t(as.matrix(met.df[,2:ncol(met.df)]))
colnames(met.data)<-met.df[,1]

met.PCA<-my.PCA(met.data,10)

# Acetylation
acet.data<-t(as.matrix(acet.df))

acet.PCA<-my.PCA(acet.data,10)

# PCA not applicable for genotype


# Sort all post-PCA matrices

exp.PCA<-my.sorter(exp.PCA)
acet.PCA<-my.sorter(acet.PCA)
met.PCA<-my.sorter(met.PCA)

# Also sort genotype matrix

geno.data<-my.sorter(geno.data)

# Find residuals to minimize effect of confounding factors


c.factors<-c.factors[c.factors[,1]%in%pt.list[,1],]
colnames(c.factors)[1]<-"projid"

c.factors.melt<-melt(c.factors[,c(1,7:12)], id.vars = c("projid"))
p<-ggplot(c.factors.melt,aes(x=value,fill=variable))
p+geom_density(alpha=0.5)+facet_wrap(~variable,scales="free")+
  theme(axis.text.x = element_text(size=13),axis.text.y=element_text(size=13),
        strip.text.x = element_text(size=22,face="bold"),axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 16, colour = "black"),legend.title=element_text(size=20))


c.factors2<-c.factors[,7:12]
rownames(c.factors2)<-c.factors[,1]

c.factors2<-my.sorter(c.factors2)                 # order to match ordering of data matrices



# check distribution of correlation coefficients for the different data types
# before and after PCA

all.rsqs<-data.frame(rbind(cbind(Dataset="Expression",PCA="Before",Correlation.Coefficient=get.stats(exp.data,c.factors2,"rsq")),
                           cbind(Dataset="Expression",PCA="After",Correlation.Coefficient=get.stats(exp.PCA,c.factors2,"rsq")),
                           cbind(Dataset="Acetylation",PCA="Before",Correlation.Coefficient=get.stats(acet.data,c.factors2,"rsq")),
                           cbind(Dataset="Acetylation",PCA="After",Correlation.Coefficient=get.stats(acet.PCA,c.factors2,"rsq")),
                           cbind(Dataset="Methylation",PCA="Before",Correlation.Coefficient=get.stats(met.data,c.factors2,"rsq")),
                           cbind(Dataset="Methylation",PCA="After",Correlation.Coefficient=get.stats(met.PCA,c.factors2,"rsq"))),stringsAsFactors=F)

all.rsqs$Correlation.Coefficient<-as.numeric(all.rsqs$Correlation.Coefficient)

write.table(all.rsqs,file="/Users/james/Desktop/jtopham_prj/data/rosmap/cc_dist_pca.tsv",sep='\t',quote=F,row.names=F)

p<-ggplot(all.rsqs,aes(x=Correlation.Coefficient,fill=PCA))
p+geom_density(alpha=0.3)+xlim(c(0,0.15))+xlab("Correlation Coefficient")+
  theme(axis.text.x = element_text(size=13),axis.text.y=element_text(size=13),
        strip.text.x = element_text(size=22,face="bold"),axis.title=element_text(size=16,face="bold"), plot.title = element_text(lineheight=.8, face="bold",size=21))+
  ggtitle("Distribution of Correlation Coefficients Before and After PCA")+
  facet_wrap(~Dataset)



## and also for pvalues, while also recording names of predictors this time


# check distribution of correlation coefficients for the different data types
# before and after PCA

all.pvals<-data.frame(rbind(cbind(Dataset="Expression",PCA="Before",feature=colnames(exp.data),p.Value=get.stats(exp.data,c.factors2,"pval")),
                           cbind(Dataset="Expression",PCA="After",feature=colnames(exp.data),p.Value=get.stats(exp.PCA,c.factors2,"pval")),
                           cbind(Dataset="Acetylation",PCA="Before",feature=colnames(acet.data),p.Value=get.stats(acet.data,c.factors2,"pval")),
                           cbind(Dataset="Acetylation",PCA="After",feature=colnames(acet.data),p.Value=get.stats(acet.PCA,c.factors2,"pval")),
                           cbind(Dataset="Methylation",PCA="Before",feature=colnames(met.data),p.Value=get.stats(met.data,c.factors2,"pval")),
                           cbind(Dataset="Methylation",PCA="After",feature=colnames(met.data),p.Value=get.stats(met.PCA,c.factors2,"pval"))),stringsAsFactors=F)

all.pvals2<-all.pvals
all.pvals2$p.Value<-as.numeric(all.pvals2$p.Value)

write.table(all.pvals2,file="/Users/james/Desktop/jtopham_prj/data/rosmap/pval_dist_pca.tsv",sep='\t',quote=F,row.names=F)

p<-ggplot(all.pvals2,aes(x=p.Value,fill=PCA))
p+geom_density(alpha=0.3)+xlim(c(0,1))+xlab("p Value")+
  theme(axis.text.x = element_text(size=13),axis.text.y=element_text(size=13),
        strip.text.x = element_text(size=22,face="bold"),axis.title=element_text(size=16,face="bold"), plot.title = element_text(lineheight=.8, face="bold",size=21))+
  ggtitle("Distribution of p Values Before and After PCA")+
  facet_wrap(~Dataset)


# Compare residuals computed using glmnet vs lm
# 
# test<-1081
# library(glmnet)
# model<-glmnet(as.matrix(c.factors2),exp.data[,test],alpha=0.5,family="gaussian")
# predictions<-predict(model,as.matrix(c.factors2),type="response")
# r<-exp.data[,test]-predictions[,which.min(colSums((exp.data[,test]-predictions)^2))]
# 
# model2<-lm(exp ~ age_death + msex + EV1 + EV2 + EV3 + studyn, data=data.frame((cbind(c.factors2,exp=exp.data[,test])),stringsAsFactors = F))
# r2<-model2$residuals
# 
# p<-ggplot(data.frame(cbind(glmnet=r,lm=r2)),aes(x=glmnet,y=lm))
# p+geom_point()+
#   theme(axis.text.x = element_text(size=13),axis.text.y=element_text(size=13),plot.title=element_text(lineheight=.8, face="bold",size=21),
#                      strip.text.x = element_text(size=22,face="bold"),axis.title=element_text(size=16,face="bold"),
#                      legend.text = element_text(size = 16, colour = "black"),legend.title=element_text(size=20))+
#   ggtitle("Comparison of Residuals Between 'glmnet' and 'lm' Packages")

colnames(exp.data)[test]
confounded.gene<-data.frame(cbind(Actual=exp.data[,test],
                                  Predicted=predictions[,which.min(colSums((exp.data[,test]-predictions)^2))]),
                                  c.factors2,stringsAsFactors=F)

p<-ggplot(confounded.gene,aes(x=Actual,y=Predicted,color=msex))
p+geom_point(size=2.5)+xlim(c(-7,10))+ylim(c(-7,10))+ggtitle("Predicted vs Actual Expression, Colored by Gender")+
  theme(axis.text.x = element_text(size=13),axis.text.y=element_text(size=13),plot.title=element_text(lineheight=.8, face="bold",size=21),
        strip.text.x = element_text(size=22,face="bold"),axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 16, colour = "black"),legend.title=element_text(size=20))+
  geom_abline(slope=1,intercept=0,color="black",linetype="longdash")+
  scale_color_gradientn(guide=guide_legend(title="Gender"),colours=c("Blue","Red"),breaks=c(0,1))


# calculate expression residuals (for all but genotype data) and save to file

exp.data.PCAresiduals<-find.residuals(exp.PCA,c.factors2)
write.table(exp.data.PCAresiduals,file="/Users/james/Desktop/jtopham_prj/data/rosmap/exp_data_residuals.tsv",sep='\t',quote=F,col.names=NA)

acet.data.PCAresiduals<-find.residuals(acet.PCA,c.factors2)
write.table(acet.data.PCAresiduals,file="/Users/james/Desktop/jtopham_prj/data/rosmap/acet_data_residuals.tsv",sep='\t',quote=F,col.names=NA)

met.data.PCAresiduals<-find.residuals(met.PCA,c.factors2)
write.table(met.data.PCAresiduals,file="/Users/james/Desktop/jtopham_prj/data/rosmap/met_data_residuals.tsv",sep='\t',quote=F,col.names=NA)


# add imputed SNPs of interest (from Broad) to genotype matrix

other.snps<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/parsed_AD_SNPs.tsv", header=T, sep='\t', stringsAsFactors=F)
other.snps2<-read.delim("/Users/james/Desktop/jtopham_prj/data/rosmap/parsed_imputed_SNPs.tsv", header=T, sep='\t', stringsAsFactors=F)
other.snps<-cbind(other.snps,other.snps2)
rownames(other.snps)<-other.snps[,1]
other.snps<-other.snps[,-1]
rownames(other.snps)<-gsub('ROS','',rownames(other.snps))
rownames(other.snps)<-gsub('MAP','',rownames(other.snps))
other.snps<-other.snps[rownames(other.snps)%in%pt.list[,1],]
other.snps<-my.sorter(other.snps)

other.snps<-other.snps[,!colnames(other.snps)%in%colnames(geno.data)]   # filter for those not already in local genotype matrix

geno.data2<-as.matrix(cbind(geno.data,other.snps))

colnames(geno.data2)<-gsub(".*__rs","rs",colnames(geno.data2))


# finally, write genotype matrix to file

write.table(geno.data2,file="/Users/james/Desktop/jtopham_prj/data/rosmap/geno_data_filtered.tsv",sep='\t',quote=F,col.names=NA)


# for saving files as RData objects, note that they were read in from the tsv file and then
# saved as binary. this means that the first column corresponds to (what were) the rownames

geno.data2<-cbind(sampleID=rownames(geno.data2),geno.data2)
rownames(geno.data2)<-NULL

for (i in 2:ncol(geno.data2)){
  geno.data2[,i]<-as.numeric(geno.data2[,i])
}

geno.data3<-fread("/Users/james/Desktop/jtopham_prj/data/rosmap/geno_data_filtered.tsv",sep='\t',header=T)
geno.data3<-data.frame(geno.data3, stringsAsFactors = F)

View(geno.data3[1:100,1:100])
class(geno.data3[,2])
geno.data<-geno.data3

save(geno.data,file="/Users/james/Desktop/jtopham_prj/data/rosmap/geno_data_filtered.RData")

load("/Users/james/Desktop/jtopham_prj/data/rosmap/geno_data_filtered.RData")
load("/Users/james/Desktop/jtopham_prj/data/rosmap/acet_data_residuals.RData")
