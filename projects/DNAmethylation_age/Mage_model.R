library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(wateRmelon)
library(glmnet)
library(caret)

clin.files <- list.files("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/data", full.names=TRUE)
clin.files <- paste(clin.files, "/Clinical/Biotab",sep="")

# Explore which features are most commonly available for the 22 datasets

for (i in 1:length(clin.files)){
  this.file<-paste(clin.files[i],grep("clinical_patient",list.files(clin.files[i]),value = T),sep="/")
  tmp.file<-read.delim(this.file,stringsAsFactors=FALSE, header=TRUE)
  deats.tmp<-as.numeric(apply(tmp.file,2,FUN=function(x){table(x=="[Not Available]")["TRUE"]}))   # assess what proportion of each feature is [Not Available]
  colcount.tmp<-nrow(tmp.file)
  deats.tmp[is.na(deats.tmp)]<-0
  if(i==1){
    feats<-colnames(tmp.file)
    deats<-deats.tmp
    colcounts<-colcount.tmp}
  else{
    feats<-c(feats,(colnames(tmp.file)))
    deats<-c(deats,deats.tmp)
    colcounts<-c(colcounts,colcount.tmp)}
}

tmp.df<-data.frame(cbind(feats,deats),stringsAsFactors = F)
tmp.df[,2]<-as.numeric(tmp.df[,2])

tmp.df2<-aggregate(deats ~ feats, data = tmp.df, sum)    # collapse data frame, summing number of NAs for each feature across the 22 data sets

feats2<-table(feats)
col1<-names(feats2)

for(i in 1:length(feats2)){   # this just grabs the values in the table, as I don't know the shortcut R syntax
  if(i==1){
    col2<-feats2[[i]][1]}
  else{
    col2<-c(col2,feats2[[i]][1])}
}

feature.summary<-data.frame(cbind(feature = col1, count = col2),stringsAsFactors = F)
feature.summary<-cbind(feature.summary,NA_prcntg=tmp.df2[,2])
feature.summary[,2]<-as.numeric(feature.summary[,2])
feature.summary[,3]<-as.numeric(feature.summary[,3])
feature.summary[,3]<-(feature.summary[,3]/sum(colcounts))*feature.summary[,2]  # third col is now proportion of 2nd col that contained NA vals!
feature.summary<-feature.summary[order(feature.summary[,2],decreasing = T),]   # sort df
feature.summary[,1]<-factor(feature.summary[,1],levels=feature.summary[,1])
feature.summary2<-subset(feature.summary,count>15)

jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
paletteSize <- 256
jBuPuPalette <- jBuPuFun(paletteSize)

p <- ggplot(feature.summary2, aes(x=feature,y=count,fill=NA_prcntg)) 
p+geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=13),axis.text.y=element_text(size=13),
        legend.text = element_text(size = 14, colour = "black"),legend.title=element_text(size=13),
        axis.title=element_text(size=14,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=17))+
  xlab("")+ylab("Number of Datasets")+
  ggtitle("Clinical Features Available from TCGA Datasets")+
  scale_fill_gradient2(guide=guide_legend(title="Number of NA vals"),midpoint=10,low = jBuPuPalette[0], high = jBuPuPalette[250])+ylim(c(0,22))


# Now construct clinical feature dataframe based off of select features

my.feature.list<-c(grep("history_other_malignancy",as.character(feature.summary[feature.summary$count==22,1]),value=T,invert=T),"birth_days_to")   # features of interest. note that birth_days_to is missing from Pheochromocytoma dataset, 
                                                                                                                                                                          # and history_other_malignancy appears twice in colorectal dataset
clin.files<-grep("edddb5a5-3382-4e37-a662-1d28ea2a4f57",clin.files,invert=T,value=T)  # remove Pheochromocytoma dataset, which only has 3 samples and does not have age recorded in clinical files

clin.data <- NULL

for (i in 1:length(clin.files)){
  this.file <- paste(clin.files[i],grep("clinical_patient",list.files(clin.files[i]),value = T),sep="/")
  that.file<-paste(clin.files[i],grep(".org_biospecimen_cqcf",list.files(clin.files[i]),value=T),sep="/")
  this.type<-fread(that.file,select = "histological_type", header=T)
  tmp.file<-fread(this.file,select = my.feature.list,header=T)
  tmp.file2<-cbind(tmp.file[3:nrow(tmp.file),],this.type[-1,])
  if (i==1){
    clin.data <- as.data.frame(tmp.file2)}
  else{
    clin.data <- rbind(clin.data,as.data.frame(tmp.file2))}
}

# Check age distribution of these 838 cancer patients

clin.sbst<-subset(clin.data,clin.data$birth_days_to!="[Not Available]")
clin.sbst$birth_days_to<-round(as.numeric(clin.sbst$birth_days_to)/-360,digits = 0)
clin.sbst<-clin.sbst[order(clin.sbst$birth_days_to,decreasing = F),]
clin.sbst$bcr_patient_barcode<-factor(clin.sbst$bcr_patient_barcode,levels=clin.sbst$bcr_patient_barcode)

p<-ggplot(clin.sbst,aes(x=birth_days_to))
p+geom_density(fill=jBuPuPalette[150],alpha=0.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=13),axis.text.y=element_text(size=13),
        axis.title=element_text(size=17,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=19))+
  xlab("Age")+ylab("Density")+
  ggtitle("Age Distribution of TCGA Samples")+
  scale_x_continuous(breaks = seq(min(clin.sbst$birth_days_to), max(clin.sbst$birth_days_to), by = 5))


# And cancer type

clin.sbst<-subset(clin.data,clin.data$histological_type!="[Not Available]")

p<-ggplot(clin.sbst,aes(x=histological_type))
p+geom_bar(fill=jBuPuPalette[250],alpha=0.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=13),axis.text.y=element_text(size=13),
        axis.title=element_text(size=17,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=19))+
  xlab("")+ylab("Count")+
  ggtitle("Cancer Type Distribution of TCGA Samples")+
  scale_y_continuous(breaks = seq(0,nrow(clin.data), by = 25))

# Construct methylation data matrix for these patients

# if previously constructed:

if  (file.exists("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/meth_data.txt")){                                 
    
  meth.data<-read.delim("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/meth_data.txt", header=T, quote = "",check.names=F)
}

# otherwise, construct df

if  (!file.exists("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/meth_data.txt")){                              

# organize list of data files to iterate through

meth.files <- list.files("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/data", full.names=TRUE)
meth.files<-grep("edddb5a5-3382-4e37-a662-1d28ea2a4f57",meth.files,invert=T,value=T)  # again, remove Pheochromocytoma dataset
meth.files <- paste(meth.files, "DNA_Methylation",sep="/")

#for(i in 1:length(meth.files)){print(list.files(meth.files[i]))}  # some files have more than 1 array

z<-1
meth.files2<-NULL
for(i in 1:length(meth.files)){
  for(j in 1:length(list.files(meth.files[i]))){
    meth.files2[z]<-paste(meth.files[i],list.files(meth.files[i])[j],sep="/")
    z<-z+1
  }
}

meth.files2<-paste(meth.files2,"Level_2",sep="/")

# import probe list to filter for (these are the 25k probes I found to be assayed in both #27 and #450 arrays)

my.probelist<-read.delim("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/my_probe_list.txt",header=T)
my.probelist[,1]<-as.character(my.probelist[,1])

#construct df

meth.colnames<-"probeID"
meth.data<-NULL
for (i in 1:length(meth.files2)){
    for (j in 1:length(list.files(meth.files2[i]))){
      this.file <- paste(meth.files2[i],list.files(meth.files2[i])[j],sep="/")
      this.pt<-strsplit(list.files(meth.files2[i])[j],"-")
      this.pt2<-paste(gsub("[0-9].","",this.pt[[1]][3]),this.pt[[1]][4],this.pt[[1]][5],sep="-") 
      tmp.file<-data.frame(fread(this.file,header=T, select = c(1,2,3)),stringsAsFactors = F)[-1,]                  # column 2 is methylated val, column 3 is unmethylated
      tmp.file<-tmp.file[tmp.file[,1]%in%my.probelist[,1],]
      tmp.file[,2]<-as.numeric(tmp.file[,2])
      tmp.file[,3]<-as.numeric(tmp.file[,3])
      tmp.file2<-cbind(tmp.file[,1],apply(tmp.file,1,FUN=function(x){round(max(as.numeric(x[2]),0)/(max(as.numeric(x[2]),0) + max(as.numeric(x[3]),0) + 100),digits=5)}))        # calculate beta value
    if (i==1&&j==1){
      meth.data <- data.frame(cbind(tmp.file2[,1],as.numeric(tmp.file2[,2])),stringsAsFactors = F)}
    else{
      meth.data <- data.frame(cbind(meth.data,as.numeric(tmp.file2[,2])),stringsAsFactors = F)}
    meth.colnames<-c(meth.colnames,this.pt2)
  }
}

colnames(meth.data)<-meth.colnames    # this somehow is changing what View() produces (prior to setting colnames, every column past 2 is a duplicate of col 2)
meth.data[,2]<-as.numeric(meth.data[,2])  # second col is always character vector, bandage it here

write.table(meth.data,file="/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/meth_data.txt", col.names=T, row.names=F, sep='\t',quote = FALSE)

}

# Explore meth dataset

meth.data[,1]<-as.character(meth.data[,1])         # clean up some nuances from read.delim
meth.data[,1]<-gsub("\"","",meth.data[,1])
colnames(meth.data)<-gsub("\"","",colnames(meth.data))

# quantile normalize? check mean beta vals

meth.data2<-data.frame(t(meth.data),stringsAsFactors = F)
meth.data2<-data.frame(cbind(rownames(meth.data2),meth.data2),stringsAsFactors = F)
meth.data2[,1]<-as.character(meth.data2[,1]) 
colnames(meth.data2)<-meth.data2[1,]
meth.data2<-meth.data2[-1,]
rownames(meth.data2)<-NULL

# for some reason, TCGA-FZ-5919 is not in clinical feature df. remove

meth.data2<-meth.data2[-grep("TCGA-FZ-5919",meth.data2[,1]),]

# match ordering of both dfs

clin.data<-clin.data[order(clin.data[,2]),]
meth.data2<-meth.data2[order(meth.data2[,1]),]

# add 'type' variable from clin features df
# more stuff can be added here later!

meth.data2<-cbind(type=clin.data$histological_type,meth.data2)

# convert character cols back to num

for(i in 3:ncol(meth.data2)){                     
  meth.data2[,i]<-as.numeric(meth.data2[,i])}

#

rownames(meth.data2)<-NULL
colnames(meth.data2)[2]<-"sampleID"


# aggregate to find mean probe beta values across different tissue types

if(file.exists("/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data.txt")){
    agg.Data2<-read.delim("/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data.txt", sep='\t', header=T)
}

if(!file.exists("/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data.txt")){

agg.Data<-apply(meth.data2[,3:ncol(meth.data2)],2,FUN=function(x){aggregate(x,list(type=meth.data2$type),mean)})

agg.Data2<-NULL
for (i in 1:length(agg.Data)){
  if(i==1){agg.Data2<-agg.Data[[i]]}
  else{agg.Data2<-cbind(agg.Data2,agg.Data[[i]][2])}
}

colnames(agg.Data2)<-colnames(meth.data2[,-2])

write.table(agg.Data2,file="/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data.txt",sep='\t',quote=F,row.names=F)

}

agg.Data3<-melt(agg.Data2)

# plot beta val dist prior to normalization

p <- ggplot(agg.Data3, aes(x = value, colour = type))
p + geom_density()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=13),axis.text.y=element_text(size=13),
        legend.position="",
        axis.title=element_text(size=16,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=17))+
  ggtitle("Beta Value Distribution Prior to Normalization")


# quantile normalize using wateRmelon package's betaqn() function

met.data.norm<-cbind(meth.data2[,1:2],t(betaqn(t(as.matrix(meth.data2[,3:ncol(meth.data2)])))))

#write.table(met.data.norm,file="/Users/james/Desktop/jtopham_prj/data/normalized_meth_data_matrix.txt", sep='\t', quote=F,row.names=F)

# and aggregate again

if(file.exists("/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data_norm.txt")){
  agg.Data.norm2<-read.delim("/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data_norm.txt", sep='\t', header=T)
}

if(!file.exists("/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data_norm.txt")){
  
  agg.Data.norm<-apply(met.data.norm[,3:ncol(met.data.norm)],2,FUN=function(x){aggregate(x,list(type=met.data.norm$type),mean)})
  
  agg.Data.norm2<-NULL
  for (i in 1:length(agg.Data.norm)){
    if(i==1){agg.Data.norm2<-agg.Data.norm[[i]]}
    else{agg.Data.norm2<-cbind(agg.Data.norm2,agg.Data.norm[[i]][2])}
  }
  
  colnames(agg.Data.norm2)<-colnames(met.data.norm[,-2])
  
  write.table(agg.Data.norm2,file="/Users/james/Desktop/jtopham_prj/data/aggregated_TCGA_data_norm.txt",sep='\t',quote=F,row.names=F)
  
}

# and plot again, this time faceting

agg.Data4<-melt(agg.Data.norm2)
agg.Data4<-cbind(method="After Normalization",agg.Data4)
agg.Data3<-cbind(method="Before Normalization",agg.Data3)
agg.Data5<-rbind(agg.Data3, agg.Data4)

p <- ggplot(agg.Data5, aes(x = value, colour = type))
p + geom_density()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=13),axis.text.y=element_text(size=13),
        legend.position="",strip.text.x = element_text(size=12,face="bold"),
        axis.title=element_text(size=16,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=17))+
  ggtitle("") + facet_wrap(~method)


# normalize probe beta values across the different tissues, by subtracting mean probe
# values (across sample's tissue) from each sample

# tst <- data.frame(cbind(type=c("a","a","a","b","b","c"), sampleID=c("jim","bob","joe","sally","mary", "paul"),cg001=c(1,2,3,4,5,6), cg002=c(1,1,1,1,2,4)),stringsAsFactors = F)
# tst[,1]<-factor(tst[,1],levels=unique(tst[,1]))
# tst[,3]<-as.numeric(tst[,3])
# tst[,4]<-as.numeric(tst[,4])
# met.data.norm<-tst

tissue.means<-aggregate(.~type,data=met.data.norm[,-2],mean)

for (i in 1:length(levels(met.data.norm[,1]))){
  dims<-dim(met.data.norm[met.data.norm[,1]==tissue.means[i,1],3:ncol(met.data.norm)])
  met.data.norm[met.data.norm[,1]==tissue.means[i,1],3:ncol(met.data.norm)]<-round(met.data.norm[met.data.norm[,1]==tissue.means[i,1],3:ncol(met.data.norm)]-as.numeric(matrix(rep(tissue.means[i,2:ncol(tissue.means)],dims[1]),ncol=dims[2],byrow=T)),digits=5)
}


# Filter data set for those with age metadata available

clin.sbst<-subset(clin.data,clin.data$birth_days_to!="[Not Available]")
clin.sbst$birth_days_to<-round(as.numeric(clin.sbst$birth_days_to)/-360,digits = 0)
colnames(clin.sbst)[5]<-"age"

met.data.norm2<-met.data.norm[met.data.norm$sampleID%in%clin.sbst$bcr_patient_barcode,]  # filter for those with age data available
rownames(met.data.norm2)<-NULL



# as per Sara, only use data from tissues that have at least 10 samples
# note that I checked to ensure still good age distribution with this subset

good.types<-names(which(table(clin.sbst$histological_type)>=10))

clin.sbst<-clin.sbst[clin.sbst$histological_type%in%good.types,]

met.data.norm2<-met.data.norm2[met.data.norm2$type%in%good.types,]



# 10 fold CV - create folds that balance age variable

folds<-createFolds(clin.sbst$age,k=10)

means<-NULL     # sanity check to see that each fold has similar mean age 
for(i in 1:10){
  if(i==1){means<-mean(clin.sbst[folds[[i]],]$age)}
  else{means<-c(means,mean(clin.sbst[folds[[i]],]$age))}
}



## now perform parameter tuning to identify values of alpha and lambda
## that minimize MSE

# read in data if already ran (it takes ~30min to test the 100 alphas)

if  (file.exists("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/alpha_accuracies.txt")){                                 
  
  alpha_accuracies<-read.delim("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/alpha_accuracies.txt", header=T, quote = "",check.names=F)
}

# otherwise, construct df. this takes quite a while

if  (!file.exists("/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/alpha_accuracies.txt")){ 

set.seed(1991)
alpha_accuracies<-NULL
for (my.alpha in seq(0,1,(1/100))){
  for (i in 1:10){
    test.indeces<-as.numeric(unlist(folds[i]))
    train.indeces<-as.numeric(unlist(folds[-i]))
    model<-glmnet(as.matrix(met.data.norm2[train.indeces,3:ncol(met.data.norm2)]),clin.sbst$age[train.indeces],alpha=my.alpha,family="gaussian")
    predictions<-predict(model,as.matrix(met.data.norm2[test.indeces,3:ncol(met.data.norm2)]),type="response")
    res<-apply(predictions, 2, FUN=function(x){mean((x-clin.sbst$age[test.indeces])^2)})                                # mean-squared error for each of 100 lambda values for this fold
    alpha_accuracies<-rbind(alpha_accuracies,c(alpha=my.alpha,best.lambda=model$lambda[which.min(res)],MSE=min(res)))      # append results
  }
}  

alpha_accuracies<-data.frame(alpha_accuracies,stringsAsFactors = F)
write.table(alpha_accuracies,file="/Users/james/Desktop/jtopham_prj/data/methylation_pubdata/alpha_accuracies.txt", col.names=T, row.names=F, sep='\t',quote = FALSE)

}

alpha_accuracies2<-alpha_accuracies[alpha_accuracies$alpha>=0.2,]
rownames(alpha_accuracies2)<-NULL

agg.lambdalpha<-aggregate(cbind(best.lambda,MSE)~alpha, data=alpha_accuracies2, mean)


this.min<-as.numeric(summary(agg.lambdalpha$best.lambda)[1])
this.max<-as.numeric(summary(agg.lambdalpha$best.lambda)[6])
this.by<-(this.max-this.min)/7


p<-ggplot(agg.lambdalpha, aes(x=alpha,y=MSE,fill=best.lambda))
p+ geom_bar(stat="identity",alpha=0.7,color="black")+ coord_cartesian(ylim=c(50,65))+ scale_fill_gradientn(guide=guide_legend(title="Lambda"),colours=rainbow(5),breaks=seq(this.min,this.max,by=this.by))+
  ylab("CV Mean-squared Error")+xlab("Alpha")+ggtitle("Elastic Net Parameter Tuning")+
  theme(axis.text.x = element_text(size=14),axis.text.y=element_text(size=14),
  strip.text.x = element_text(size=14,face="bold"),legend.title=element_text(size=15),
  axis.title=element_text(size=18,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))



#my.std <- function(x) sd(x)/sqrt(length(x))
#A.errors<-apply(alpha_accuracies,2,my.std)
#limits<-aes(ymax = A.means[,2] + A.errors, ymin= A.means[,2] - A.errors)
#... + geom_errorbar(limits)


# best values for alpha/lambda (which minimizes MSE)
# note that these will be used when constructing final model

best.params<-agg.lambdalpha[which.min(agg.lambdalpha$MSE),]

best.alpha<-best.params$alpha
best.lambda<-best.params$best.lambda



## now run 10 fold CV (with nested 9 fold CV) to determine accuracy on training dataset
## for each round of the 10 fold CV, this will run 9 fold CV on the training folds 
## to determine best parameters to use on the test fold


set.seed(1991)
model.results<-NULL
for (i in 1:10){
  subset_alpha_accuracies<-NULL          # clear temp results df for parameter selection that will take place for each fold
  
  test.indeces<-as.numeric(unlist(folds[i]))
  train.indeces<-as.numeric(unlist(folds[-i]))
      
  subset.folds<-folds[-i]                # perform 9 fold CV on training folds to determine best parameters to use on test fold
  for (this.alpha in seq(0,1,(1/20))){   # test 20 different alphas, 100 lambdas
    for (j in 1:9){
      subset.test.indeces<-as.numeric(unlist(subset.folds[j]))
      subset.train.indeces<-as.numeric(unlist(subset.folds[-j]))
      model<-glmnet(as.matrix(met.data.norm2[subset.train.indeces,3:ncol(met.data.norm2)]),clin.sbst$age[subset.train.indeces],alpha=this.alpha,family="gaussian")
      predictions<-predict(model,as.matrix(met.data.norm2[subset.test.indeces,3:ncol(met.data.norm2)]),type="response") 
      res<-apply(predictions, 2, FUN=function(x){mean((x-clin.sbst$age[subset.test.indeces])^2)})                       # mean-squared error for each of 100 lambda values for this alpha for this fold
      subset_alpha_accuracies<-rbind(subset_alpha_accuracies,c(alpha=this.alpha,best.lambda=model$lambda[which.min(res)],MSE=min(res)))      # append results
  }
}

  this.agg.results<-aggregate(cbind(best.lambda,MSE)~alpha, data=subset_alpha_accuracies, mean)      # aggregate results to find average best lambda for each value of alpha
 
  this.best.params<-this.agg.results[which.min(this.agg.results$MSE),]

  this.best.alpha<-this.best.params$alpha
  this.best.lambda<-this.best.params$best.lambda

  print(paste0('Alpha: ',this.best.alpha))     # stdout sanity check
  print(paste0('Lambda: ',this.best.lambda)

  model<-glmnet(as.matrix(met.data.norm2[train.indeces,3:ncol(met.data.norm2)]),clin.sbst$age[train.indeces],alpha=this.best.alpha,family="gaussian")
  predictions<-predict(model,as.matrix(met.data.norm2[test.indeces,3:ncol(met.data.norm2)]),type="response",s=this.best.lambda)
  res<-cbind(gender=clin.sbst$gender[test.indeces],race=clin.sbst$race[test.indeces],type=clin.sbst$histological_type[test.indeces],actual=clin.sbst$age[test.indeces],predicted=predictions[,1])
  model.results<-rbind(model.results,res)  # actual vs predicted age
}

colnames(model.results)[1:3]<-c("Gender","Race","Type")    # capitalize for aesthetic plotting purposes

model.results<-data.frame(model.results,stringsAsFactors = F)
model.results[,1]<-factor(model.results[,1],levels=unique(model.results[,1]))
model.results[,2]<-factor(model.results[,2],levels=unique(model.results[,2]))
model.results[,3]<-factor(model.results[,3],levels=unique(model.results[,3]))
model.results[,4]<-as.numeric(model.results[,4])
model.results[,5]<-as.numeric(model.results[,5])


# for plotting purposes, include r^2 on plot

df<-data.frame(cbind(x=model.results$actual, y=model.results$predicted),stringsAsFactors = F)
lm.sum<-lm(y ~ x, df)
r.sq<-paste0("r^2 = ",round(summary(lm.sum)$r.squared,digits=4))

r.sq

# plot result, color by TYPE

p<-ggplot(model.results, aes(x=actual,y=predicted,color=Type))
p+geom_point()+xlab("Actual Age")+ylab("Predicted Age")+ggtitle(expression(paste("Actual vs. Predicted Age (Colored by Type), ", r^{2} ,' = 0.6144')))+
  theme(axis.text.x = element_text(size=14),axis.text.y=element_text(size=14),
  strip.text.x = element_text(size=14,face="bold"),legend.title=element_text(size=15),
        axis.title=element_text(size=18,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))+
  ylim(c(0,100))+xlim(c(0,100))+geom_smooth(aes(group = 1),method="lm",fullrange=T,fill=jBuPuPalette[200],alpha=0.2,color="black")


# plot result, color by GENDER

p<-ggplot(model.results, aes(x=actual,y=predicted,color=Gender))
p+geom_point()+xlab("Actual Age")+ylab("Predicted Age")+ggtitle(expression(paste("Actual vs. Predicted Age (Colored by Gender), ", r^{2} ,' = 0.6144')))+
  theme(axis.text.x = element_text(size=14),axis.text.y=element_text(size=14),
        strip.text.x = element_text(size=14,face="bold"),legend.title=element_text(size=15),
        axis.title=element_text(size=18,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))+
  ylim(c(0,100))+xlim(c(0,100))+geom_smooth(aes(group = 1),method="lm",fullrange=T,fill=jBuPuPalette[200],alpha=0.2,color="black")

# plot result, color by RACE

p<-ggplot(model.results, aes(x=actual,y=predicted,color=Race))
p+geom_point()+xlab("Actual Age")+ylab("Predicted Age")+ggtitle(expression(paste("Actual vs. Predicted Age (Colored by Race), ", r^{2} ,' = 0.6144')))+
  theme(axis.text.x = element_text(size=14),axis.text.y=element_text(size=14),
        strip.text.x = element_text(size=14,face="bold"),legend.title=element_text(size=15),
        axis.title=element_text(size=18,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))+
  ylim(c(0,100))+xlim(c(0,100))+geom_smooth(aes(group = 1),method="lm",fullrange=T,fill=jBuPuPalette[200],alpha=0.2,color="black")


## now show that model is predicting age rather than type by
## removing individual types from training and testing on that type
## note that this is now a 13-fold CV

new.folds<-lapply(as.list(good.types),FUN=function(x){which(met.data.norm2$type==x)})

set.seed(1991)
model.results2<-NULL
for (i in 1:length(new.folds)){
  test.indeces<-as.numeric(unlist(new.folds[i]))
  train.indeces<-as.numeric(unlist(new.folds[-i]))
  model<-glmnet(as.matrix(met.data.norm2[train.indeces,3:ncol(met.data.norm2)]),clin.sbst$age[train.indeces],alpha=best.alpha,family="gaussian")
  predictions<-predict(model,as.matrix(met.data.norm2[test.indeces,3:ncol(met.data.norm2)]),type="response",s=best.lambda)
  res<-cbind(type=clin.sbst$histological_type[test.indeces],actual=clin.sbst$age[test.indeces],predicted=predictions[,1])
  model.results2<-rbind(model.results2,res)  # actual vs predicted age
}

# take a quick peak at MSE across the different types

model.results<-cbind(method="Normal",model.results[,3:5],diff_sqrd=(as.numeric(model.results[,4])-as.numeric(model.results[,5]))^2)
model.results2<-cbind(method="Removal",model.results2,diff_sqrd=(as.numeric(model.results2[,2])-as.numeric(model.results2[,3]))^2)
model.results2<-data.frame(model.results2,stringsAsFactors = F)
model.results2[,1]<-factor(model.results2[,1],levels=unique(model.results2[,1]))
model.results2[,2]<-factor(model.results2[,2],levels=unique(model.results2[,2]))
model.results2[,3]<-as.numeric(model.results2[,3])
model.results2[,4]<-as.numeric(model.results2[,4])
model.results2[,5]<-as.numeric(model.results2[,5])

colnames(model.results)[2]<-"type"    # convert back to lower case

MSE1<-aggregate(diff_sqrd ~ type, data=model.results, mean)
MSE1<-cbind(method="Normal",MSE1)

MSE2<-aggregate(diff_sqrd ~ type, data=model.results2, mean)
MSE2<-cbind(method="Removal",MSE2)

MSE3<-rbind(MSE1,MSE2)

p<-ggplot(MSE3,aes(x=type,y=diff_sqrd,fill=method))
p+geom_bar(stat="identity",alpha=0.5,position="stack")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=13),axis.text.y=element_text(size=13),
        legend.text = element_text(size = 14, colour = "black"),legend.title=element_text(size=13),
        axis.title=element_text(size=17,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))+
  xlab("")+ylab("Mean-squared Error")+
  ggtitle("MSE Upon Removal of Different Types")


# now do scatter plots, faceting on each type being removed
# construct custom facet labels (to include r2 vals)

t<-tapply(model.results2$actual,model.results2$type,FUN=function(x){x})
t2<-tapply(model.results2$predicted,model.results2$type,FUN=function(x){x})

r.sqs<-NULL
for (i in 1:length(t)){
  lm.sum<-lm(t2[[i]] ~ t[[i]])
  r.sqs<-c(r.sqs,summary(lm.sum)$r.squared)
}

type.names<-names(t)

final.titles<-paste0(type.names,", r2 = ", round(r.sqs,digits=4))

model.results2$facet.labels<-apply(model.results2,1,FUN=function(x){final.titles[which(type.names==x[2])]})

p<-ggplot(model.results2, aes(x=actual,y=predicted,color=type))
p+geom_point()+xlab("Actual Age")+ylab("Predicted Age")+ggtitle("Actual vs. Predicted Age, Upon Individual Type Removal")+
  theme(axis.text.x = element_text(size=14),axis.text.y=element_text(size=14),
        strip.text.x = element_text(size=14,face="bold"),legend.position="",
        axis.title=element_text(size=18,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))+
  ylim(c(0,120))+xlim(c(0,120))+facet_wrap(~facet.labels)+stat_smooth(method="lm",fullrange = T,ncol=4)+ 
  theme(strip.text.x = element_text(size = 8,face="bold"))


# explore reason for Stomach and Kidney underperformance
# assess age distribution across different types

jBuPuFun2 <- colorRampPalette(brewer.pal(n = 9, "Blues"))
jBuPuPalette2 <- jBuPuFun2(paletteSize)

p<-ggplot(model.results2, aes(x=actual))
p+geom_density(fill=jBuPuPalette2[150],alpha=0.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=13),axis.text.y=element_text(size=13),
        axis.title=element_text(size=17,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=19))+
  xlab("Age")+ylab("Density")+
  ggtitle("Age Distribution Across Cancer Types")+
  scale_x_continuous(breaks = seq(min(clin.sbst$age), max(clin.sbst$age), by = 5))+
  facet_wrap(~type)+ 
  theme(strip.text.x = element_text(size = 12))


# plot stdev of age against MSE gain

MSE4<-MSE3[order(as.character(MSE3[,2])),]

MSE.diff<-MSE4[MSE4$method=="Removal",3]-MSE4[MSE4$method=="Normal",3]
std.age<-aggregate(actual~type,data=model.results2,sd)
std.age<-std.age[order(as.character(std.age[,1])),]

std.age<-cbind(std.age,MSE.diff)
colnames(std.age)[1:2]<-c("Type","stdev.age")

df<-data.frame(cbind(x=std.age$stdev.age, y=std.age$MSE.diff),stringsAsFactors = F)
lm.sum<-lm(y ~ x, df)
r.sq<-paste0("r^2 = ",round(summary(lm.sum)$r.squared,digits=4))


p<-ggplot(std.age,aes(x=stdev.age,y=MSE.diff,color=Type))
p+geom_point(size=5)+xlab("Standard Deviation of Age")+ylab("Increase in MSE")+ggtitle(expression(paste("Standard Deviation of Age vs. Increase in MSE, ", r^{2} ,' = 0.0161')))+
  theme(axis.text.x = element_text(size=14),axis.text.y=element_text(size=14),
        strip.text.x = element_text(size=14,face="bold"),legend.text = element_text(size = 14, colour = "black"),legend.title=element_text(size=13),
        axis.title=element_text(size=18,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))+
  geom_smooth(aes(group = 1),method="lm",fullrange=T,color="black",se=F)


# plot sample number vs MSE gain

nvsg<-data.frame(type=as.character(met.data.norm2$type),stringsAsFactors=F)
nvsg<-cbind(nvsg,count=1)
nvsg<-aggregate(count~type,data=nvsg,sum)
nvsg<-nvsg[order(nvsg[,1]),]
nvsg<-cbind(nvsg,MSE.diff=std.age$MSE.diff)
nvsg[,1]<-factor(nvsg[,1],levels=nvsg[,1])


df<-data.frame(cbind(x=nvsg$count, y=nvsg$MSE.diff),stringsAsFactors = F)
lm.sum<-lm(y ~ x, df)
r.sq<-paste0("r^2 = ",round(summary(lm.sum)$r.squared,digits=4))


p<-ggplot(nvsg,aes(x=count,y=MSE.diff,color=type))
p+geom_point(size=5)+xlab("Number of Samples")+ylab("Increase in MSE")+ggtitle(expression(paste("Number of Samples vs. Increase in MSE, ", r^{2} ,' = 0.0523')))+
  theme(axis.text.x = element_text(size=14),axis.text.y=element_text(size=14),
        strip.text.x = element_text(size=14,face="bold"),legend.text = element_text(size = 14, colour = "black"),legend.title=element_text(size=13),
        axis.title=element_text(size=18,face="bold"),plot.title = element_text(lineheight=.8, face="bold",size=20))+
  geom_smooth(aes(group = 1),method="lm",fullrange=T,color="black",se=F)


# heatmap
# 
# library(pheatmap)
# 
# hm.sbst<-t(met.data.norm2[,3:1500])
# colnames(hm.sbst)<-met.data.norm2[,2]
# 
# pheatmap(hm.sbst,annotation_col=data.frame(met.data.norm2[,1],row.names=colnames(hm.sbst)),legend=F,show_rownames=F,show_colnames=F)







