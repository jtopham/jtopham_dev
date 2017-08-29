#!/usr/bin/env Rscript


## This script is the second and final module of the SNP Explorer pipeline
##
## It takes as input:
##	1. SNP ID
##	2. filename for bayesNet results
##	3. filename for QTL results
##	4. adjusted p value threshold
##	5. disease name for this SNP
##
## It requires 3 files that have the names of nearby genes, methylation probes, and acetylation peaks,
## titled z_rs123_genes.txt , z_rs123_probes.txt and z_rs123_peaks.txt, for example.
##
## It will identify all QTLs (above given threshold value) for the given SNP and output to results file. 
## If the SNP has all 3 types of QTLs, it will attempt to construct a bayesian network model by 
## exhaustively searching through 27 possible structures

options(warn=-1) # suppress warnings

suppressMessages(library(bnlearn))


# set up functions

## this function takes as input a data matrix
## and removes any columns with 10 or more NA values
## and outputs the filtered data matrix


my.NAremover <- function(input.matrix){
  output.matrix <- input.matrix[, colSums(is.na(input.matrix)) <= 10]
  return(output.matrix)
}



## this function takes the following as input:
## 	1. input matrix, where first column is sampleIDs, rest are data for features (ie. probes, peaks etc)
##	2. input vector, to correlate against all features in matrix (ie. genotype vector for a SNP)
##	3. QTL.type, a string specifying the type of QTL investigated. must be one of: 'aQTL', 'eQTL' or 'mQTL'
##	4. disease, a string specifying the disease name
##	5. feature.list, the list of nearby features that is produced by module 1 (python)
##	6. file.name, user-specified name of QTL results file (NOT bayesNet results file)
##	7. threshold, user-specified threshold for adjusted p value at which to call significance
## the function writes information pertaining to any QTLs to the results tsv
## and returns a nx2 matrix with sample IDs and data for the most significant feature, if one exists
## (if not, it will return -1)

my.correlate <- function(input.matrix,input.vector, QTL.type, disease, feature.list,file.name,threshold){

	# if only one feature is nearby:
	if(ncol(input.matrix)==2){
		pVal <- cor.test(as.numeric(input.matrix[,2]),input.vector[,2],method="spearman")$p.value
		if (pVal > threshold){
			print(paste0('Adjusted p value: ',pVal, ' does not meet threshold value ', '(', threshold, ')'))
			return(-1)
    		}
		print(paste("Best (and only) gene/probe/peak is: ",names(pVal), " with adj. p-value of: ", pVal,sep=""))
	
		if(QTL.type=="eQTL"){
			line<-c(disease,colnames(input.vector)[2],QTL.type,feature.list[which(feature.list[,1]==colnames(input.matrix)[2]),2],pVal,1,pVal,feature.list[which(feature.list[,1]==colnames(input.matrix)[2]),ncol(feature.list)])	
		}
		if(QTL.type!="eQTL"){
			line<-c(disease,colnames(input.vector)[2],QTL.type,colnames(input.matrix)[2],pVal,1,pVal,feature.list[which(feature.list[,1]==colnames(input.matrix)[2]),ncol(feature.list)])
		}	
		line <- as.matrix(t(line))
		write.table(line,file=file.name,row.names=F,quote=F,col.names=F,sep='\t', append=T)
		return(input.matrix)
  	}


	# otherwise, analyze all features
	pVals <- apply(input.matrix[,2:ncol(input.matrix)],2,FUN=function(x){return(cor.test(as.numeric(x),input.vector[,2],method="spearman")$p.value)})        # correlate cols of matrix with vector of interest
	pVals <- sort(pVals)
	pVals.adj <- p.adjust(pVals, method="bonferroni")  
  

	# if best feature is not significant:
	if(pVals.adj[1]>threshold){
		print(paste0('Adjusted p value: ',pVals.adj[1], ' does not meet threshold value ', '(', threshold, ')'))
    		return(-1)
	}
  
	# otherwise, see if there are more than one QTL
	print(paste("Best gene/probe/peak is: ",names(pVals.adj[1]), " with adj. p-value of: ", pVals.adj[1],sep=""))
	sig.names <- names(pVals.adj[pVals.adj<=threshold])                              							            	  # grab names of features that meet significance threshold
  
 	# if there is only one QTL
	if (length(sig.names)==1){
		if(QTL.type=="eQTL"){
                        line<-c(disease,colnames(input.vector)[2],QTL.type,feature.list[which(feature.list[,1]==sig.names),2],pVals[1],nrow(feature.list),pVals.adj[1],feature.list[which(feature.list[,1]==sig.names),ncol(feature.list)])             
                }
                if(QTL.type!="eQTL"){
                        line<-c(disease,colnames(input.vector)[2],QTL.type,sig.names,pVals[1],nrow(feature.list),pVals.adj[1],feature.list[which(feature.list[,1]==sig.names),ncol(feature.list)])
                }
                line <- as.matrix(t(line))
                write.table(line,file=file.name,row.names=F,quote=F,col.names=F,sep='\t',append=T)
		output <- data.frame(cbind(input.matrix[,1],input.matrix[,match(sig.names,colnames(input.matrix))]),stringsAsFactors=F)
		colnames(output)<-c("sampleID",sig.names)
		return(output)
	}

	# otherwise, append info for all QTLs to file
	for (i in 1:length(sig.names)){
		if(QTL.type=="eQTL"){
                	line<-c(disease,colnames(input.vector)[2],QTL.type,feature.list[which(feature.list[,1]==sig.names[i]),2],pVals[i],nrow(feature.list),pVals.adj[i],feature.list[which(feature.list[,1]==sig.names[i]),ncol(feature.list)])  
                }
        	if(QTL.type!="eQTL"){
                	line<-c(disease,colnames(input.vector)[2],QTL.type,sig.names[i],pVals[i],nrow(feature.list),pVals.adj[i],feature.list[which(feature.list[,1]==sig.names[i]),ncol(feature.list)])
                }
		line <- as.matrix(t(line))
		write.table(line,file=file.name,row.names=F,quote=F,col.names=F,sep='\t',append=T)
	}

	# finally, output matrix of sampleIDs and data vector for most significant feature
	output <- data.frame(cbind(input.matrix[,1],input.matrix[,match(sig.names[1],colnames(input.matrix))]),stringsAsFactors=F)
	#output[,2]<-as.numeric(output[,2])
 	colnames(output)<-c("sampleID", sig.names[1])
  	return(output)
}




# get SNP name

args <- commandArgs(TRUE)
this.snp <- args[1]
results.filename <- args[2]
results.filename2<- args[3]
usr.threshold<-as.numeric(args[4])
disease.name<-args[5]

# static filenames and dynamic filenames (latter is dependent on SNP ID)

exp.data.filename <- "/Users/james/Desktop/jtopham_prj/data/rosmap/exp_data_residuals.RData"  # expression data
acet.data.filename <- "/Users/james/Desktop/jtopham_prj/data/rosmap/acet_data_residuals.RData"  # acetylation data
met.data.filename <- "/Users/james/Desktop/jtopham_prj/data/rosmap/met_data_residuals.RData"  # methylation data
geno.data.filename <- "/Users/james/Desktop/jtopham_prj/data/rosmap/geno_data_filtered.RData"  # genotype data

exp.list.filename <- paste0("/Users/james/Desktop/jtopham_prj/scripts/rosmap/", "z_", this.snp, "_genes.txt")
met.list.filename <- paste0("/Users/james/Desktop/jtopham_prj/scripts/rosmap/", "z_", this.snp, "_probes.txt")
acet.list.filename <- paste0("/Users/james/Desktop/jtopham_prj/scripts/rosmap/", "z_", this.snp, "_peaks.txt")

# retrieve expression matrix

load(exp.data.filename)
exp.feats<-read.delim(exp.list.filename, header=F, sep='\t', stringsAsFactors=F)

if (sum(colnames(exp.data)%in%exp.feats[,1])==0){
	print('Error: no nearby genes found in data matrix')
	stop('No nearby genes found in data matrix')
}
exp.matrix<-exp.data[,c(1,which(colnames(exp.data)%in%exp.feats[,1]))]

# retreive acetylation matrix

load(acet.data.filename)
acet.feats<-read.delim(acet.list.filename, header=F, sep='\t', stringsAsFactors=F)

if (sum(colnames(acet.data)%in%acet.feats[,1])==0){
	print('Error: no nearby peaks found in data matrix')
        stop('No nearby peaks found in data matrix')
}
acet.matrix<-acet.data[,c(1,which(colnames(acet.data)%in%acet.feats[,1]))]


# retrieve methylation matrix

load(met.data.filename)
met.feats<-read.delim(met.list.filename, header=F, sep='\t', stringsAsFactors=F)

if (sum(colnames(met.data)%in%met.feats[,1])==0){
	print('Error: no nearby probes found in data matrix')
        stop('No nearby probes found in data matrix')
}
met.matrix<-met.data[,c(1,which(colnames(met.data)%in%met.feats[,1]))]


# retrieve genotype matrix

load(geno.data.filename)

if (sum(colnames(geno.data)==this.snp)==0){
	print('Error: SNP not found in genotype matrix')
	stop('SNP not found in genotype matrix')
}
geno.matrix<-geno.data[,c(1,which(colnames(geno.data)==this.snp))]


# correlations are done below
# QTL results are appended to QTL results file, within function
# if there exists at least one significant QTL, info is parsed
# that will be appended to bayesNet results file - if one of
# each QTL exists
# switch is used to control whether script continues to bayesNet
# analysis, dependent on whether SNP correlates with each of 3 features

switch <- 0

# correlate SNPs with expression

bestExp.vector <- my.correlate(exp.matrix, geno.matrix, "eQTL",disease.name,exp.feats,results.filename2,usr.threshold)

# if one of genes meets p value threshold, grab info for bayesNet results tsv file

if (class(bestExp.vector) != "numeric"){
	exp.pv<-cor.test(bestExp.vector[,2],geno.matrix[,2],method="spearman")$p.value
	exp.pv.a<-exp.pv*nrow(exp.feats)
	exp.name<-colnames(bestExp.vector)[2]
	gene.name <- exp.feats[which(exp.feats[,1]==exp.name),2]
	gene.dist <- exp.feats[which(exp.feats[,1]==exp.name),3]
}
if (class(bestExp.vector) == "numeric"){ switch <- -1 }

# correlate methylation probes with best SNP

bestProbe.vector <- my.correlate(met.matrix, geno.matrix, "mQTL",disease.name,met.feats,results.filename2,usr.threshold)

if (class(bestProbe.vector) != "numeric"){
	met.pv<-cor.test(bestProbe.vector[,2],geno.matrix[,2],method="spearman")$p.value
	met.pv.a<-met.pv*nrow(met.feats)
	try(met.exp.pv<-cor.test(bestProbe.vector[,2],bestExp.vector[,2],method="spearman")$p.value,silent=T)  # catch error if no eQTL
	met.name<-colnames(bestProbe.vector)[2]
	met.dist<-met.feats[which(met.feats[,1]==met.name),2]
}

if (class(bestProbe.vector) == "numeric"){ switch <- -1 }

# correlate acetylation probes with SNP

bestPeak.vector <- my.correlate(acet.matrix, geno.matrix, "aQTL",disease.name,acet.feats,results.filename2,usr.threshold)

if (class(bestPeak.vector) != "numeric"){
	acet.pv<-cor.test(bestPeak.vector[,2],geno.matrix[,2],method="spearman")$p.value
	acet.pv.a<-acet.pv*nrow(acet.feats)
	try(acet.exp.pv<-cor.test(bestPeak.vector[,2],bestExp.vector[,2],method="spearman")$p.value,silent=T)  # catch error if no eQTL
	try(acet.met.pv<-cor.test(bestPeak.vector[,2],bestProbe.vector[,2],method="spearman")$p.value,silent=T) # catch error if no mQTL
	acet.name<-colnames(bestPeak.vector)[2]
	acet.dist<-acet.feats[which(acet.feats[,1]==acet.name),2]
}

if (class(bestPeak.vector) == "numeric"){ switch <- -1 }

# sanity check plot, can be switched on if needed

#pdf(file = "plot.pdf")
#plot(bestSNP_vector[,2],bestPeak.vector[,2])
#dev.off()

#pdf(file = "plot2.pdf")
#plot(bestPeak.vector[,2], bestProbe_vector[,2])
#dev.off()

# bayesNet

# do not run bayesNet if genotype does not correlate with at least one feature of each type (met, exp, acet)

if (switch == -1){
	print('Aborting: SNP does not correlate with any features of one or more feature types')
	stop('SNP does not correlate with any features of one or more feature types')
}

rats<-data.frame(cbind(geno.matrix[,2], bestExp.vector[,2], bestProbe.vector[,2], bestPeak.vector[,2]),stringsAsFactors = FALSE)
colnames(rats)<-c("Genotype","Expression","Methylation", "Acetylation")


# hold off on this while using imputed SNPs
#rats[,1]<-as.factor(rats[,1])  # make genotype a discrete variable

if (sum(is.na(rats[,1])) > ((nrow(rats))/2)){
	print('Error: genotype vector has more than 50% NA values')
	stop('Genotype vector has more than 50% NA values')
}

rats<-na.omit(rats)

people_with_snp <- sum(rats[,1]>=1)

# exhaustively create every possible graph model (27)

model1<-empty.graph(names(rats))   # G to E
arcs(model1)<-matrix(c("Genotype", "Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model2<-empty.graph(names(rats))   # G to M
arcs(model2)<-matrix(c("Genotype", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model3<-empty.graph(names(rats))   # G to A
arcs(model3)<-matrix(c("Genotype", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))


model4<-empty.graph(names(rats))   # G to E and E to M
arcs(model4)<-matrix(c("Genotype", "Expression", "Expression", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model5<-empty.graph(names(rats))   # G to E and E to A
arcs(model5)<-matrix(c("Genotype", "Expression", "Expression", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model6<-empty.graph(names(rats))   # G to E and E to M and E to A
arcs(model6)<-matrix(c("Genotype", "Expression", "Expression", "Methylation", "Expression", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model7<-empty.graph(names(rats))   # G to E and E to M and M to A
arcs(model7)<-matrix(c("Genotype", "Expression", "Expression", "Methylation", "Methylation", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

model8<-empty.graph(names(rats))   # G to M and M to E
arcs(model8)<-matrix(c("Genotype", "Methylation", "Methylation", "Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model9<-empty.graph(names(rats))   # G to M and M to A
arcs(model9)<-matrix(c("Genotype", "Methylation", "Methylation", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model10<-empty.graph(names(rats))   # G to M and M to E and M to A
arcs(model10)<-matrix(c("Genotype", "Methylation", "Methylation", "Expression", "Methylation", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model11<-empty.graph(names(rats))   # G to M and M to E and E to A
arcs(model11)<-matrix(c("Genotype", "Methylation", "Methylation", "Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model12<-empty.graph(names(rats))   # G to M and M to A and A to E
arcs(model12)<-matrix(c("Genotype", "Methylation", "Methylation", "Acetylation", "Acetylation", "Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

model13<-empty.graph(names(rats))   # G to A and A to M
arcs(model13)<-matrix(c("Genotype", "Acetylation", "Acetylation", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model14<-empty.graph(names(rats))   # G to A and A to E
arcs(model14)<-matrix(c("Genotype", "Acetylation", "Acetylation", "Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model15<-empty.graph(names(rats))   # G to A and A to M and M to E
arcs(model15)<-matrix(c("Genotype", "Acetylation", "Acetylation", "Methylation", "Methylation", "Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model16<-empty.graph(names(rats))   # G to A and A to E and E to M
arcs(model16)<-matrix(c("Genotype", "Acetylation", "Acetylation", "Expression", "Expression", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model17<-empty.graph(names(rats))   # G to A and A to E and A to M
arcs(model17)<-matrix(c("Genotype", "Acetylation", "Acetylation", "Expression", "Acetylation", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

model18<-empty.graph(names(rats))   # G to E and G to M
arcs(model18)<-matrix(c("Genotype", "Expression", "Genotype", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model19<-empty.graph(names(rats))   # G to E and G to M and E to A
arcs(model19)<-matrix(c("Genotype", "Expression", "Genotype", "Methylation", "Expression", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model20<-empty.graph(names(rats))   # G to E and G to M and E to M
arcs(model20)<-matrix(c("Genotype", "Expression", "Genotype", "Methylation", "Expression", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))


model21<-empty.graph(names(rats))   # G to E and G to A
arcs(model21)<-matrix(c("Genotype", "Expression", "Genotype", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model22<-empty.graph(names(rats))   # G to E and G to A and A to M
arcs(model22)<-matrix(c("Genotype", "Expression", "Genotype", "Acetylation", "Acetylation", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model23<-empty.graph(names(rats))   # G to E and G to A and E to M
arcs(model23)<-matrix(c("Genotype", "Expression", "Genotype", "Acetylation", "Expression", "Methylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

model24<-empty.graph(names(rats))   # G to M and G to A
arcs(model24)<-matrix(c("Genotype", "Methylation", "Genotype", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model25<-empty.graph(names(rats))   # G to M and G to A and A to E
arcs(model25)<-matrix(c("Genotype", "Methylation", "Genotype", "Acetylation","Acetylation","Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))
model26<-empty.graph(names(rats))   # G to M and G to A and M to E
arcs(model26)<-matrix(c("Genotype", "Methylation", "Genotype", "Acetylation","Methylation","Expression"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))

model27<-empty.graph(names(rats))   # G to E and G to M and G to A
arcs(model27)<-matrix(c("Genotype", "Expression", "Genotype", "Methylation", "Genotype", "Acetylation"),ncol = 2, byrow = TRUE,dimnames = list(NULL, c("from", "to")))


model.list<-list(model1,model2,model3,model4,model5,model6,model7,model8,model9,model10,model11,model12,model13,model14,model15,model16,model17,model18,model19,model20,model21,model22,model23,model24,model25,model26,model27) # combine all possible models in a list

best.model<-model.list[which.max(unlist((lapply(model.list,FUN = function(x){score(x,rats)}))))][[1]]    # find top scoring model given the data

# append several params, header can be found in control panel (bash script)
line <- c(disease.name,people_with_snp,this.snp,
	  exp.name,gene.name,gene.dist,exp.pv,nrow(exp.feats),exp.pv.a,
          met.name,met.dist,met.pv,nrow(met.feats),met.pv.a,
          acet.name,acet.dist,acet.pv,nrow(acet.feats),acet.pv.a,
	  as.character(score(best.model,rats)),as.character(which.max(unlist((lapply(model.list,FUN = function(x){score(x,rats)})))), sep='    '),met.exp.pv, acet.exp.pv, acet.met.pv)

line <- as.matrix(t(line))

# write results to tsv file, under the user-supplied filename
write.table(line, file=results.filename, row.names = FALSE, quote=F, col.names = FALSE, sep='\t',  append=TRUE)


