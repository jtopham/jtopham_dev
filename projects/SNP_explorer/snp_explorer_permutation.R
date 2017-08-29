#!/usr/bin/env Rscript


## The purpose of this script is to validate the bayesian network
## analysis results using a permutation test
## It takes as input a SE_BN_results.tsv file
## and outputs a BN results file of the same format
##
## Example execution:
## Rscript snp_explorer_permutation.R /SE_BN_results.tsv

options(warn=-1) # suppress warnings

suppressMessages(library(bnlearn))


# write header to permutation results tsv file

results.filename<-"SE_BN_permutation_results.tsv"

header<-data.frame(t(c("SNP.ID","Ensembl.ID","Gene.Name","Geno.Exp.pvalue","Probe.ID","Geno.Met.pvalue",
          "Peak.Name","Geno.Acet.pvalue","Model.Score","Model.Structure")))

write.table(header,file=results.filename,row.names = FALSE, quote=F, col.names = FALSE, sep='\t')

# get name of SNP Explorer BN results file

args <- commandArgs(TRUE)
usr.filepath <- args[1]

#BN.res<-read.delim("SE_BN_concat_results.tsv",header=T,sep='\t',stringsAsFactors=F) # testing purposes

# retrieve SE BN results data

BN.res<-read.delim(usr.filepath, header=T, sep='\t', stringsAsFactors=F)

# retrieve expression matrix

load("/Users/james/Desktop/jtopham_prj/data/rosmap/exp_data_residuals.RData")

exp.matrix<-exp.data[,which(colnames(exp.data)%in%BN.res$Ensembl.ID)]

exp.matrix<-exp.matrix[,match(BN.res$Ensembl.ID,colnames(exp.matrix))] # preserve ordering

# retreive acetylation matrix

load("/Users/james/Desktop/jtopham_prj/data/rosmap/acet_data_residuals.RData")

acet.matrix<-acet.data[,which(colnames(acet.data)%in%BN.res$Peak.Name)]

acet.matrix<-acet.matrix[,match(BN.res$Peak.Name,colnames(acet.matrix))] # preserve ordering

# retrieve methylation matrix

load("/Users/james/Desktop/jtopham_prj/data/rosmap/met_data_residuals.RData")

met.matrix<-met.data[,which(colnames(met.data)%in%BN.res$Probe.ID)]

met.matrix<-met.matrix[,match(BN.res$Probe.ID,colnames(met.matrix))] # preserve ordering

# retrieve genotype matrix

load("/Users/james/Desktop/jtopham_prj/data/rosmap/geno_data_filtered.RData")

geno.matrix<-geno.data[,which(colnames(geno.data)%in%BN.res$SNP.ID)]

geno.matrix<-geno.matrix[,match(BN.res$SNP.ID,colnames(geno.matrix))] # preserve ordering

# perform permutation testing

for (i in 1:1000){
  geno.matrix<-geno.matrix[sample(c(1:nrow(geno.matrix))),] # permute rows of genotype matrix

  for (j in 1:ncol(geno.matrix)){
    rats<-data.frame(cbind(geno.matrix[,j], exp.matrix[,j], met.matrix[,j], acet.matrix[,j]),stringsAsFactors = FALSE)
    colnames(rats)<-c("Genotype","Expression","Methylation", "Acetylation")
    rats<-na.omit(rats)

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
    line <- c(colnames(geno.matrix)[j],
              colnames(exp.matrix)[j],
              BN.res$Gene.Name[j],
              cor.test(geno.matrix[,j],exp.matrix[,j],method="spearman")$p.value,
              colnames(met.matrix)[j],
              cor.test(geno.matrix[,j],met.matrix[,j],method="spearman")$p.value,
              colnames(acet.matrix)[j],
              cor.test(geno.matrix[,j],acet.matrix[,j],method="spearman")$p.value,
    	  as.character(score(best.model,rats)),as.character(which.max(unlist((lapply(model.list,FUN = function(x){score(x,rats)})))), sep='    '))
    
    line <- as.matrix(t(line))
    
    # write results to tsv file
    write.table(line, file=results.filename, row.names = FALSE, quote=F, col.names = FALSE, sep='\t',  append=TRUE)
  }
}
    
