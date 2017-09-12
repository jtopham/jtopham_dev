#!/usr/bin/env Rscript

## This script provides a homemade function to perform threshold-based GSEA (using phyper).
## It is employs multiple hypergeometric tests, corrected for multiple hypotheses (BH).

my.GSEA.hyper <- function(input.genes, map.file, all.genes){  
  ## This function takes as input:
  ## (1) a gene list (orderless; vector)
  ## (2) a mapping file (msigID to genes), ie 'all_msig.RData'
  ## (3) a list of all genes assayed in the experiment (ie all genes that were assigned
  ## a read count in RNA-seq data
  ## Output is dataframe of q-values for each gene set tested, along with other info.
  res <- NULL
  overlap.size <- NULL
  
  for (i in 1:nrow(map.file)){
    this.msig <- unlist(unname(map.file[i, 1]))
    
    # ensure genes in gene set were assayed for
    genes.in.set <- unname(unlist(map.file[i ,2]))[unname(unlist(map.file[i ,2])) %in% all.genes]
    
    # find how many genes overlap between genes in gene set and genes of interest
    overlap <- sum(input.genes %in% genes.in.set)
    
    if(overlap > 1){
      res[i] <- phyper(overlap - 1, length(input.genes), length(all.genes) - length(input.genes),
                       length(genes.in.set),lower.tail = F)
      overlap.size[i] <- overlap
    }
    if(overlap <= 1){
      res[i] <- 1
      overlap.size[i] <- overlap
    }

    # be verbose and indicate when halfway done
    if(i == floor(nrow(map.file)/2)){print('Halfway done!')}
  }
  res2 <- data.frame(cbind(msig.annotation = unlist(unname(map.file[, 1])) , p.val = res,
                           num.overlap = overlap.size, num.input = length(input.genes)), 
                     stringsAsFactors = F)
  res2$p.val <- as.numeric(res2$p.val)
  res2 <- cbind(res2, padj = p.adjust(res2$p.val, method="BH"))
  res2 <- res2[order(res2$padj), ]
  return(res2)
}
