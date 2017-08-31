## RNA-sequencing analysis

> Here I've included some scripts relevant in the analysis of RNA-sequencing data. I find analysis of this NGS data type the most trivial, and therefore don't have many scripts I think would prove useful to people. However, in the near future I'm planning to construct a complete RNA-seq analysis (DEA) pipeline in CWL/Makefile as a tool capable of automating the entire process for other folks.

### [Assigning read counts to genes [bash]](https://github.com/jtopham/jtopham_dev/blob/master/projects/RNAseq_analysis/featCount_script.sh)

> After processing RNA-seq BAM files, the first step in any
> analysis pipeline often consists of constructing a
> *count matrix* for every sample. In this step, read counts
> are assigned to each gene in the specified genome. For this
> I prefer using *featureCounts* from Subread (at least 
> compared to HTSeq) - as I've found that it picks up lowly-
> expressed and lncRNAs better more effectively

<br>


<br>

***


