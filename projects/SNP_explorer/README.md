## xQTLs and SNP Explorer

> Quantitative trait loci [(QTLs)](https://en.wikipedia.org/wiki/Quantitative_trait_locus) are discrete genomic loci that correlate with some phenotypic trait, with typical examples including expression of a particular gene or 5-hMC of a particular locus. Identification of disease QTLs, or disease-related SNPs that have some measurable effect on a trait, are a popular research topic, especially in the field of [GWAS](https://en.wikipedia.org/wiki/Genome-wide_association_study).

I completed this project during the same bioinformatics rotation in which I worked on the [DNA-methylation-Age project](https://github.com/jtopham/jtopham_dev/tree/master/projects/DNAmethylation_age), since it was during that time that I had access to the very rich (and at the time new) [ROSMAP](https://www.synapse.org/#!Synapse:syn3219045) dataset, which included genotype (SNP-chip), DNA methylation (Illumina 450k), expression (RNA-seq) and H3K27ac (ChIP-seq) data.

This project marked my first exposure to constructing a high-throughput pipeline, which was run on the computing resources provided by the Broad Institute. I created the pipeline using a combination of python and R, with over-arching *control panels* in bash. The final tool is fully functional and includes error handling, therefore making it readily useable by others!

The purpose of the SNP explorer pipeline is to identify xQTLs (eQTLs, mQTLs and aQTLs) given a list of SNPs as well as the respective data matrices (RNA-seq, DNA methylation and ChIP-seq). For all SNPs with all three types of QTLs, a [Bayesian network model](https://en.wikipedia.org/wiki/Bayesian_network) will be constructed using the R package *bnlearn*.

The pipeline consists of a control panel and 3 modules:

### [Control panel [bash]](https://github.com/jtopham/jtopham_dev/blob/master/projects/SNP_explorer/snp_explorer_cp.sh)

> This script interactively parses user inputs for several parameters
> such as filenames, window size in which to search for QTLs, and
> significance thresholds to use

### [Module Zero [bash]](https://github.com/jtopham/jtopham_dev/blob/master/projects/SNP_explorer/snp_explorer_m0.sh)

> This script serves in facilitating a threading process, in which the
> input list of SNPs is split into 10 subsets which are processed in 
> parallel

### [Module One [python]](https://github.com/jtopham/jtopham_dev/blob/master/projects/SNP_explorer/snp_explorer_m1.sh)

> This script finds all types of QTLs within the user-specified 
> window, and outputs intermediate lists of QTLs to local files

### [Module Two [R]](https://github.com/jtopham/jtopham_dev/blob/master/projects/SNP_explorer/snp_explorer_m2.sh)

> Finally, this script collects all the discovered QTLs and, where
> possible, constructs a bayesian network (BN) model for each SNP. The 
> final output is a list of SNPs and their respective QTLs, the
> significance levels of each, and the score of the best BN model

<br>


<br>

***

![pipeline](https://github.com/jtopham/jtopham_dev/blob/master/projects/SNP_explorer/tutorial/SE_outline.png?raw=true "SNP Explorer pipeline overview")

