## Gene set enrichment analysis tool (R)

> NGS analyses often produce, at some point, a list of genes of interest. Common examples include genes found to be differentially expressed in some condition, genes bound by a transcription factor of interest, or genes adjacent to enhancer regions of interest. In any case, many analyses avoid focusing on these genes individually, and instead look to infer something about them in the context of a biological pathway/system. 

Gene set enrichment analysis (GSEA) is a common approach to the above statement, as it tests whether known sets of related genes (ie. genes involved in cellular metabolism) are statistically enriched in the list of genes of interest. GSEA becomes particularly powerful in an exploratory-type analysis when a comprehensive list of gene sets are considered, such as those provided by the [Molecular Signatures Database](http://software.broadinstitute.org/gsea/msigdb) (mSigDB), which houses over 18,000 gene sets.

GSEA can be performed using several statistical frameworks. Consider a **sample**, of size *n*, drawn from a **population** of size *N*. If the population has *K* occurrences of interest, what is the chance that the sample contains *k* of these occurrences? The answer to this question can be determined for different values of *n*, *N*, *k* and *K* using the hypergeometric test (implemented in base R), which is derived from the [hypergeometric distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution).

Here, I've implemented hypergeometric tests to determine the level of significance for multiple gene sets individually, in terms of how enriched each gene set is given a list of genes of interest. Consider the mSigDB gene set that corresponds to [canonical genes involved in hypoxia](http://software.broadinstitute.org/gsea/msigdb/cards/HALLMARK_HYPOXIA.html). In this case, our **population** is all genes in the genome (*N* is roughly 20,000, then), and *K* is how many genes are in this gene set (200). If our sample corresponds to genes found to be differentially expressed (say, *n* = 150 genes), perhaps 20 of the genes in our sample are found in the hypoxia-related gene set. The hypergeometric test would provide a level of significance for this observation, in the form of a raw p value.

P-values are corrected using the [Benjamini-Hochberg procedure](https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini.E2.80.93Hochberg_procedure).

<br>

[GSEA tool [R]](https://github.com/jtopham/jtopham_dev/blob/master/utilities/GSEA/my_gsea.R)

> This will perform the hypergeometric-based GSEA on a list of genes of interest using user-supplied gene sets,
and can be invoked from within R. The script is reasonably fast, taking approx 30 seconds to test approximately
18,000 gene sets.
