## DNA methylation and age

> "DNA methylation age measures the cumulative effect of an epigenetic maintenance system. This novel epigenetic clock can be used to address a host of questions in developmental biology, cancer and aging research" (Horvath 2013)

Associations between 5-hydroxymethylcytosine [(DNA methylation)](https://en.wikipedia.org/wiki/5-Hydroxymethylcytosine) and ageing were established in a number of different studies, which generally have pointed to a [link between chronological age and DNA hypermethylation](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-015-0118-4). In 2013, Steve Horvath leveraged a predictive machine learning framework to model chronological age based on DNA methylation levels (from Illumina probe-based assays) in [this paper](https://www.ncbi.nlm.nih.gov/pubmed/24138928).

Interestingly, once the model has been trained on healthy samples, it can be applied in disease settings to investigate whether changes in apparant age (coined *molecular age*) could be linked to different diseases such as neurodegenerative disorders and cancer.

During the third bioinformatics rotation during my grad studies, I had at my disposal the The Religious Orders Study and Memory and Aging Project (ROSMAP) Study [(ROS/MAP)](https://www.synapse.org/#!Synapse:syn3219045) dataset, which had Illumina 450k data for several hundred patients, along with clinical data which included age. I constructed my own machine learning system to model age based on probe intensities at a number of representative loci along the genome. This project encompassed a few really interesting techniques such as:

### [Normalization by linear regression (*residualization*)](https://github.com/jtopham/jtopham_dev/blob/master/projects/DNAmethylation_age/AD_DF_construction_PCA_residuals.R#L63-L68)

### [Nested cross-validation](https://github.com/jtopham/jtopham_dev/blob/master/projects/DNAmethylation_age/Mage_model.R#L423-L441)

<br>

-[Data pre-processing script [R]](https://github.com/jtopham/jtopham_dev/blob/master/projects/DNAmethylation_age/AD_DF_construction_PCA_residuals.R)

> The first 10 principal components are removed in a general attempt to remove confounding factors (such 
> as gender). For the methylation data, *residualization* is performed as an additional normalization 
> step. Briefly: a linear model (limma) is constructed to predict beta values at each probe given 
> information regarding known/recorded confounding factors (gender, ethnicity, study number etc.).
> This linear model is then applied to predict the beta value of every probe. Each prediction is
> therefore the amount of variation able to be explained by confounding factors alone, so the 
> residual (the remaining part that is not explained by confounding factors) is kept as the new
> (transformed) data point

<br>

-[Main analysis script [R]](https://github.com/jtopham/jtopham_dev/blob/master/projects/DNAmethylation_age/Mage_model.R)

> The model is built off matched-normal ("healthy") data from 838 cancer patients (TCGA project).
> Methylation data is quantile normalized using wateRmelon package's *betaqn()* function, and
> further normalized to account for tissue-specific differences by mean centering samples by
> tissue type. An [elastic net](https://en.wikipedia.org/wiki/Elastic_net_regularization) model is
> constructed, with parameter tuning performed to select appropriate values for alpha and lambda.
> This is done via nested cross-validation, in which a 9-fold CV is performed to select optimal
> parameters to be used within each fold of the outer 10-fold CV 


