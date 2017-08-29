#!/bin/bash

printf ' _____ _   _ ______   _______   ________ _     ___________ ___________
/  ___| \ | || ___ \ |  ___\ \ / /| ___ \ |   |  _  | ___ \  ___| ___ \
\ `--.|  \| || |_/ / | |__  \ V / | |_/ / |   | | | | |_/ / |__ | |_/ /
 `--. \ . ` ||  __/  |  __| /   \ |  __/| |   | | | |    /|  __||    /
/\__/ / |\  || |     | |___/ /^\ \| |   | |___\ \_/ / |\ \| |___| |\ \
\____/\_| \_/\_|     \____/\/   \/\_|   \_____/\___/\_| \_\____/\_| \_|


'

## this script acts as a control panel for the SNP Explorer pipeline

## it takes as argument a filepath to a tsv file that contains disease
## names (must be one word, ie. 'alzheimers_disease') in it's first column
## and SNP IDs (ie. rs12345) in it's second column. these two fields are
## to be separated by a tab, with each SNP ID followed by a newline '\n'

## if the input file consists of more than 100 lines, the file is split
## into 10 smaller subset files. each are then sent through the pipeline
## in parallel (essentially 10 threads on single CPU)

## the pipeline will output progress (and error handling) to STDOUT. the final output
## (two tsv files containing results on QTLs and successfully-constructed
## bayesian networks) files are named 'SE_QTL_results.tsv' and 'SE_BN_results.tsv'

## example of script execution:
## bash snp_explorer_cp.sh | tee stdout.txt


# use first cmd argument as filepath containing list of \n-delim SNP IDs

printf '\nEnter filepath for tsv file containing disease names (first column) and SNP IDs (second column), and press [ENTER]:\n'

read filepath

printf '\nEnter window size (in bp) for searches and press [ENTER]:\n'

read window

printf '\nEnter threshold adjusted p value for significance calls (numbers only, ie. 0.005) and press [ENTER]:\n'

read threshold


# print header to tsv file for QTL results

resultsNameQTL='SE_QTL_results.tsv'

printf 'Disease Name\tSNP ID\tQTL Type\tFeature Name\tp value unadj\tNearby Features\tp value adj\tDistance to SNP\n' > $resultsNameQTL


# print header to tsv file for bayes net results

resultsNameBN='SE_BN_results.tsv'

printf 'Disease Name\tSNP Prevalence\tSNP ID\tEnsembl ID\tGene Name\tGene Distance to SNP\tGeno-Exp pvalue\tNearby Genes\tGeno-Exp pvalue adj\tProbe ID\tProbe Distance to SNP\tGeno-Met pvalue\tNearby Probes\tGeno-Met pvalue adj\tPeak Name\tPeak Distance to SNP\tGeno-Acet pvalue\tNearby Peaks\tGeno-Acet pvalue adj\tModel Score\tModel Structure\tMet-Exp pvalue\tAcet-Exp pvalue\tMet-Acet pvalue\n' > $resultsNameBN


# calculate size (number of lines) of 10 subset files

subset_size=($(wc -l $filepath))
subset_size=$((subset_size / 10))

# run on single file if less than 100 lines
# otherwise, split file into 10 subsets and
# run pipeline on each, in parallel

if [ "$subset_size" -lt "10" ]; then
        bash snp_explorer_m0.sh $filepath $window $window $threshold
else
        split -l $subset_size $filepath subset_SNPs_
	for this_subset in ./subset_SNPs_*
	do
		bash snp_explorer_m0.sh "$this_subset" $window $window $threshold &
	done 
fi

# wait for all threads to complete before proceeding

wait

# remove temporary subset files

rm subset_SNPs_a*
