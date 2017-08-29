#!/bin/bash

## this script is the zeroth module of the SNP Explorer pipeline
## it serves in facilitating a threading process, in which the
## list of SNPs is split into 10 subsets which are analyzed in
## parallel

## it takes cmd line args as input:
##	1. filepath of subset of input tsv file
##	2. window for probes
##	3. window for peaks
##	4. threshold pvalue for significance calls

## example execution:
## bash snp_explorer_m0.sh inputSNPs.txt 1000000 1000000 0.05


resultsNameQTL='SE_QTL_results.tsv'

resultsNameBN='SE_BN_results.tsv'


# send SNPs through pipe, utilizing user-input and disease/SNP names from input tsv file

cat $1 | while read disease snp
do

if python snp_explorer_m1.py $snp $2 $3; then                           	    # enter search windows here
	Rscript snp_explorer_m2.R $snp $resultsNameBN $resultsNameQTL $4 $disease
	rm '/Users/james/Desktop/jtopham_prj/scripts/rosmap/z_'$snp'_'*		    # remove tmp files from py script to clean up workspace
fi

done

