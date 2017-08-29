#!/usr/bin/env python


# this script finds all nearby SNPs, methylation probes and acetylation peaks within some window of a gene start site
# the input is a gene name and windows for SNps / peaks. ie.    python genome_algebra.py LSM1 200000 200000
# this script outputs data matrices for nearby probes and peaks, as this can be done faster in python
# it ouputs names of nearby SNPs, for which data matrix will be constructed using R, which is faster than python
# the rationale being that specific rows and columns can be read from large data matrices faster in python and R (respectively)

import sys
import numpy
import csv

# define function for checking min value of windows (relative to gene position)

def minChecker(min):
	if min < 0:
		return 0
	else:
		return(int(min))


# file paths for my Broad account directory

filename1 = "/home/unix/jtopham/data/ensembl_mapper2.txt" 		      # maps gene names to  genomic position
filename2 = "/home/unix/jtopham/data/snp_mapper/snpArrayAffy6.txt"            # maps SNP ids to genomic position
filename3 = "/home/unix/jtopham/data/methylation/filtered_meth_mapper.txt"    # maps CpG sites to genomic position
filename4 = "/home/unix/jtopham/data/acetylation/peakInfo.csv"                # maps ChipSeq peaks to genomic position 
filename5 = "/home/unix/jtopham/data/methylation/probe_list.txt"              # list of probes, in order, from ill450k array
filename6 = "/home/unix/jtopham/data/methylation/ill450kMeth_all_740_imputed.txt"  # methylation data matrix
filename7 = "/home/unix/jtopham/data/acetylation/peakList.txt"		      # list of peaks, in order, from chipSeq assay
filename8 = "/home/unix/jtopham/data/acetylation/chipSeqResiduals.csv"        # acetylation data matrix

# boiler plate:

args = sys.argv     		 #  arguments should be: gene name, SNP window, probe window
this_gene = args[1]
nearby_snps = []
nearby_snps.append("SNP\tLocation")
nearby_probes = []
nearby_peaks = []


# get start position of gene

gene_position = -1
with open(filename1, "r") as ins:    			                      
	for line in ins:
		if line.split()[1] == this_gene:
			gene_position = int(line.split()[3])
			gene_chr = "chr" + line.split()[2]
			break

if gene_position == -1:
	sys.exit('Error: gene maps to no positon')


# arrange windows for feature search space

g_window = int(args[2])
m_window = int(args[3])                                                      # set methylation window here
a_window = 1000000                                                           # set acetylation window here (1mb)

m_min_range = gene_position - m_window/2
m_max_range = gene_position + m_window/2

a_min_range = gene_position - int(a_window)/2
a_max_range = gene_position + int(a_window)/2

g_min_range = gene_position - g_window/2
g_max_range = gene_position + g_window/2

g_min_range = minChecker(g_min_range)
m_min_range = minChecker(m_min_range)
a_min_range = minChecker(a_min_range)


# find SNPs within window of gene start position

with open (filename2, "r") as ins:                                  
	for line in ins:
		if line.split()[1] == gene_chr:
			if int(line.split()[2]) >= g_min_range and int(line.split()[2]) <= g_max_range:
				nearby_snps.append(line.split()[8] + '\t' + line.split()[1])

if len(nearby_snps) == 0:
	sys.exit('Error: no SNPs found within supplied window')



# find methylation probes within window of gene start position

with open(filename3, "r") as ins:
	for line in ins:
		if line.split()[1].isdigit():
			if int(line.split()[1]) >= m_min_range and int(line.split()[1]) <= m_max_range:
				nearby_probes.append(line.split()[0])

if len(nearby_probes) == 0:
	sys.exit('Error: no methylation probes found within supplied window')



# find acetyltion peaks within some window of gene start position

with open(filename4, "r") as ins:
	for line in ins:
		if line.split()[1].isdigit():
			this_mean = (int(line.split()[1]) + int(line.split()[2])) / 2
			if this_mean >= a_min_range and this_mean <= a_max_range:                       # find acetylation peaks within chosen window of SNP
				nearby_peaks.append(line.split()[0])

if len(nearby_peaks) == 0:
	sys.exit('Error: no acetylation peaks found within 1mb window')


# std.out summary of what was found


print('Window search summary for ', this_gene, ': ')
print('Probes found:\t', len(nearby_probes))
print('Peaks found:\t', len(nearby_peaks))
print('SNPs found:\t', len(nearby_snps))



# retrieve methylation matrix faster than R

## determine indeces of nearby probes within data matrix
## read-in only these indeces, write to tsv file


mIndeces = []
fp = open(filename5)
for i, line in enumerate(fp):
	if line.split()[0] in nearby_probes:
		mIndeces.append(i)
fp.close()

if len(mIndeces)<1:
	sys.exit('No nearby probes were found in methylation array')

mIndeces.sort()

meth_matrix = []
fp = open(filename6)
for i, line in enumerate(fp):
	if i == 0:
		meth_matrix.append(line.split())
	if i in mIndeces:
		meth_matrix.append(line.split())
	if i > mIndeces[-1]:
		break
fp.close()

mFile = "z_" + this_gene + "_methylation.txt"                                # write methylation matrix to file
with open(mFile, "wt") as f:
	writer = csv.writer(f, delimiter = '\t')
	writer.writerows(meth_matrix)



# retrieve acetylation matrix faster than R

## retrieve indeces of nearby peaks in data matrix
## read-in only these indeces and write to tsv file


pIndeces = []
fp = open (filename7)
for i, line in enumerate(fp):
	if line.split()[0] in nearby_peaks:
		pIndeces.append(i)
fp.close()

if len(pIndeces)<1:
	sys.exit('No nearby peaks were found in acetylation array')

pIndeces.sort()

acet_matrix = []
fp = open(filename8)
for i, line in enumerate(fp):
	if i == 0:
		acet_matrix.append(line.split())
	if i in pIndeces:
		acet_matrix.append(line.split())
	if i > pIndeces[-1]:
		break
fp.close()

aFile = "z_" + this_gene + "_acetylation.txt"
with open(aFile,"wt") as f:
	writer = csv.writer(f,delimiter = '\t')
	writer.writerows(acet_matrix)


# output names of nearby SNPs to file
# data matrix for these will be retrieved by R (faster)


gFile = "z_" + this_gene + "_genotype.txt" 				  				     

with open(gFile, "w") as f:
	f.write("\n".join(nearby_snps))
	f.write("\n")
	f.close()

