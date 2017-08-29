#!/usr/bin/env python


# this script finds all nearby gene TSSs, methylation probes and acetylation peaks within some window of a SNP
# the input is a SNP ID and windows for genes / probes. ie.  python snp_explorer.py rs9629043 200000 200000
# it ouputs names of nearby genes, probes and peaks, for which data matrix will be constructed using R

import sys
import numpy
import csv

# define function for checking min value of windows (relative to SNP position)

def minChecker(min):
	if min < 0:
		return 0
	else:
		return(int(min))


# file paths for local machine


filename1 = "/Users/james/Desktop/jtopham_prj/data/rosmap/snpArrayAffy6.txt"           # maps SNP ids to genomic position
filename2 = "/Users/james/Desktop/jtopham_prj/data/rosmap/ensembl_mapper3.txt" 	       # maps gene names to  genomic position
filename3 = "/Users/james/Desktop/jtopham_prj/data/rosmap/filtered_CpG_map.tsv"        # maps CpG sites to genomic position
filename4 = "/Users/james/Desktop/jtopham_prj/data/rosmap/peakInfo.tsv"                # maps ChipSeq peaks to genomic position

# boiler plate:

args = sys.argv     		 #  arguments should be: gene name, SNP window, probe window
this_snp = args[1]
nearby_genes = []
nearby_probes = []
nearby_peaks = []


# get start position of SNP

snp_position = -1
with open(filename1, "r") as ins:    			                      
	for line in ins:
		if (line.split()[4] == this_snp or line.split()[8] == this_snp):
			snp_position = int(line.split()[2])
			snp_chr = line.split()[1]
			break


if snp_position == -1:
	sys.exit('Error ' + '(' + this_snp + '): unable to map SNP to position')


snp_chr = snp_chr[3:]        # use only chromosome number, without preceding 'chr' string


# arrange windows for feature search space

g_window = int(args[2])                                                      # set gene window here
m_window = int(args[3])                                                      # set methylation window here
a_window = 1000000                                                           # set acetylation window here (1mb default)

m_min_range = snp_position - m_window/2
m_max_range = snp_position + m_window/2

a_min_range = snp_position - int(a_window)/2
a_max_range = snp_position + int(a_window)/2

g_min_range = snp_position - g_window/2
g_max_range = snp_position + g_window/2

g_min_range = minChecker(g_min_range)                                        # ensure no negative window values
m_min_range = minChecker(m_min_range)
a_min_range = minChecker(a_min_range)


# find genes within window of snp start position

with open (filename2, "r") as ins:                                  
	for line in ins:
		if line.split()[2] == snp_chr:
			if int(line.split()[3]) >= g_min_range and int(line.split()[3]) <= g_max_range:
				this_gene = "ENSG" + line.split()[0] + '\t' + line.split()[1]
				nearby_genes.append(this_gene)

if len(nearby_genes) == 0:
	sys.exit('Error ' + '(' + this_snp + '): no genes found within designated window')



# find methylation probes within window of snp start position

with open(filename3, "r") as ins:
	for line in ins:
		if line.split()[1] == snp_chr:
			if int(line.split()[2]) >= m_min_range and int(line.split()[2]) <= m_max_range:
				nearby_probes.append(line.split()[0])

if len(nearby_probes) == 0:
	sys.exit('Error ' + '(' + this_snp + '): no methylation probes found within designated window')



# find acetyltion peaks within some window of gene start position

with open(filename4, "r") as ins:
	for line in ins:
		if line.split()[3] == snp_chr:
			this_mean = (int(line.split()[1]) + int(line.split()[2])) / 2
			if this_mean >= a_min_range and this_mean <= a_max_range:                       # find acetylation peaks within chosen window of SNP
				nearby_peaks.append(line.split()[0])

if len(nearby_peaks) == 0:
	sys.exit('Error ' + '(' + this_snp + '): no acetylation peaks found within 1mb window')


# std.out summary of what was found

print(' ')
print('Window search summary for ' + this_snp + ': ')
print('Nearby probes:     ' +  str(len(nearby_probes)))
print('Nearby peaks:      ' + str(len(nearby_peaks)))
print('Nearby gene TSS:   ' + str(len(nearby_genes)))
print(' ')

# output names of nearby features to individual files
# data matrix for these will be retrieved by R in next
# module (myBayesNet.R)


gFile = "z_" + this_snp + "_genes.txt"

with open(gFile, "w") as f:
	f.write("\n".join(nearby_genes))
	f.write("\n")
	f.close()

mFile = "z_" + this_snp + "_probes.txt"                                # write methylation matrix to file

with open(mFile, "w") as f:
	f.write("\n".join(nearby_probes))
	f.write("\n")
	f.close()

aFile = "z_" + this_snp + "_peaks.txt"

with open(aFile, "w") as f:
	f.write("\n".join(nearby_peaks))
	f.write("\n")
	f.close()
