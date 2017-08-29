#!/usr/bin/env python

# the purpose of this script is to parse a list of SNPs from the 753 imputed genotype files
# from the Broad server.

# the script takes as input a list of \n-delim SNP IDs (ie. rs1234..) and returns a
# file (parsed_imputed_SNPs.tsv) containing the sample IDs (first col) along with
# vector of genotype data for each SNP. note that the input exists as a filename
# specified by user as first argument

# example execution:
# python snp_grabber.py /home/unix/jtopham/data/my_snp_list.txt



import glob
import numpy
import csv
import sys


args = sys.argv                  #  argument should be path to file containing list of \n-delim SNP ids

usr_file = args[1]
snp_list = []

with open(usr_file, "r") as ins:
	for line in ins:
		snp_list.append(line.strip('\n'))


snp_matrix = [[] for _ in range(len(snp_list) + 1)]

genoFiles = glob.glob("/home/unix/jtopham/data/genotype/*")



# place sampleIDs into first column of output matrix

fp = open (genoFiles[0])
for i, line in enumerate(fp):
	if i==0:
		snp_matrix[0].append('sampleID')
	else:
		snp_matrix[0].append(line.split()[0])
fp.close()


# grab data vectors for each SNP out of the 753 imputed genotype files

z = 0
for file in genoFiles:
	this_header = []
	these_indeces = []
	fp = open (file)
	for i, line in enumerate(fp):
		if i==0:
			for snp in line.split():
				this_header.append(snp)			 # parse header only for this file
		else:
			break
	fp.close()

	for snp in snp_list:
		try:
			snpIndex = this_header.index(snp)		 # attempt to index header for snps of interest
			these_indeces.append(snpIndex+1)		 # account for missing 'sampleID' header
		except:
			continue


	if len(these_indeces) > 0:
		for index in these_indeces:
			z+=1
			this_vector = []
			fp = open (file)
			for i, line in enumerate(fp):
				if i==0:
					this_vector.append(line.split()[index-1])	# get correct column name
				else:
					this_vector.append(line.split()[index])		# get vector of data
			fp.close()

			for val in this_vector:
				snp_matrix[z].append(val)


# error handling

if z!=len(snp_list):
	print('Warning: only ', z, ' of ', len(snp_list), ' found.')
	for rm in range(z+1, len(snp_matrix)):						# remove blank columns from matrix
		snp_matrix.pop()


# transpose matrix

snp_matrix2 = numpy.transpose(snp_matrix)


# write output matrix to tsv file

gFile = "parsed_imputed_SNPs.tsv"
with open(gFile, "wt") as f:
	writer = csv.writer(f,delimiter = '\t')
	writer.writerows(snp_matrix2)





