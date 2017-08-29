#!/usr/bin/env python

#the following script takes a list of gene names and finds the nearby SNP that is most correlated with expression of that gene

import sys
import scipy
import glob
from scipy.stats.stats import pearsonr

filename = "/home/unix/jtopham/data/snp_mapper/snpArrayAffy6.txt"            # maps SNP ids to genomic position
filename2 = "/home/unix/jtopham/data/ensembl_mapper2.txt" 		     # maps gene names to  genomic position
filename3 = "/home/unix/jtopham/data/expression/residual_gene_expression_expressed_genes_2FPKM100ind.txt"  # expression data

args = sys.stdin.readlines()        					     # list of \n-delim gene names supplied by stdin ( ie. more genelist.txt | python thisscript.pl ))

for arg in args:
	this_gene = arg.strip('\n')
	exp_vector = []
	exp_colname = -1
	exp_pts = []

	with open(filename3,"r") as ins:
		j = 0
		for line in ins:
			this_split = line.split()
			tmp_pt = line.split()[0]
			if j == 0:							     # check only header for gene name
				i = -1
				for item in this_split:
					i += 1
					if item.split(sep=":")[0] == this_gene:
						exp_colname = i				     # determine which column to extract expression values from
			else:
				exp_pts.append(tmp_pt.split(sep=":")[0])		     # construct vector of patient IDs
			j += 1

	if exp_colname == -1:
		#print('Warning: ', this_gene, ' not found in expression data matrix')
		continue

	with open(filename3, "r") as ins:
		i = 0
		for line in ins:
			if i == 0:
				i += 1
				continue
			exp_vector.append(line.split()[(exp_colname + 1)])		      # construct vector of expression values for this gene correcting for column-shift issue

	gene_position = []
	nearby_snps = []

	with open(filename2, "r") as ins:    			             # get start position of gene
		for line in ins:
			if line.split()[1] == this_gene:
				gene_position.append(int(line.split()[3]))

	#if len(gene_position)>1:
		#print('Warning:', this_gene, 'maps to more than one position!')

	if len(gene_position)==0:
		#print(this_gene, " maps to no positon")
		continue

	window = 8000			                                     # set window here
	min_range = int(gene_position[0]) - int(window)/2
	max_range = int(gene_position[0]) + int(window)/2

	with open (filename, "r") as ins:				     # find SNPs that are near the start of this gene, within some window
		for line in ins:
			if int(line.split()[2]) >= min_range and int(line.split()[2]) <= max_range:
				nearby_snps.append(line.split()[8])

	if len(nearby_snps)==0:
		#print('Warning: no SNPs found within ', window, ' bp of ', this_gene)
		continue

	# now determine SNP most highly correlated with expression of this gene

	genoFiles = glob.glob("/home/unix/jtopham/data/genotype/*.txt")

	nearby_snps_corrs = []
	nearby_snps2 = []								     # this is used to keep track of SNPs that are found on the SNP chip- for indexing purposes
	snp_locations = []
	snp_colnames = []

	for this_snp in nearby_snps:
		colname_coord = -1
		for file in genoFiles:                                                       # find snp id amongst 744 file headers
			with open(file,"r") as ins:
				for line in ins:
					i = -1
					for name in line.split():
						i += 1
						if name == this_snp:
							colname_coord = i                    # record which column (if any, in this file) corresponds to SNP of interest
					break                                                # take only header
			if colname_coord >= 0:						     # stop searching files if SNP has been located
				true_habitat = file
				break

		if colname_coord == -1:							     # skip this SNP if it was not included on SNP Chip
			continue

		nearby_snps2.append(this_snp)						     # keep track of SNPs that will have correlation coefficient calculated
		snp_locations.append(true_habitat)					     # keep track of SNP file location and column
		snp_colnames.append(colname_coord)

		# get vectors of genotype data and patient IDs

		geno_vector = []
		geno_pts = []
		with open(true_habitat, "r") as ins:
			i = 0
			for line in ins:
				if i==0:
					i += 1
					continue					     # skip first line
				tmp_pt = line.split()[0]			             
				geno_pts.append(''.join(list(tmp_pt)[3:]))	             # store vector of patient ids (removing 'ROS' and 'MAP')
				geno_vector.append(line.split()[(colname_coord + 1)])	     # correct for colname-shifting issue when indexing rows (see May 26th entry in log for details)

		
		intersect = [val for val in exp_pts if val in geno_pts]			     # get list of patients with both expression and genotype data 

		new_geno_vector = []
		new_exp_vector = []

		for patient in intersect:						     # get list of indeces for the two data vectors' intersection patients
			new_geno_vector.append(float(geno_vector[geno_pts.index(patient)]))
			new_exp_vector.append(float(exp_vector[exp_pts.index(patient)]))

		# calculate spearman correlation using scipy
			
		correlation = pearsonr(new_geno_vector, new_exp_vector)[0]
		#print(this_gene, '\t and \t', this_snp, '\t',"Pearson Correlation:",'\t', correlation)
		
		nearby_snps_corrs.append(float(abs(correlation)))

	bestsnp = nearby_snps_corrs.index(max(nearby_snps_corrs))			     # this is the index of the best SNP within the list of SNPs that had correlations calculated
	final_geno_vector = []

	#print('Best Candidate for', this_gene, ' is: ', nearby_snps2[bestsnp], 'with absolute score: ', max(nearby_snps_corrs))

	with open(snp_locations[bestsnp], "r") as ins:					     # store genotype vector for best SNP
		i = 0
		for line in ins:
			if i == 0:
				i += 1
				continue
			final_geno_vector.append(line.split()[(snp_colnames[bestsnp] + 1)])

	gFile = nearby_snps2[bestsnp] + "_genotypevector.txt" 				     # write genotype vector to file
	genoFile = open(gFile,'w+')
	genoFile.write("\n".join(final_geno_vector))
	genoFile.close()

	eFile = nearby_snps2[bestsnp] + "_expressionvector.txt"				     # write expression vector to file
	expFile = open(eFile,'w+')
	expFile.write("\n".join(exp_vector))
	expFile.close()

	pFile = nearby_snps2[bestsnp] + "_patientvector.txt"				     # write list of patients corresponding to genotype vector to file
	ptFile = open(pFile,'w+')
	ptFile.write("\n".join(intersect))
	ptFile.close()


	print(nearby_snps2[bestsnp])							     # output SNP id to STDOUT, for piping purposes
