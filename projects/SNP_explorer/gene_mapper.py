#!/usr/bin/env python

#the following script takes a list of gene names and finds the nearby SNP that is most correlated with expression of that gene

import sys
import glob
import csv
import numpy

def myFinder(anItem,aList):
	try:
		index = aList.index(anItem)
	except ValueError:
		return -1
	return index

# file paths for my Broad account directory
filename1 = "/home/unix/jtopham/data/snp_mapper/snpArrayAffy6.txt"            # maps SNP ids to genomic position
filename2 = "/home/unix/jtopham/data/ensembl_mapper2.txt" 		      # maps gene names to  genomic position
filename3 = "/home/unix/jtopham/data/expression/residual_gene_expression_expressed_genes_2FPKM100ind.txt"  # expression data
filename4 = "/home/unix/jtopham/data/phenotype/finalptlist.txt"               # list of patients that have all 4 data types available
filename5 = "/home/unix/jtopham/data/methylation/filtered_meth_mapper.txt"    # maps CpG sites to genomic position
filename6 = "/home/unix/jtopham/data/acetylation/peakInfo.csv"                # maps ChipSeq peaks to genomic position 
filename7 = "/home/unix/jtopham/data/methylation/myProbeList.txt"             # list of probes- taken from methylation data matrix
filename8 = "/home/unix/jtopham/data/methylation/ill450kMeth_all_740_imputed.txt" # FULL methylation data matrix
filename9 = "/home/unix/jtopham/data/phenotype/p3.txt"                        # patient mapping file
filename10 = "/home/unix/jtopham/data/acetylation/peakList.txt"  	      # list of peaks - taken from acetylation data matrix
filename11 = "/home/unix/jtopham/data/acetylation/chipSeqResiduals.csv"	      # FULL acetylation data matrix

args = sys.argv     					     #  arguments should be: gene name, SNP window, probe window

# boiler plate:
exp_vector = []
exp_colname = -1
exp_pts = []
nearby_snps = []
snp_chrs = []
final_pt_list = []
new_exp_vector = [[], []]
meth_probes = []
acetyl_peaks = []
meth_probes2 = []
good_pts1 = []
meth_pts1 = []
meth_pts2 = []
meth_intersect = []
final_meth_vector = []
peak_list = []

this_gene = args[1]

with open(filename4,"r") as ins:
	for line in ins:
		final_pt_list.append(line.strip('\n'))

with open(filename10, "r") as ins:
	for line in ins:
		peak_list.append(line)

with open(filename3,"r") as ins:             # find which column to extract for expression vector, and construct list of patients from expression data
	j = 0
	for line in ins:
		tmp_pt = line.split()[0]
		if j==0:                             # only search header
			this_split = line.split()
			i = -1
			for item in this_split:
				i += 1
				if item.split(sep=":")[0] == this_gene:
					exp_colname = i
		else:
			exp_pts.append(tmp_pt.split(sep=":")[0])
		j += 1

if exp_colname == -1:
	sys.exit('Error: gene not found in data matrix')


with open(filename3, "r") as ins:           # construct the expression vector
	i = 0
	for line in ins:
		if i == 0:                          # skip header
			i += 1
			continue
		exp_vector.append(line.split()[(exp_colname + 1)])		      # correct for column-shift issue

exp_intersect = [val for val in exp_pts if val in final_pt_list]		  # find patients of interest

new_exp_vector[1].append(this_gene)

for patient in exp_intersect:                                             # filter expression vector
	new_exp_vector[1].append(exp_vector[exp_pts.index(patient)])

exp_intersect.insert(0, " ")

for pt in exp_intersect:
	new_exp_vector[0].append(pt)

new_exp_vector2 = numpy.transpose(new_exp_vector)

gene_position = -1
with open(filename2, "r") as ins:    			                      # get start position of gene
	for line in ins:
		if line.split()[1] == this_gene:
			gene_position = int(line.split()[3])
			break

if gene_position == -1:
	sys.exit('Error: gene maps to no positon')


myProbeList = []
with open(filename7, "r") as ins:
	for line in ins:
		myProbeList.append(line.strip('\n'))                         # read in list of probe IDs

pt_map1 = []
pt_map2 = []
with open(filename9, "r") as ins:                                            # read in patient mapper
	for line in ins:
		pt_map1.append(line.split()[0])
		pt_map2.append(line.split()[1])

g_window = int(args[2])
m_window = int(args[3])                                                      # set methylation window here
a_window = 50000                                                             # set acetylation window here


m_min_range = gene_position - m_window/2
m_max_range = gene_position + m_window/2

a_min_range = gene_position - int(a_window)/2
a_max_range = gene_position + int(a_window)/2

g_min_range = gene_position - g_window/2
g_max_range = gene_position + g_window/2

def minChecker(min):
	if min < 0:
		return 0
	else:
		return(int(min))

g_min_range = minChecker(g_min_range)
m_min_range = minChecker(m_min_range)
a_min_range = minChecker(a_min_range)

# find SNPs within window of gene start position
with open (filename1, "r") as ins:                                  
	for line in ins:
		if int(line.split()[2]) >= g_min_range and int(line.split()[2]) <= g_max_range:
			nearby_snps.append(line.split()[8])
			snp_chrs.append(line.split()[1])

if len(nearby_snps) == 0:
	sys.exit('Error: no SNPs found within supplied window')

# find methylation probes within window of gene start position
with open(filename5, "r") as ins:
	for line in ins:
		if line.split()[1].isdigit():
			if int(line.split()[1]) >= m_min_range and int(line.split()[1]) <= m_max_range:
				meth_probes.append(line.split()[0])

if len(meth_probes) == 0:
	sys.exit('Error: no methylation probes found within supplied window')

# find acetyltion peaks within some window of gene start position
with open(filename6, "r") as ins:
	for line in ins:
		if line.split()[1].isdigit():
			this_mean = (int(line.split()[1]) + int(line.split()[2])) / 2
			if this_mean >= a_min_range and this_mean <= a_max_range:                       # find acetylation probes within chosen window of SNP
				acetyl_peaks.append(line.split()[0])

print('probes found: ', len(meth_probes))
print('peaks found: ', len(acetyl_peaks))
print('SNPs found: ', len(nearby_snps))

meth_matrix = [[] for _ in range(len(meth_probes) + 1)]

z = 0
for probe in meth_probes:
	meth_vector2 = []
	try:
		mIndex = myProbeList.index(probe)
	except:
		continue                                             # skip probe if it is not in data matrix

	meth_probes2.append(probe)                                   # keep track of probes that will have correlation calculated       
	i = -1
	with open (filename8, "r") as ins:
		for line in ins:
			i += 1
			if i == 0 and z == 0:                        # only grab list of patients when iterating through data matrix the first time
				good_pts2 = []
				for pt in line.split():
					meth_pts1.append(pt)
					if pt in pt_map1:
						good_pts1.append(pt)                            # keep list of pts that can be mapped
						good_pts2.append(pt_map2[pt_map1.index(pt)])    # mapping is done here
				for pt in good_pts2:
					tmp = pt.replace('ROS','')                              # process patient IDs (removing 'ROS' and 'MAP')
					meth_pts2.append(tmp.replace('MAP', ''))
				meth_intersect = [val for val in meth_pts2 if val in final_pt_list]
			if i == mIndex:
				meth_vector1 = line.split()          # extract  methylation data vector for this probe
				break

	for pt in meth_intersect:
		meth_vector2.append(float(meth_vector1[meth_pts1.index(good_pts1[meth_pts2.index(pt)])]))       # filter the methylation data vector for intersect pts
	meth_vector2.insert(0,probe)
	for val in meth_vector2:
		meth_matrix[z].append(val)
	z += 1

meth_intersect.insert(0, " ")

for pt in meth_intersect:
	meth_matrix[z].append(pt)

meth_matrix[z], meth_matrix[0] = meth_matrix[0], meth_matrix[z]    # place rownames into first column of data matrix

for rm in range(z+1, len(meth_matrix)):		# remove blank columns from matrix
	meth_matrix.pop()

print('dim of meth matrix: ', len(meth_matrix),  len(meth_matrix[0]))

meth_matrix2 = numpy.transpose(meth_matrix)

nearby_snps_corrs = []
snp_locations = []
snp_colnames = []
snp_matrix = [[] for _ in range(len(nearby_snps) + 1)]
geno_pts = []                            
                            
z = 0
for i in range(0,len(nearby_snps)):                  # find column index and file of ith SNP
	genoFiles = glob.glob("/home/unix/jtopham/data/genotype/" + snp_chrs[i] + "*.txt")
	colname_coord = -1
	for file in genoFiles:
		with open(file,"r") as ins:
			for line in ins:
				colname_coord = myFinder(nearby_snps[i], line.split())
				break                                                # take only header
		if colname_coord >= 0:						     # stop searching files if SNP has been located
			true_habitat = file
			break

	if colname_coord == -1:							     # skip this SNP if it was not included on SNP Chip
		continue

	snp_locations.append(true_habitat)					     # keep track of SNP file location and column
	snp_colnames.append(colname_coord)

# get vectors of genotype data and patient IDs

	geno_vector = []

	with open(true_habitat, "r") as ins:
		p = 0
		for line in ins:
			if p==0:
				p += 1
				continue					     # skip first line
			if z == 0:
				tmp_pt = line.split()[0]
				geno_pts.append(''.join(list(tmp_pt)[3:]))	    # store vector of patient ids (removing 'ROS' and 'MAP')

			geno_vector.append(line.split()[(colname_coord + 1)])	     # correct for colname-shifting issue when indexing rows (see May 26th entry in log for details)
                            
	geno_intersect = [val for val in geno_pts if val in final_pt_list]		 # find patients of interest

	new_geno_vector = []

	for patient in geno_intersect:						     # filter genotype vector for patients of interest
		new_geno_vector.append(float(geno_vector[geno_pts.index(patient)]))

	new_geno_vector.insert(0,nearby_snps[i])
	for val in new_geno_vector:
		snp_matrix[z].append(val)    #try appending vals individually rather than whole vector 
	z += 1

for rm in range(z+1, len(snp_matrix)):		# remove blank columns from matrix
	snp_matrix.pop()

geno_intersect.insert(0," ")

for pt in geno_intersect:                       # add rowname vector to matrix
	snp_matrix[z].append(pt)

snp_matrix[z], snp_matrix[0] = snp_matrix[0], snp_matrix[z]      # place rownames (pts) into first column

print('dim of geno matrix:', len(snp_matrix), len(snp_matrix[0]))

snp_matrix2 = numpy.transpose(snp_matrix)

# output acetylation vector/matrix, if any peaks found

acetyl_vector1 = []
acetyl_vector2 = []
acetyl_pts = []
acetyl_matrix = [[] for _ in range(len(acetyl_peaks)+1)]
acetyl_intersect = []

z = 0
if len(acetyl_peaks) != 0:
	for peak in acetyl_peaks:
		try:
			pIndex = peak_list.index(peak)
		except:
			continue
		i = -1
		with open(filename11, "r") as ins:
			for line in ins:
				i += 1
				if i == 0 and z == 0:
					acetyl_pts = line.split()
				if i == pIndex:
					acetyl_vector1 = line.split()
					break

		# need to filter vector here
		acetyl_intersect = [val for val in acetyl_pts if val in final_pt_list]
		for pt in acetyl_intersect:
			acetyl_vector2.append(acetyl_vector1[acetyl_pts.index(pt)])

		for val in acetyl_vector2:
			acetyl_matrix[z].append(val)
		z += 1

	if z >= 1:						   # if at least one peak mapped
		acetyl_intersect.insert(0, " ")
		for pt in acetyl_intersect:
			acetyl_matrix[z].append(pt)
		acetyl_matrix[0], acetyl_matrix[z] = acetyl_matrix[z], acetyl_matrix[0]

		print('acetyl intersect \t aceyl vector2 ')
		for i in range(0,10):
			print(acetyl_intersect[i], '\t', acetyl_vector2[i])

		acetyl_matrix2 = numpy.transpose(acetyl_matrix)
		aFile = this_gene + "acetylation.txt"
		with open(aFile, "wt") as f:
			writer = csv.writer(f, delimiter = '\t')
			writer.writerows(acetyl_matrix2)

gFile = this_gene + "_genotype.txt" 				     # write genotype vector to file
with open(gFile, "wt") as f:
	writer = csv.writer(f,delimiter = '\t')
	writer.writerows(snp_matrix2)

eFile = this_gene + "_expression.txt"				     # write expression vector to file
with open(eFile, "wt") as f:
	writer = csv.writer(f, delimiter = '\t')
	writer.writerows(new_exp_vector2)

mFile = this_gene + "_methylation.txt"				     # write methylation matrix to file
with open(mFile, "wt") as f:
	writer = csv.writer(f, delimiter = '\t')
	writer.writerows(meth_matrix2)

