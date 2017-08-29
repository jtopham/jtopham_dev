#!/usr/bin/env python

#the following script takes a list of SNP ids and maps them to nearby methlyation probes and acetylation peaks

import sys
import scipy
from scipy.stats.stats import pearsonr

filename = "/home/unix/jtopham/data/snp_mapper/snpArrayAffy6.txt"            # maps SNP ids to genomic position
filename2 = "/home/unix/jtopham/data/methylation/filtered_meth_mapper.txt"   # maps CpG sites to genomic position
filename3 = "/home/unix/jtopham/data/acetylation/peakInfo.csv"               # maps ChipSeq peaks to genomic position 
filename5 = "/home/unix/jtopham/data/methylation/myProbeList.txt"	     # list of probes- taken from methylation data matrix
filename6 = "/home/unix/jtopham/data/methylation/ill450kMeth_all_740_imputed.txt" # FULL methylation data matrix
filename7 = "/home/unix/jtopham/data/phenotype/p3.txt"			     # patient mapping file

snp_mapper = []
with open(filename, "r") as ins:     					     # read in file contining list of SNP ids and their corresponding genomic positions
	for line in ins:
        	snp_mapper.append(line.strip('\n'))

myProbeList = []
with open(filename5, "r") as ins:
	for line in ins:
		myProbeList.append(line.strip('\n'))                 	     # read in list of probe IDs

pt_map1 = []
pt_map2 = []
with open(filename7, "r") as ins:                                            # read in patient mapper
	for line in ins:
		pt_map1.append(line.split()[0])
		pt_map2.append(line.split()[1])

args = sys.stdin.readlines()        					     # this is a list of \n-delim SNPs supplied by stdin ( ie. more snplist.txt | python thisscript.pl )

for arg in args:
	this_snp = arg.strip('\n')

	i = 0
	for line in snp_mapper:
		if line.split()[8] == this_snp:
			this_position = line.split()[3]
			i += 1
			break
	if i == 0:
		print(this_snp, 'dont map')
		continue # snp not mappable

	with open(filename3, "r") as ins:
		for window in range(500,2000000,2000):
			acetyl_peaks = []
			a_min_range = int(this_position) - int(window)/2
			a_max_range = int(this_position) + int(window)/2
			if a_min_range < 0:
				a_min_range = 0
			for line in ins:
				if line.split()[1].isdigit():
					this_mean = (int(line.split()[1]) + int(line.split()[2])) / 2
					if this_mean >= a_min_range and this_mean <= a_max_range:                       # find acetylation probes within chosen window of SNP
						acetyl_peaks.append(line.split()[0])
			print(window, '\t', len(acetyl_peaks))

	continue#only test acetylation window


	geno_vector = []
	geno_vector2 = []
	filename4 = "/home/unix/jtopham/scripts/" + this_snp + "_genotypevector.txt"   		# name of genotype vector for this SNP,  output by previous pipe module

	with open(filename4, "r") as ins:
		for line in ins:
			geno_vector.append(line.strip('\n'))

	geno_pts = []
	filename8 = "/home/unix/jtopham/scripts/" + this_snp + "_patientvector.txt"		# vector of patient IDs for genotype vector, output by previous pipe module 

	with open(filename8, "r") as ins:
		for line in ins:
			geno_pts.append(line.strip('\n'))

	this_position = []

	for line in snp_mapper:     					     # identify position of this SNP
		if line.split()[8] == this_snp:
			this_position.append(line.split()[3])

	#if len(this_position)>1:
		#print('Warning:', this_snp, 'maps to more than one position!')

	if len(this_position)==0:
   	        #print('Warning:', this_snp, 'does not map to any position')
		unmapped_snps_counter += 1
		continue

	m_window = 3000                                			     # set methylation window here
	a_window = 100000			      			     # set acetylation window here

	m_min_range = int(this_position[0]) - int(m_window)/2
	m_max_range = int(this_position[0]) + int(m_window)/2

	a_min_range = int(this_position[0]) - int(a_window)/2
	a_max_range = int(this_position[0]) - int(a_window)/2

	if m_min_range < 0:          					     # ensure min_range doesn't fall below zero (likely unnecessary)
		m_min_range = 0

	if a_min_range < 0:
		a_min_range = 0

	meth_probes = []
	acetyl_peaks = []

	with open(filename2, "r") as ins:     				     # read in file (real-time) containing list of CpG sites and corresponding genomic positions
		for line in ins:
			if line.split()[1].isdigit():
				if int(line.split()[1]) >= m_min_range and int(line.split()[1]) <= m_max_range: # find meth probes within chosen window of SNP
					meth_probes.append(line.split()[0])

	with open(filename3, "r") as ins:
		for line in ins:
			if line.split()[1].isdigit():
				this_mean = (int(line.split()[1]) + int(line.split()[2])) / 2
				if this_mean >= a_min_range and this_mean <= a_max_range:                       # find acetylation probes within chosen window of SNP
					acetyl_peaks.append(line.split()[0])

	if len(meth_probes)==0:       # skip snp if no meth probes in window (for now)
		continue

	# determine meth probe most correlated to that SNP

	probe_scores = []
	meth_probes2 = []
	good_pts1 = []
	meth_pts1 = []
	meth_pts2 = []
	meth_vector2 = []
	intersect = []
	final_meth_vector = []

	z = 0
	for probe in meth_probes:
		try:
			mIndex = myProbeList.index(probe)
		except:
			continue					     # skip probe if it is not in data matrix
		
		meth_probes2.append(probe)				     # keep track of probes that will have correlation calculated	
		i = -1
		with open (filename6, "r") as ins:
			for line in ins:
				i += 1
				if i == 0 and z == 0:			     # keep list of patients
					good_pts2 = []
					for pt in line.split():
						meth_pts1.append(pt)
						if pt in pt_map1:
							good_pts1.append(pt)		         	# keep list of pts that can be mapped
							good_pts2.append(pt_map2[pt_map1.index(pt)])    # mapping is done here
					for pt in good_pts2:
						tmp = pt.replace('ROS','')				# process patient IDs (removing 'ROS' and 'MAP')
						meth_pts2.append(tmp.replace('MAP', ''))
						intersect = [val for val in meth_pts2 if val in geno_pts]  # determine intersection between methylation and genotype patients
				if i == mIndex:
					meth_vector1 = line.split()	     # construct methylation data vector
					break
		z += 1

		for pt in intersect:
			meth_vector2.append(float(meth_vector1[meth_pts1.index(good_pts1[meth_pts2.index(pt)])]))	# filter the methylation data vector for intersect pts
			geno_vector2.append(float(geno_vector[geno_pts.index(pt)]))
		
		probe_scores.append(abs(pearsonr(geno_vector2, meth_vector2)[0]))
	
	if len(meth_probes2)>=1:						     # requires at least on methylation probe to be in the window and included in data matrix
		bestprobe = meth_probes2[probe_scores.index(max(probe_scores))]
		mIndex = myProbeList.index(bestprobe)

		with open(filename6, "r") as ins:				     # store data vector for most correlated probe
			i = -1
			for line in ins:
				i += 1
				if i == mIndex:
					j = 0
					for val in line.split():
						if j >= 1:			     # skip first val, which is the rowname (probeID) of the data matrix
							final_meth_vector.append(val)
						j += 1

		mFile = this_snp + "_methylationvector.txt"                  # print methylation vector to file
		methFile = open(mFile,'w+')
		methFile.write("\n".join(final_meth_vector))
		methFile.close()

	else:
		print("Warning: no probes to correlate!")

	#if len(acetyl_peaks)>=1:
	#	aFile = this_snp + "_acetylationpeaks.txt"                   # print acetylation peaks to file
	#	acetylFile = open(aFile,'w+')
	#	acetylFile.write("\n".join(acetyl_peaks))
	#	acetylFile.close()

