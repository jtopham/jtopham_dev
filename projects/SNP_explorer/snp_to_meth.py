#!/usr/bin/env python

import sys

filename = "/home/unix/jtopham/data/snp_mapper/snpArrayAffy6.txt"   # maps SNP ids to genomic position
filename2 = "/home/unix/jtopham/data/methylation/filtered_meth_mapper.txt"   # maps CpG sites to genomic position

snp_mapper = []

with open(filename, "r") as ins:     # read in file contining list of SNP ids and their corresponding genomic positions
    for line in ins:
        snp_mapper.append(line.strip('\n'))

args = sys.stdin.readlines()         # this is a list of \n-delim SNPs supplied by stdin ( ie. more snplist.txt | python thisscript.pl )

for arg in args:
    this_snp = arg.strip('\n')
        
    this_position = []
        
    for line in snp_mapper:    # identify position of this SNP
        if line.split()[8] == this_snp:
            this_position.append(line.split()[3])
    
    if len(this_position)>1:
        print('Warning:', this_snp, 'maps to more than one position!')
                
    if len(this_position)==0:
        #print('Warning:', this_snp, 'does not map to any position')
        continue
                    
    min_range = int(this_position[0]) - 1000     # use 2kb window
    max_range = int(this_position[0]) + 1000
                            
    if min_range < 0:          # ensure min_range doesn't fall below zero (likely unnecessary)
        min_range = 0
                                    
    meth_probes = []
                                        
    with open(filename2, "r") as ins:     # read in file (real-time) containing list of CpG sites and corresponding genomic positions
        for line in ins:
            if line.split()[1].isdigit():
            if int(line.split()[1]) >= min_range and int(line.split()[1]) <= max_range: # find meth probes within 2kb window of SNP
                meth_probes.append(line.split()[0])
                                                            
                                                            
    for probe in meth_probes:  # print list of methylation probes to std out
        print(probe)


#this part will be moved to different file later - input to this probably ("for each snp"-> function here)

filename = "/Users/james/Desktop/jtopham_prj/Data/mock_geno_data.txt"   # contains list of SNPs

snp_list = []
stop = 1

with open(filename, "r") as ins:
    for line in ins:
        snp_list.append(line.strip('\t'))
        stop += 1
        if stop==2:      # take only header
            break

snp_list2 = snp_list[0]

snp_list_split = snp_list2.split()  # this contains 3197 SNP ids, some of which are false ie. 'chr9: 281234'

#print('\n')
#print(snp_list_split[0])