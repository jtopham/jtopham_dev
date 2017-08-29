#!/bin/bash

## This script functions to streamline performing featureCount on multiple bash files 
## within the current directory

# local:

#export PATH=/home/jtopham/Software/subread/subread-1.4.6-p5-source/bin:$PATH

# marrahost01:

export PATH=/subread-1.4.6-p5-source/bin:$PATH

GTF=/projects/jtopham_prj/Thesis_prj/data/RNA_seq/Homo_sapiens.GRCh37.75.gtf
EXPTNAME=HEK293
CPUS=8    # number of threads
#MAPQ=10   # minimum mapping quality score a read must satisfy in order to be counted
GENEMX=${EXPTNAME}_featCount.txt

# -p specifies paired-end reads, so that paired reads are not counted twice
# -B specifies to NOT count fragments (ie. discard singletons)
# -C specifies to NOT count chimeric reads

featureCounts -pC -T $CPUS -t 'gene' -a $GTF -o $GENEMX /projects/jtopham_prj/Thesis_prj/data/RNA_seq/*.bam


