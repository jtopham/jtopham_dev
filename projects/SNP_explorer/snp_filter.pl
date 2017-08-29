#!\strawberry\perl\bin

## the purpose of this script is to filter the large (11GB) SNP annotation file for the
## 600k SNPs that were included in the Affymetrix genotype assay, using a hash table.
## program runs out of memory when attempting to make one large hash table, so annotation
## data was subsetted into 3 parts - a hash is constructed and deleted on each subset


use strict;
use warnings;


my %hash = ();

# header for output, as per http://ucscbrowser.genap.ca/cgi-bin/hgTables?hgsid=1187568_HiKf4SC849x9IecdaxddpCBaaQEX&hgta_doSchemaDb=hg19&hgta_doSchemaTable=snp141

my $output = "bin\tchrom\tchromStart\tchromEnd\tname\tscore\tstrand\trefNCBI\trefUCSC\tobserved\tmolType\tclass\tvalid\tavHet\tavHetSE\tfunc\tlocType\tweight\texceptions\tsubmitterCount\tsubmitters\talleleFreqCount\talleles\talleleNs\talleleFreqs\tbitfields\n";

# filepaths for 3 files that are each a unique subset of the large 
# list of annotated SNPs (64 million, from dbSNP)
# and also list of SNPs to search for (from assay)

my $dir_annotated_SNPs1 = '/Users/james/Desktop/jtopham_prj/data/rosmap/snp_annotation/snp141_1.txt';
my $dir_annotated_SNPs2 = '/Users/james/Desktop/jtopham_prj/data/rosmap/snp_annotation/snp141_2.txt';
my $dir_annotated_SNPs3 = '/Users/james/Desktop/jtopham_prj/data/rosmap/snp_annotation/snp141_3.txt';

my @file_list = ($dir_annotated_SNPs1, $dir_annotated_SNPs2, $dir_annotated_SNPs3);

my $dir_assay_SNPs = '/Users/james/Desktop/jtopham_prj/data/rosmap/assayed_snps.txt';


# construct three hash tables, using genotyped-SNP list and the three
# subsets from the large annotation list. after each hash, mutual SNPs
# are recorded and hash is cleared to save on memory


for (my $i=0; $i <= (scalar @file_list); $i++){

    # insert all annotated SNPs (\t-delimited) into hash table
    
    open(my $SNP_list, '<:encoding(UTF-8)', $file_list[$i])
        or die "Could not open file '$file_list[$i]' $!";

    while (my $this_row = <$SNP_list>) {
        my @this_row_elements = split(/\t/, $this_row);
        my $this_SNP = $this_row_elements[4];
        chomp $this_SNP;
        $hash{$this_SNP}{"annotation"} = $this_row;
    }

    # insert assayed SNPs (\n-delimited) into hash table

    open(my $assay_list, '<:encoding(UTF-8)', $dir_assay_SNPs)
        or die "Could not open file '$dir_assay_SNPs' $!";

    while (my $this_row = <$assay_list>) {
        chomp $this_row;
        $hash{$this_row}{"assay"} = "present";
    }

    # determine which SNPs are present in both files

    foreach my $snp (keys %hash){
        while (my ($origin, $value) = each %{ $hash{$snp} } ) {
            if (keys %{ $hash{$snp} } == 2){
                $output = $output . $hash{$snp}{"annotation"};
            }
            last;
        }
    }
    
    # clear hash

    for (keys %hash){delete $hash{$_};}

}


# output to file

my $newFile = "annotated_SNPs_filtered.tsv";
open(my $fh, '>', $newFile);
print $fh $output;
close $fh;

