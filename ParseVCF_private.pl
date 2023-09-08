#!/bin/perl
use strict;
use warnings;

# parses a multi-sample gVCF file and, given a focal species with several locations (as identified by a key),
# plus at least one other species of interest, returns a tab-delineated text file, with CHROM\tPOS for each site
# where there is a private ALLELE in the focal SPECIES not present in the others (excluding putative hybrids),
# the identity of that allele, followed by one column for each location giving a Y/N for presence/absence of that
# private allele.
#
# The assumed format of the "key" file corresponds to the Myzo_WGS metadata files (full metadata), with the
# (currently!) important columns the 3rd (sample location), 4th (species) and 7th (sample name).
#
# Also requires as input the minimum number of individuals in each species (note: NOT location sample) needed in
# order to make the call as to whether an allele is private.
#
# As written, the program ignores any individual with the "species" name "Hybrid" in the key file. Name accordingly!
# Alternatively, if you give the program an incomplete key, it will ignore any sample in the vcf that it has not
# already seen in the key, so you could simply delete all hybrids (and possibly the "MO" cardinalis sample) from
# the key being used.
#
# PRINTS TO STDOUT!!!!!!!!

# Usage: perl ParseVCF_private.pl MtR_chr13_hmd_varsite.vcf Cardinalis Myzo_WGS_meta_corrected.csv 5 > private_card_chr13.txt


my $sample_vcf = $ARGV[0]; #input vcf (can be gzipped)
my $focal_species = $ARGV[1]; # name of focal species (must match format of the key)
my $indiv_key = $ARGV[2]; # csv file, formatted as Myzo metadata, with the location, species, and sample name in the 3rd,4th,7th columns
my $min_indiv = $ARGV[3]; # minimum number of individuals in a group needed to make a call re presence/absence

# we will get these from the key file
my %loc_hash = ();

my %focalsp_indiv_hash = (); # holds the nqmes individuals in the focal species
my %nonfocalsp_indiv_hash = ();
my %focalsp_nums = (); # holds the indices of the vcf for individuals in the focal species
my %nonfocalsp_nums = ();
my %focalsp_loc_nums = (); # a hash of hashes, one for each location, each similar to focalsp_nums


open(FILE, "$indiv_key") || print ("Can't open $indiv_key!\n");
foreach my $line(<FILE>)
{
    chomp($line);
    my @indiv_info = split(",",$line);
    if ($indiv_info[3] eq $focal_species){
        $focalsp_indiv_hash{$indiv_info[6]} = $indiv_info[2];
        $loc_hash{$indiv_info[2]} = 1;
    }
    else {
        if (substr($indiv_info[3],0,6) ne 'Hybrid'){
            $nonfocalsp_indiv_hash{$indiv_info[6]} = $indiv_info[2];
        }
    }
}
close(FILE);

# 'convert' location hash to array; for consistent indexing order in output

my @loc_list = ();
foreach my $key (keys %loc_hash){
    push @loc_list, $key;
    $focalsp_loc_nums{$key} = ();
}
my $num_loc = @loc_list;


print "CHROM\tPOS\t",$focal_species,"_private_allele";
for (my $i = 0; $i < $num_loc; $i++){
    print "\t",$loc_list[$i];
}
print "\n";

#open(INFILE, "$sample_vcf")|| print ("Can't open $sample_vcf!\n");

if ($sample_vcf =~ /.gz$/) {
    open(INFILE, "gunzip -c $sample_vcf |") || die "can't open $sample_vcf";
}
else {
    open(INFILE, $sample_vcf) || die "can't open $sample_vcf";
}

foreach my $line(<INFILE>)
{
    chomp($line);
    if (substr($line,0,2) ne '##'){
        my @vcf_info = split("\t", $line);
        #get the column indices for individuals in the focal and non-focal species. EXCLUDE KNOWN HYBRIDS
        if ($vcf_info[0] eq '#CHROM'){
            for (my $i = 9; $i < @vcf_info; $i++){
                if (exists $focalsp_indiv_hash{$vcf_info[$i]}){
                    $focalsp_nums{$i} = $vcf_info[$i];
                    my $cur_loc = $focalsp_indiv_hash{$vcf_info[$i]};
                    $focalsp_loc_nums{$cur_loc}{$i} = 1;
                }
                else {
                    if (exists $nonfocalsp_indiv_hash{$vcf_info[$i]}){
                        $nonfocalsp_nums{$i} = $vcf_info[$i];
                    }
                }
            }
        }
        else{
            #check for simple (no adjacent indels) bi-allelic sites
            if ((length( $vcf_info[4] ) == 1) && ($vcf_info[4] ne '*')){
                #count ref (0) and alt (1) alleles in each group
                my $focal_count_0 = 0;
                my $focal_count_1 = 0;
                
                my $non_focal_count_0 = 0;
                my $non_focal_count_1 = 0;
                
                my $tot_focal = 0;
                my $tot_non_focal = 0;
                
                my $ident_private = 'NA';
                my $nuc_private = 'N';
                
                for (my $i = 9; $i < @vcf_info; $i++){
                    if (exists $focalsp_nums{$i}) {
                        if (substr($vcf_info[$i],0,1) eq '0') {
                            $focal_count_0++;
                        }
                        if (substr($vcf_info[$i],0,1) eq '1') {
                            $focal_count_1++;
                        }
                        if (substr($vcf_info[$i],2,1) eq '0') {
                            $focal_count_0++;
                        }
                        if (substr($vcf_info[$i],2,1) eq '1') {
                            $focal_count_1++;
                        }
                        
                    }
                    if (exists $nonfocalsp_nums{$i}) {
                        if (substr($vcf_info[$i],0,1) eq '0') {
                            $non_focal_count_0++;
                        }
                        if (substr($vcf_info[$i],0,1) eq '1') {
                            $non_focal_count_1++;
                        }
                        if (substr($vcf_info[$i],2,1) eq '0') {
                            $non_focal_count_0++;
                        }
                        if (substr($vcf_info[$i],2,1) eq '1') {
                            $non_focal_count_1++;
                        }
                    }
                }
                
                $tot_focal = $focal_count_0 + $focal_count_1;
                $tot_non_focal = $non_focal_count_0 + $non_focal_count_1;
                
                if (($tot_focal >= $min_indiv * 2) && ($tot_non_focal >= $min_indiv * 2)) {
                    if (($focal_count_0 > 0) && ($non_focal_count_0 == 0)) {
                        $ident_private = 0;
                        $nuc_private = $vcf_info[3];
                    }
                    if (($focal_count_1 > 0) && ($non_focal_count_1 == 0)){
                        $ident_private = 1;
                        $nuc_private = $vcf_info[4];
                    }
                }
                # go through each location being considered, check if the private allele is present
                if ($ident_private ne 'NA'){
                    print $vcf_info[0],"\t",$vcf_info[1],"\t",$nuc_private;
                    for (my $i = 0; $i < $num_loc; $i++){
                        my $is_present = 'N';
                        my $cur_loc = $loc_list[$i];
                        foreach my $index (keys %{$focalsp_loc_nums{$cur_loc}}){
                            my $cur_geno = $vcf_info[$index];
                            if ((substr($cur_geno,0,1) eq $ident_private) || (substr($cur_geno,2,1) eq $ident_private)){
                                $is_present = 'Y';
                            }
                        }
                        print "\t",$is_present;
                    }
                    print "\n";
                }
            }
        }
    }
}
close(INFILE);
