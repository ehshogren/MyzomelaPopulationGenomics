#!/bin/perl
use strict;
use warnings;
#use IO::File;

# parses a multi-sample gVCF file and, for two groups of individuals (as identified by a key),
# returns two tab-delineated text files, with CHROM\tPOS of 1) fixed differences and 2) shared
# polymorphisms. Also requires as input the minimum number of individuals in each group needed
# in order to make the call (fixed/poly), and a prefix for output files.

# Usage: perl ParseVCF_fixed_shared.pl MtR_chr13_hmd_AlloIndvs_varsite.vcf card_allo,tris_allo myzo_indiv_key.txt 5 chr13_allo


my $sample_vcf = $ARGV[0]; #input vcf (can be gzipped)
my $pop_name_list = $ARGV[1]; #comma-separated list of groups to compare (card_allo,tris_allo))
my $indiv_key = $ARGV[2]; #tab-delinated text file, each line giving the name and group for a single individual
my $min_indiv = $ARGV[3]; #minimum number of individuals in a group needed to call a site "fixed"
my $out_prefix = $ARGV[4]; #prefix of all output files



my @pops = split(",",$pop_name_list);
my $pop1 = $pops[0];
my $pop2 = $pops[1];

my $outfile_fixed = $out_prefix.'_fixed.txt';
my $outfile_shared = $out_prefix.'_shared.txt';
my $outfile_private1 = $out_prefix.'_private1_'.$pop1.'-'.$pop2.'.txt';
my $outfile_private2 = $out_prefix.'_private2_'.$pop1.'-'.$pop2.'.txt';

my %indiv_pop_hash = ();
my %indiv_nums_pop1 = ();
my %indiv_nums_pop2 = ();

open(OUTFILE_FIXED, ">$outfile_fixed")|| print ("Can't open $outfile_fixed!\n");
print OUTFILE_FIXED "#CHROM\tPOS\n";
close(OUTFILE_FIXED);

open(OUTFILE_SHARED, ">$outfile_shared")|| print ("Can't open $outfile_shared!\n");
print OUTFILE_SHARED "#CHROM\tPOS\n";
close(OUTFILE_SHARED);

open(OUTFILE_PRIVATE1, ">$outfile_private1")|| print ("Can't open $outfile_private1!\n");
print OUTFILE_PRIVATE1 "#CHROM\tPOS\n";
close(OUTFILE_PRIVATE1);

open(OUTFILE_PRIVATE2, ">$outfile_private2")|| print ("Can't open $outfile_private2!\n");
print OUTFILE_PRIVATE2 "#CHROM\tPOS\n";
close(OUTFILE_PRIVATE2);


open(FILE, "$indiv_key") || print ("Can't open $indiv_key!\n");
foreach my $line(<FILE>)
{
    chomp($line);
    my @indiv_info = split("\t",$line);
    $indiv_pop_hash{$indiv_info[0]} = $indiv_info[1];
}
close(FILE);

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
            #get the columns for individuals in the two pops
            if ($vcf_info[0] eq '#CHROM'){
                for (my $i = 9; $i < @vcf_info; $i++){
                    if (exists $indiv_pop_hash{$vcf_info[$i]}){
                        if ($indiv_pop_hash{$vcf_info[$i]} eq $pop1){
                            $indiv_nums_pop1{$i} = $vcf_info[$i];
                        }
                        else {
                            if ($indiv_pop_hash{$vcf_info[$i]} eq $pop2){
                                $indiv_nums_pop2{$i} = $vcf_info[$i];
                            }
                        }
                    }
                }
            }
            #check for simple (no adjacent indels) bi-allelic sites
            if ((length( $vcf_info[4] ) == 1) && ($vcf_info[4] ne '*')){
                #count ref (0) and alt (1) alleles in each pop
                my $pop1_count_0 = 0;
                my $pop1_count_1 = 0;
                my $pop2_count_0 = 0;
                my $pop2_count_1 = 0;
                my $tot_pop1 = 0;
                my $tot_pop2 = 0;
                 
                for (my $i = 9; $i < @vcf_info; $i++){
                    if (exists $indiv_nums_pop1{$i}) {
                        if (substr($vcf_info[$i],0,1) eq '0') {
                            $pop1_count_0++;
                        }
                        if (substr($vcf_info[$i],0,1) eq '1') {
                            $pop1_count_1++;
                        }
                        if (substr($vcf_info[$i],2,1) eq '0') {
                            $pop1_count_0++;
                        }
                        if (substr($vcf_info[$i],2,1) eq '1') {
                            $pop1_count_1++;
                        }

                    }
                    if (exists $indiv_nums_pop2{$i}) {
                        if (substr($vcf_info[$i],0,1) eq '0') {
                            $pop2_count_0++;
                        }
                        if (substr($vcf_info[$i],0,1) eq '1') {
                            $pop2_count_1++;
                        }

                        if (substr($vcf_info[$i],2,1) eq '0') {
                            $pop2_count_0++;
                        }
                        if (substr($vcf_info[$i],2,1) eq '1') {
                            $pop2_count_1++;
                        }
                    }
                }
                
                $tot_pop1 = $pop1_count_0 + $pop1_count_1;
                $tot_pop2 = $pop2_count_0 + $pop2_count_1;

                if (($tot_pop1 >= $min_indiv * 2) && ($tot_pop2 >= $min_indiv * 2)) {

                    if (($pop1_count_0 > 0) && ($pop1_count_1 > 0) && ($pop2_count_0 > 0) && ($pop2_count_1 > 0)) {
                        open(OUTFILE_SHARED, ">>$outfile_shared");
                        print OUTFILE_SHARED $vcf_info[0],"\t",$vcf_info[1],"\n";
                        close(OUTFILE_SHARED);
                   }
                    if ((($pop1_count_0 == 0) && ($pop2_count_1 == 0)) || (($pop1_count_1 == 0) && ($pop2_count_0 == 0))) {
                        open(OUTFILE_FIXED, ">>$outfile_fixed");
                        print OUTFILE_FIXED $vcf_info[0],"\t",$vcf_info[1],"\n";
                        close(OUTFILE_FIXED);
                    }
                    if (($pop1_count_0 > 0) && ($pop1_count_1 > 0)) {
                        if (($pop2_count_0 == 0) || ($pop2_count_1 == 0)) {
                            open(OUTFILE_PRIVATE1, ">>$outfile_private1");
                            print OUTFILE_PRIVATE1 $vcf_info[0],"\t",$vcf_info[1],"\n";
                            close(OUTFILE_PRIVATE1);
                        }
                    }
                    if (($pop2_count_0 > 0) && ($pop2_count_1 > 0)) {
                        if (($pop1_count_0 == 0) || ($pop1_count_1 == 0)) {
                            open(OUTFILE_PRIVATE2, ">>$outfile_private2");
                            print OUTFILE_PRIVATE2 $vcf_info[0],"\t",$vcf_info[1],"\n";
                            close(OUTFILE_PRIVATE2);
                        }
                    }
                }
            }
        }
}

close(INFILE);

