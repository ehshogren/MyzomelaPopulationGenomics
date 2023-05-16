# Gene flow in recently sympatric honeyeaters with neo-sex chromosomes
This is an overview of the analyses and code used in Shogren et al. to understand the evolutionary context and degree and direction of introgression between *Myzomela cardinalis* and *Myzomela tristrami*. 

The software and programs used throughout the pipeline include:
- trimgalore
	- FastQC
	- cutadapt
	- pigz
	- multiQC
- bwa
- samtools 
- GATK (v4.1.7.0)
	- picard (v2.12.0)
- qualimap
- PSMC
- vcf2phylip
- popArt
- python3
- vcftools
- bcftools
- pixy
- perl
- plink1.9
- ADMIXTURE
- R (v4.1.1)
	- snpRelate
	- tidyverse
- Dsuite

## Contents
- raw read processing and alignment
- aligned read filtering 
- calling variants 
- filtering variants
- PSMC
- mtDNA
- pi, dxy, Fst, Tajima's D
- private/shared/fixed alleles
- PCA
- ADMIXTURE
- ABBA-BABA

## raw read processing and alignment
samples were sequenced over 1-3 lanes. To maximize parallelization, I ran analyses in batches of about 20 individuals, waiting to concatenate files for individuals until later in the pipeline. Thus, many of the processing scripts loop over each sample/cell/lane ID. 

### trimming reads, fastQC
Starting with removing adapters using trimgalore, and generating a fastQC report. Adapter sequences were provided by the sequencer. Programs needed for this step include trimgalore, fastQC, cutadapt, and pigz

```
sample_names_path=/CardAllo_scl.txt
reads_path=/raw_data/
fastqc_out=/FastQC_trimmed
out_path=/trimgalore


cat $sample_names_path |
while read -r sample cell lane; do 
trim_galore --paired --retain_unpaired --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --stringency 5 --phred33 --fastqc_args "--outdir $fastqc_out" --basename ${sample}_${cell}_${lane} --cores 12 --output_dir $out_path $reads_path/$sample/${sample}_*_${cell}_${lane}_1.fq.gz $reads_path/$sample/${sample}_*_${cell}_${lane}_2.fq.gz ;
done
```

Once FastQC has generated a report for each fastq file, can run multiqc in the directory with the fastQC output files. It will automatically detect the relevant files and generate a summary report for the dataset. Just navigate to the fastqc_out directory and run: 
```
multiqc .
```
### bwa-mem + samtools sort
align trimmed reads to the reference genome using bwa-mem, then immediately sort and convert to bam files using samtools. For this step need bwa and samtools

First, need to create a bwa database index for the reference genome
```
bwa index /Mt_v1.0_MAIN.fa
```
Once index files are generated and in the same directory as the reference genome file, can continue with the alignment and sorting step.
```
sample_names_path=/CardAllo_scl.txt
database_path=/Mt_v1.0_MAIN.fa
trimmed_reads_path=/trimgalore
out_path=/aligned_sorted

cat $sample_names_path | 
while read -r sample cell lane; do 
bwa mem -t 12 -M $database_path $trimmed_reads_path/${sample}_${cell}_${lane}_R1_val_1.fq.gz $trimmed_reads_path/${sample}_${cell}_${lane}_R2_val_2.fq.gz | samtools view -bS | samtools sort --threads 12 -o $out_path/${sample}_${cell}_${lane}_sorted.bam;
done
``` 
## aligned read filtering
bam files will now be cleaned up prior to calling variants, using GATK and picard.

### add read groups
adding information on read groups to presumably help the algorithm later on in the variant calling pipeline
Included RGID (sampleID), RGPU (cell.lane), RSM (sampleID again), RGPL (Illumina) for all samples.
For RGLB, first sequencing run = MyzoWGS1, second sequencing run = MyzoWGS2, and M. pulchella = KU 
```
sample_names_path=/CardAllo_scl.txt
sorted_bam_path=/aligned_sorted
rg_bam_path=/RG_added

cat $sample_names_path | 
while read sample cell lane; do 
gatk AddOrReplaceReadGroups -I $sorted_bam_path/${sample}_${cell}_${lane}_sorted.bam -O $rg_bam_path/${sample}_${cell}_${lane}_sorted_rg.bam -RGID ${sample} -RGPU ${cell}.${lane} -RGSM ${sample} -RGPL Illumina -RGLB MyzoWGS1 -SO coordinate;
done
```
### mark duplicates
at this point in addition to marking duplicates, will combine samples that are split across multiple bam files so that we have a single bam per individual. To accomplish this, group and run samples based on the number of lanes they were split across, and provide the appropriate number of INPUT arguments to concatenate. That is, if an individual had two fastq files (and subsequently two bam files) will include two INPUT arguments that list each of those bam files. Output will be a single <sample>.bam file. 
```
sample_names_path=/Mt_v1_MarkDup_TwoLnSamples_2.txt
rg_bam_path=/RG_added
marked_bam_path=/marked_dup
metrics_path=/marked_dup_metrics

cat $sample_names_path |
while read sample cell1 lane1 cell2 lane2; do 
java -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
INPUT=$rg_bam_path/${sample}_${cell1}_${lane1}_sorted_rg.bam \
INPUT=$rg_bam_path/${sample}_${cell2}_${lane2}_sorted_rg.bam \
OUTPUT=$marked_bam_path/${sample}_sorted_rg_marked.bam \
METRICS_FILE=$metrics_path/${sample}_metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000;
done
```
### fixmate information
although there is now a single bam file per individual, still parallelizing by splitting up into approximately 20 sample groups, so switching from sample/cell/lane (scl.txt) files to just sampleID.txt files
```
sample_names_path=/CardAllo_SampleID.txt
marked_bam_path=/marked_dup
fixmate_bam_path=/fixed

cat $sample_names_path |
while read sample; do 
gatk FixMateInformation -INPUT $marked_bam_path/${sample}_sorted_rg_marked.bam -OUTPUT $fixmate_bam_path/${sample}_sorted_rg_marked_fixmate.bam -SO coordinate -CREATE_INDEX true;
done
```
### Qualimap
generates a report for each sample with statistics for mapping, read numbers, duplicates, etc. Requires program qualimap
```
sample_names_path=/CardAllo_SampleID.txt
bam_path=/fixed
out_path=/qualimap


cat $sample_names_path | 
while read sample; do 
qualimap bamqc -bam $bam_path/${sample}_sorted_rg_marked_fixmate.bam -sd -c -nw 400 -hm 3 -outdir $out_path/${sample} --java-mem-size=120G;
done
```
once you have run qualimap for all samples, you can generate a multi-sample qualimap
first, need a sample file that is a tab-separated, with columns for sample ID, path to qualimap output for that sample, and the species/sampling locality (pop) of the sample.
```
SAMPLEPOPFILE=/Mt_v1_SampleID_pop.txt
QUALIMAP=/qualimap/

cat ${SAMPLEPOPFILE} | while read sampleid pop; do 
echo -e ${sampleid}'\t'${QUALIMAP}${sampleid}'\t'${pop};
done > multibam_qualimap_samplefile.txt
```
using multibam_qualimap_samplefile.txt, can run multibamqc
```
sample_names_path=/multibam_qualimap_samplefile.txt
out_path=/multibam_qualimap

qualimap multi-bamqc -d $sample_names_path -outdir $out_path
```
## calling variants
Using GATK4 to call variants using the now cleaned up bam files. Will need GATK, picard, and samtools

### reference index
need to generate both faidx and dict files for the reference genome, leaving them in the same directory as the reference genome fasta file
for the .dict file, need GATK and picard
```
java -jar /software/picard/2.12.0/picard.jar CreateSequenceDictionary \
R=/Mt_v1.0_MAIN.fa \
O=/Mt_v1.0_MAIN.dict
```
for the .fai files need samtools (I used v1.7)
```
samtools faidx /Mt_v1.0_MAIN.fa
```
### repeat masking
to speed up the haplotype call step (or at least make it not as slow), use a bed file with coordinates of repetitive regions, generated using RepeatMasker. 
may need to change extension of file to .bed and remove header, then index so GATK can work with it
```
gatk --java-options "-Xmx4g -Xms4g" IndexFeatureFile --input /Mt_v1.0_MAIN_RM_sites_to_filter.bed
```
### HaplotypeCaller
This generates a per-individual g.vcf.gz file, preliminary to genotyping/calling variants across multiple samples.
Because we are particularly interested in sex chromosome variation and want to keep coverage consistent for the algorithm, calling autosomal and pseudo-autosomal regions (diploid in all individuals) separately from sex chromosomes. For sex chromosomes, running haplotype caller separately on males (diploid for Z) and females (haploid for Z and W sequence). This is accomplished using interval lists (-L argument) and sampleID files separated by sex, based on coverage of Z sequence from qualimap (females with half coverage on Z sequence). However, we maintain diploid (default) settings for all runs. 

autosomes and neo-pseudoautosomal region (neoPAR) haplotype call
```
sample_names_path=/CardAllo_SampleID.txt
reference_path=/Mt_v1.0_MAIN.fa
fixmate_bam_path=/fixed
out_path=/haplotype_call

cat $sample_names_path |
while read sample; do 
gatk --java-options "-Xmx64G" HaplotypeCaller \
-R $reference_path \
-I $fixmate_bam_path/${sample}_sorted_rg_marked_fixmate.bam \
-O $out_path/${sample}_auto_neoPAR.g.vcf.gz \
-L /IntervalList_files/Mt_v1.0_MAIN_autos_neoPAR_interval.list \
-XL /Mt_v1.0_MAIN_RM_sites_to_filter.bed \
-ERC GVCF;
done
```
female sex chromosomes haplotype call
```
sample_names_path=/Mt_v1_females.txt
reference_path=/Mt_v1.0_MAIN.fa
fixmate_bam_path=/fixed
out_path=/haplotype_call

cat $sample_names_path |
while read sample; do 
gatk --java-options "-Xmx64G" HaplotypeCaller \
-R $reference_path \
-I $fixmate_bam_path/${sample}_sorted_rg_marked_fixmate.bam \
-O $out_path/${sample}_sexchr.g.vcf.gz \
-L /chrZ1_Z2_W1_W2_exclPAR_interval.list \
-XL /Mt_v1.0_MAIN_RM_sites_to_filter.bed \
-ERC GVCF;
done
```
male sex chromosomes haplotype call
```
sample_names_path=/Mt_v1_males.txt
reference_path=/Mt_v1.0_MAIN.fa
fixmate_bam_path=/fixed
out_path=/haplotype_call

cat $sample_names_path |
while read sample; do 
gatk --java-options "-Xmx64G" HaplotypeCaller \
-R $reference_path \
-I $fixmate_bam_path/${sample}_sorted_rg_marked_fixmate.bam \
-O $out_path/${sample}_sexchr.g.vcf.gz \
-L /chrZ1_Z2_interval.list \
-XL /Mt_v1.0_MAIN_RM_sites_to_filter.bed \
-ERC GVCF;
done
```
### combine GVCFs
combines individual g.vcfs to genomic region g.vcfs that include all individuals sequenced in that region. Splitting genome up using interval lists again to speed things along. 
need a sample map that is simply a .list file with a path to the sample.g.vcf.gz files on each line.
```
SAMPLEPOPFILE=/Mt_v1_SampleID_pop.txt
GVCF=/haplotype_call/

cat ${SAMPLEPOPFILE} | while read sampleid pop; do 
echo -e ${GVCF}${sampleid}_auto_neoPAR.g.vcf.gz;
done > Mt_v1_genomicsDB_samplemap_auto_neoPAR.list
```
using sample map, run CombineGVCFs using GATK and picard
```
gatk --java-options "-Xmx84g" CombineGVCFs \
-R /Mt_v1.0_MAIN.fa \
-L chr1 \
-XL /Mt_v1.0_MAIN_RM_sites_to_filter.bed \
-V /Mt_v1_genomicsDB_samplemap_auto_neoPAR.list \
-O /Mt_v1_chr1_cohort.g.vcf.gz
```
### genotype GVCFs
after combining to region-specific g.vcf.gzs, genotype to generate an unfiltered, all sites vcf.gz file
```
gatk --java-options "-Xmx64g" GenotypeGVCFs \
-R /Mt_v1.0_MAIN.fa \
-all-sites -L chr1 \
-XL /Mt_v1.0_MAIN_RM_sites_to_filter.bed \
-V /Mt_v1_chr1_cohort.g.vcf.gz \
-O /Mt_v1_chr1_unfilt_allsites.vcf.gz
```
## filtering variants
