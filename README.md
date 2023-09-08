# MyzomelaPopulationGenomics
Pipeline and scripts associated with population genomic study of Myzomela honeyeaters
# Gene flow in recently sympatric honeyeaters with neo-sex chromosomes
This is an overview of the analyses and code used in Shogren et al. to understand the evolutionary context and degree and direction of introgression between *Myzomela cardinalis* and *Myzomela tristrami*. 

The software and programs used throughout the pipeline include:
- [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
	- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
	- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
	- [pigz](https://zlib.net/pigz/)
	- [multiQC](https://multiqc.info/docs/)
- [bwa](https://bio-bwa.sourceforge.net/)
- [samtools/bcftools](http://samtools.github.io/bcftools/bcftools.html) 
- [GATK](https://gatk.broadinstitute.org/hc/en-us) (v4.1.7.0)
	- picard (v2.12.0)
- [qualimap](http://qualimap.conesalab.org/)
- [PSMC](https://github.com/lh3/psmc)
- [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)
- [PopART](https://popart.maths.otago.ac.nz/)
- [python3](https://www.python.org/downloads/)
- [vcftools](https://vcftools.github.io/man_latest.html)
- [pixy](https://pixy.readthedocs.io/en/latest/about.html)
- [perl](https://www.perl.org/)
- [plink1.9](https://www.cog-genomics.org/plink/)
- [R](https://www.r-project.org/) (v4.1.1)
	- [SNPRelate](https://www.bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html)
- [ADMIXTURE](https://dalexander.github.io/admixture/index.html)
- [Dsuite](https://github.com/millanek/Dsuite)

## Contents
- [raw read processing and alignment](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#raw-read-)
- [aligned read filtering](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#aligned-read-filtering)
- [calling variants](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#calling-variants)
- [filtering variants](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#filtering-variants)
- [PSMC](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#psmc)
- [mtDNA](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#mtdna)
- [pi, dxy, Fst, Tajima's D](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#pi-dxy-fst-tajimas-d)
- [fixed/shared alleles](https://github.com/ehshogren/MyzomelaPopulationGenomics#fixedshared-alleles)
- [private alleles](https://github.com/ehshogren/MyzomelaPopulationGenomics#private-alleles)
- [PCA](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#pca)
- [ADMIXTURE](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#admixture)
- [triangle plot](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#triangle-plot)
- [ABBA-BABA](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/README.md#abba-baba)

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
at this point in addition to marking duplicates, will combine samples that are split across multiple bam files so that we have a single bam per individual. To accomplish this, group and run samples based on the number of lanes they were split across, and provide the appropriate number of INPUT arguments to concatenate. That is, if an individual had two fastq files (and subsequently two bam files) will include two INPUT arguments that list each of those bam files. Output will be a single sample.bam file. 
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
### unfiltered statistics
prior to filtering, use VCFtools to generate a variant site version of the VCF and look at statistics (average coverage per site, per individual, etc.) which will be used to hone depth filtering parameters later on. Uses vcftools and samtools (v1.7)
```
vcf_all=/Mt_v1_chr1_unfilt_allsites.vcf.gz

vcftools --gzvcf $vcf_all \
--mac 1 \
--recode --recode-INFO-all --stdout | bgzip -c > /Mt_v1_chr1_unfilt_varsite_recode.vcf.gz
tabix /Mt_v1_chr1_unfilt_varsite_recode.vcf.gz

vcf_var=/Mt_v1_chr1_unfilt_varsite_recode.vcf.gz
out=/MtR_chr1_unfilt_varsite
vcftools --gzvcf $vcf_var --freq2 --out $out --max-alleles 2
vcftools --gzvcf $vcf_var --depth --out $out
vcftools --gzvcf $vcf_var --site-mean-depth --out $out
vcftools --gzvcf $vcf_var --site-quality --out $out
vcftools --gzvcf $vcf_var --missing-indv --out $out
vcftools --gzvcf $vcf_var --missing-site --out $out
vcftools --gzvcf $vcf_var --het --out $out
```
### flagging low quality sites
Use GATK's VariantFiltration to flag sites that don't meet  minimum criteria for quality (best practices, hard filtering thresholds, non-model organisms)
```
gatk --java-options "-Xmx36g -Xms36g" VariantFiltration \
-V /Mt_v1_chr1_unfilt_allsites.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-O /Mt_v1_chr1_flag_qual_allsite.vcf.gz
```
### flagging low quality sites and removing spurious heterozygotes on sex chromosomes
Checking for erroneous heterozygote variant calls on female W and Z sequence. Using modified scripts, originally written by Drew Schield for [analyses of Z-linked genomic variation in barn swallows](https://github.com/drewschield/Z-chromosome_analysis_hirundo)
Start by extracting variant sites from the all-site (will need unzipped VCF for the python3 script), and use the list of female individuals in the dataset (in our case identified by mapped coverage on Z chromosome). This will output a bed file, which you can then index for GATK and use to flag sites in VariantFiltration
```
python3 identifyFemaleWhetSites.py /Mt_v1_females.txt /Mt_v1_chrW1_W2_unfilt_varsite_recode.vcf /Female_WhetSites.bed
```
```
gatk --java-options "-Xmx4g -Xms4g" IndexFeatureFile --input /Female_WhetSites.bed
```
```
gatk --java-options "-Xmx84g -Xms84g" VariantFiltration \
-V /Mt_v1_chrW1W2_unfilt_allsites.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-filter "SOR > 3.0" --filter-name "SOR3" \
--mask /FemaleWhetSites.bed --mask-name WHET \
-O /Mt_v1_chrW1W2_flag_qual_het_allsite.vcf.gz
```
### filtering sites for quality and depth
Use bcftools to recode all sites failing filters as missing (./.)
First, include all sites with PASS flag (quality filter)
```
bcftools filter --threads 12 -i 'FILTER="PASS"' \
--set-GTs . -O z -o /Mt_v1_chr1_q_allsite.vcf.gz \
/Mt_v1_chr1_flag_qual_allsite.vcf.gz

tabix -p vcf /Mt_v1_chr1_q_allsite.vcf.gz
```
Next, exclude any sites that fall outside of depth min/max thresholds (depth filter)
```
bcftools filter --threads 6 -e 'MEAN(FORMAT/DP)>34 || MEAN(FORMAT/DP)<10' \
--set-GTs . -O z -o /Mt_v1_chr1_qd_allsite.vcf.gz \
/Mt_v1_chr1_q_allsite.vcf.gz

tabix /Mt_v1_chr1_qd_allsite.vcf.gz
```
This is now the quality and depth filtered dataset used for calculating pi, dxy , Fst, Tajima's D, private/fixed/shared alleles, ABBA-BABA analyses
### filtering for linkage disequilibrium and minor allele frequency
at this point working with only variant sites, so easier to run for all autosomes (used GATK GatherVCFs to combine VCFs). Use vcftools and samtools
```
vcftools --gzvcf /Mt_v1_autos_qd_varsite_recode.vcf.gz \
--maf 0.05 \
--recode --recode-INFO-all --stdout | bgzip -c > /Mt_v1_autos_qd_maf05_varsite_recode.vcf.gz

tabix /Mt_v1_autos_qd_maf05_varsite_recode.vcf.gz
```
LD prune using plink. To accommodate variant naming in plink, first create annotated VCF in bcftools.
```
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' /Mt_v1_autos_qd_maf05_varsite_recode.vcf.gz -O z --output /Mt_v1_autos_qd_maf05_varsite_recode_ann.vcf.gz

plink --vcf /Mt_v1_autos_qd_maf05_varsite_recode_ann.vcf.gz --indep-pairwise 50 5 0.5 --allow-extra-chr --chr-set 80 --out /Mt_v1_autos_qd_maf05_varsite_recode_LD

plink --vcf /Mt_v1_autos_qd_maf05_varsite_recode_ann.vcf.gz --extract /Mt_v1_autos_qd_maf05_varsite_recode_LD.prune.in --make-bed --allow-extra-chr --chr-set 80 --out /Mt_v1_autos_qd_maf05_varsite_recode_LD

plink --bfile /Mt_v1_autos_qd_maf05_varsite_recode_LD --recode vcf --allow-extra-chr --chr-set 80 --out /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned

bgzip /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.vcf

tabix /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.vcf.gz
```
This dataset is now filtered for quality, depth, minor allele frequency minimum of 0.05, and pruned for linkage disequilibrium. Used in PCA, and ADMIXTURE analyses

## PSMC
followed this [tutorial](https://informatics.fas.harvard.edu/psmc-journal-club-walkthrough.html) to generate consensus sequence and run PSMC for autosomes, in only allopatric populations.
### re-calling variants
using bcftools and samtools to call variants from filtered bam files
```
reference_fasta=/Mt_v1.0_MAIN.fa
sampleIDs=/Mt_v1_Ugi.txt

cat $sampleIDs |
while read sample; do
samtools mpileup -Q 30 -q 30 -u -v \
-f $reference_fasta /${sample}_sorted_rg_marked_fixmate.bam |
bcftools call -c --threads 12 |\
bgzip -c > /${sample}_genome_bcftoolscall.vcf.gz;
done
```
tabix bcftoolscall.vcf.gz files using samtools
```
tabix <sample>_genome_bcftools.vcf.gz
```
generate a consensus sequence and restrict to autosome using a regions file
```
sampleIDs=/scratch/juy3_lab/MyzoWGS/Mt_v1_WGS/SampleID_files/Mt_v1_Ugi.txt
regions_file=/scratch/juy3_lab/MTRIS_v1.0/Mt_v1.0_autos.bed

cat $sampleIDs |
while read sample; do
bcftools view --threads 12 --regions-file $regions_file /${sample}_genome_bcftoolscall.vcf.gz | \
vcfutils.pl vcf2fq -d 10 -D 34 -Q 30 > /${sample}_autos_consensus.fq;
done
```
convert consensus sequence to a psmc fasta file
```
psmcfa_script=/psmc/utils/fq2psmcfa
sampleIDs=/scratch/juy3_lab/MyzoWGS/Mt_v1_WGS/SampleID_files/Mt_v1_Ugi.txt

cat $sampleIDs |
while read sample; do
$psmcfa_script /${sample}_autos_consensus.fq > /${sample}_autos.psmcfa;
done
```
run psmc using default settings
```
psmc_script=/psmc/psmc
sampleIDs=/Mt_v1_Ugi.txt

cat $sampleIDs |
while read sample; do
$psmc_script -p "4+25*2+4+6" -o /${sample}_autos.psmc /${sample}_autos.psmcfa;
done
```
modified workflow from [Schumer lab](https://openwetware.org/wiki/Schumer_lab:_Commonly_used_workflows#Demographic_Inference_with_PSMC) for plotting. Used generation times from Bird (2020). 
## mtDNA
filtered mitochondrial sequence slightly differently - only imposing a minimum depth filter (10x) and masking out regions of heteroplasmy
```
vcftools --gzvcf /Mt_v1_chr34_40_mt_neoPAR_q_allsite.vcf.gz \
--bed /mtDNA_include_update.bed \
--min-meanDP 10.0 \
--recode --recode-INFO-all --stdout | bgzip -c > /Mt_v1_mtDNA_MASK_qdmin_allsite.vcf.gz

tabix /Mt_v1_mtDNA_MASK_qdmin_allsite.vcf.gz
```
then, convert vcf file to phylip format using [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py).  Also using the *Myzomela pulchella* individual as an outgroup.
```
vcf2phylip=/vcf2phylip/vcf2phylip.py
input_vcf=/Mt_v1_mtDNA_MASK_qdmin_allsite.vcf.gz

python3 $vcf2phylip --input $input_vcf --outgroup KU8
```
import the resulting phylip file into the [PopART](https://popart.maths.otago.ac.nz/) GUI program, and import a traits file for labeling of nodes (Mt_v1_SampleID_pop_trait_MakUgiThreeOut_relabel.txt)
Construct a TCS network, color nodes according to the traits (populations) and arrange as needed. 

## pi, dxy, Fst, Tajima's D
### pixy
Used [pixy](https://pixy.readthedocs.io/en/latest/about.html) to calculate genome wide summary statistics (pi, dxy, Fst) in 50kb windows. Used a quality, depth filtered all sites VCF, and a tab-delimited populations file, first column with sample ID, second column with sample population.
```
pixy --stats pi fst dxy \
--vcf /Mt_v1_chr1_qd_allsite.vcf.gz \
--populations /Mt_v1_SampleID_pop_card_tris.txt \
--window_size 50000 \
--n_cores 24
```
After running pixy, filtered output by removing windows with fewer than 10,000 sites sequenced. This prevented estimation of statistics for regions that were largely masked/missing data due to repetitive sequence, etc. plotted output in ggplot, modifying developer's example scripts.

### Tajima's D
For each population (--keep argument), estimated Tajima's D in 50kb windows using vcftools.
```
vcftools --gzvcf /Mt_v1_autos_mt_neoPAR_qd_varsite.vcf.gz \
--keep /CardAllo_SampleID.txt \
--TajimaD 50000 \
--out /Mt_v1_autos_mt_neoPAR_50kb_cardUgi
```
## fixed/shared alleles
Used [custom perl script](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/ParseVCF_fixed_shared.pl) to count the number of alleles fixed and shared between populations/species. Then ran for each pairwise comparison of population/species.

Perl script - ParseVCF_fixed_shared.pl - parses a multi-sample gVCF file and, for two groups of individuals (as identified by a key), returns two tab-delineated text files, with CHROM\tPOS of 1) fixed differences and 2) shared polymorphisms. Also requires as input the minimum number of individuals in each group needed in order to make the call (fixed/poly), and a prefix for output files. Key file is a two column tab-delimited file with sample ID and population code.

```
module load perl

input_vcf=/Mt_v1_chr1_qd_varsite.vcf.gz
key_file=/Mt_v1_SampleID_pop_card_tris.txt
prefix=chr1

perl /ParseVCF_fixed_shared.pl \
$input_vcf Tristrami_High,Cardinalis_SYM $key_file 5 $prefix
```

## private alleles
Used [custom perl script](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/ParseVCF_private.pl) to identify alleles private to parses a multi-sample gVCF file and, given a focal species with several locations (as identified by a key), plus at least one other species of interest, returns a tab-delineated text file, with CHROM\tPOS for each site where there is a private ALLELE in the focal SPECIES not present in the others (excluding putative hybrids), the identity of that allele, followed by one column for each location giving a Y/N for presence/absence of that private allele. The assumed format of the "key" file corresponds to the Myzo_WGS metadata files (full metadata), with the (currently!) important columns the 3rd (sample location), 4th (species) and 7th (sample name). Also requires as input the minimum number of individuals in each species (note: NOT location sample) needed in order to make the call as to whether an allele is private. As written, the program ignores any individual with the "species" name "Hybrid" in the key file. Name accordingly! Alternatively, if you give the program an incomplete key, it will ignore any sample in the vcf that it has not already seen in the key, so you could simply delete all hybrids (and possibly the "MO" cardinalis sample) from the key being used. Prints to stdout. 

```
module load perl

input_vcf=/Mt_v1_chr1_qd_varsite_recode.vcf.gz
meta_file=/Mt_v1_WGS_card_tris_meta_corrected.csv
prefix=private_card_chr1

perl /scratch/juy3_lab/MyzoWGS/Mt_v1_WGS/fixed_shared_private/fixed_shared_private_scripts/ParseVCF_private.pl \
$input_vcf Cardinalis $meta_file 5 > $prefix.txt
```
## PCA
Used R v4.1.1 with packages tidyverse, gdsfmt, and SNPRelate
Nice tutorial used to build my Rscript: https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html

## ADMIXTURE
Further filtering of quality/depth variant sites VCFs to impost minor allele frequency threshold of 0.05
```
module load vcftools
module load samtools/1.7
module load bcftools

vcftools --gzvcf /Mt_v1_autos_qd_varsite_recode.vcf.gz \
--maf 0.05 \
--recode --recode-INFO-all --stdout | bgzip -c > /Mt_v1_autos_qd_maf05_varsite_recode.vcf.gz

tabix /Mt_v1_autos_qd_maf05_varsite_recode.vcf.gz
```

Then remove sites in linkage disequilibrium using plink - an extra initial step with bcftools to annotate/rename variants in the VCF file
```
module load plink/1.9
module load vcftools
module load samtools/1.7
module load bcftools


bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' /Mt_v1_autos_qd_maf05_varsite_recode.vcf.gz -O z --output /Mt_v1_autos_qd_maf05_varsite_recode_ann.vcf.gz

plink --vcf /Mt_v1_autos_qd_maf05_varsite_recode_ann.vcf.gz --indep-pairwise 50 5 0.5 --allow-extra-chr --chr-set 80 --out /Mt_v1_autos_qd_maf05_varsite_recode_LD

plink --vcf /Mt_v1_autos_qd_maf05_varsite_recode_ann.vcf.gz --extract /Mt_v1_autos_qd_maf05_varsite_recode_LD.prune.in --make-bed --allow-extra-chr --chr-set 80 --out /Mt_v1_autos_qd_maf05_varsite_recode_LD

plink --bfile /Mt_v1_autos_qd_maf05_varsite_recode_LD --recode vcf --allow-extra-chr --chr-set 80 --out /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned

bgzip /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.vcf

tabix /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.vcf.gz
```
Then, use an sh script to rename autosomes to remove any characters as these will not work with admixture program.
```
#!/bin/bash

sed -i.bak -e 's/chr1A/60/g' /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.bim
sed -i.bak -e 's/chr4A/70/g' /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.bim
sed -i.bak -e 's/chr5_remnant/50/g' /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.bim
sed -i.bak -e 's/chr5_W2seg4-neoPAR_fusion/54/g' /Mt_v1_autos_qd_maf05_varsite_recode_LDpruned.bim
```
Finally, ready for ADMIXTURE analysis
```
module load admixture

for K in 1 2 3 4 5 6 7 8 9 10; \
do admixture --cv /Mt_v1_autos_qd_maf05_varsite_recode_LD.bed $K | tee /Mt_v1_autos_log${K}.out; done
```
plotted results in R v4.1.1 using packages tidyverse, ggplot2, forcats, ggthemes, and patchwork

## triangle plot
First, identify SNPs fixed between species using allopatric individuals and vcftools
```
module load vcftools
module load samtools/1.7

vcftools --gzvcf /Mt_v1_autos_UgiThreeMak_qd_maf05_bi_nomiss_varsite_recode.vcf.gz \
--weir-fst-pop /TrisAllo_SampleID.txt \
--weir-fst-pop /Mt_v1_UgiThree.txt \
--out /triangle/Mt_v1_autos_UgiThree_High_qd_maf05_bi_nomiss_varsite
```

then, use vcftools to make a VCF file with only sympatric individuals and the snps fixed between the two species in allopatry
```
module load vcftools
module load samtools/1.7
module load bcftools

vcftools --gzvcf /Mt_v1_autos_UgiThreeMak_qd_maf05_bi_nomiss_varsite_recode.vcf.gz \
--keep /Mt_v1_Sympatric.txt \
--positions /diagSNPs_Fstfixed_UgiThreeHigh_qd_maf05_bi_nomiss.txt \
--recode --recode-INFO-all --stdout | bgzip -c > /Mt_v1_autos_Sym_qd_maf05_bi_nomiss_fixedSNP_recode.vcf.gz

tabix /Mt_v1_autos_Sym_qd_maf05_bi_nomiss_fixedSNP_recode.vcf.gz
```
Next, pull VCF into R using vcfR and use introgress package to calculate interspecific heterozygosity and hybrid index. See [Mt_v1_triangle.Rmd](https://github.com/ehshogren/MyzomelaPopulationGenomics/blob/main/Mt_v1_triangleplot.Rmd) file for details.

## ABBA-BABA
Used Dsuite program to assess introgression using ABBA-BABA analysis, *D* statistics and *f4* admixture ratios
```
module load gcc/4.9.4
module load zlib

program_path=/Dsuite/Build/
input_vcf=/Mt_v1_autos_qd_varsite_recode.vcf.gz
pop_sets=/Mt_v1_Dsuite_SampleID_TrisHigh_TrisSym_CardSym_Out.tsv
 
$program_path/Dsuite Dtrios -o /Mt_v1_autos_tA_tS_cS_Dtrios $input_vcf $pop_sets 
```
looped over chromosomes to calculate statistics for each.
```
module load gcc/4.9.4
module load zlib
module load bcftools

program_path=/Dsuite/Build/
input_vcf=/Mt_v1_autos_qd_varsite_recode.vcf.gz
pop_sets=/Mt_v1_Dsuite_SampleID_TrisHigh_TrisSym_CardSym_Out.tsv
tree_path=/tA_tS_cS_tree.nwk
chrom_names=/Mt_v1.0_autos_chrom.txt
 
cat $chrom_names |
while read chrom; do
NUMLINES=$(bcftools view -r ${chrom} $input_vcf | wc -l)
bcftools view -r ${chrom} $input_vcf | $program_path/Dsuite Dtrios -l $NUMLINES -c -t $tree_path -o /Mt_v1_${chrom}_tA_tS_cS_Dtrios stdin $pop_sets;
done 
```
