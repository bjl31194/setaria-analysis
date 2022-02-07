#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 12.6.2021
# Description: This script runs QC on pooled Setaria reads from 3 phenotyped groups, maps it to the Yugu1 reference genome, calls variants (SNPS/indels), and jointly genotypes samples using GATK :)
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=setaria_call_variants_GATK	              # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task
#SBATCH --mem=32gb			                                # Total memory for job
#SBATCH --time=3-00:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/bjl31194/setaria/log.%x_%j			    # Standard output and error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

##set input and output directory variables
OUTDIR="/work/gene8940/bjl31194/setaria"
DATADIR="/work/gene8940/instructor_data"

##if output directory doesn't exist, create it
if [ ! -d ${OUTDIR} ]
then
    mkdir -p ${OUTDIR}
fi

##load modules
ml FastQC/0.11.9-Java-11 BWA/0.7.17-GCC-8.3.0 SAMtools/1.10-GCC-8.3.0 BCFtools/1.10.2-GCC-8.3.0 GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

##download data: this is the NCBI RefSeq FTP link, but for this analysis I will be using the Phytozome version 2.2 downloaded locally to OUTDIR
#curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/263/155/GCF_000263155.2_Setaria_italica_v2.0/GCF_000263155.2_Setaria_italica_v2.0_genomic.fna.gz | gunzip > ${OUTDIR}/GCF_000263155.2.fna

##make indexed fasta of ref with only first 9 scaffolds (chromosomes)
samtools faidx ${OUTDIR}/Sitalica_312_v2.fa
samtools faidx ${OUTDIR}/Sitalica_312_v2.fa scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8 scaffold_9 > ${OUTDIR}/Sitalica_312_v2_1-9.fa
samtools faidx ${OUTDIR}/Sitalica_312_v2_1-9.fa
samtools dict ${OUTDIR}/Sitalica_312_v2_1-9.fa

##run qc on reads
fastqc ${OUTDIR}/SI*.fastq.gz

##index reference
bwa index ${OUTDIR}/Sitalica_312_v2_1-9.fa

##map to ref, convert to BAM, and sort, + create .bai index
bwa mem -t 8 -R "@RG\tID:ALB\tSM:SI_1716_A\tPL:illumina" ${OUTDIR}/Sitalica_312_v2_1-9.fa ${OUTDIR}/SI_1716_A_R1.fastq.gz ${OUTDIR}/SI_1716_A_R2.fastq.gz | samtools view - -O BAM -o ${OUTDIR}/SI_1716_A_1-9.bam

bwa mem -t 8 -R "@RG\tID:GRN\tSM:SI_1716_G\tPL:illumina" ${OUTDIR}/Sitalica_312_v2_1-9.fa ${OUTDIR}/SI_1716_G_R1.fastq.gz ${OUTDIR}/SI_1716_G_R2.fastq.gz | samtools view - -O BAM -o ${OUTDIR}/SI_1716_G_1-9.bam

bwa mem -t 8 -R "@RG\tID:YEL\tSM:SI_1716_Y\tPL:illumina" ${OUTDIR}/Sitalica_312_v2_1-9.fa ${OUTDIR}/SI_1716_Y_R1.fastq.gz ${OUTDIR}/SI_1716_Y_R2.fastq.gz | samtools view - -O BAM -o ${OUTDIR}/SI_1716_Y_1-9.bam

bwa mem -t 8 -R "@RG\tID:575\tSM:SI_575\tPL:illumina" ${OUTDIR}/Sitalica_312_v2_1-9.fa ${OUTDIR}/SI_575_R1.fastq.gz ${OUTDIR}/SI_575_R2.fastq.gz | samtools view - -O BAM -o ${OUTDIR}/SI_575_1-9.bam

samtools sort --threads 8 ${OUTDIR}/SI_1716_A_1-9.bam -o ${OUTDIR}/SI_1716_A_1-9_sorted.bam
samtools sort --threads 8 ${OUTDIR}/SI_1716_G_1-9.bam -o ${OUTDIR}/SI_1716_G_1-9_sorted.bam
samtools sort --threads 8 ${OUTDIR}/SI_1716_Y_1-9.bam -o ${OUTDIR}/SI_1716_Y_1-9_sorted.bam
samtools sort --threads 8 ${OUTDIR}/SI_575_1-9.bam -o ${OUTDIR}/SI_575_1-9_sorted.bam

samtools index ${OUTDIR}/SI_1716_A_1-9_sorted.bam
samtools index ${OUTDIR}/SI_1716_G_1-9_sorted.bam
samtools index ${OUTDIR}/SI_1716_Y_1-9_sorted.bam
samtools index ${OUTDIR}/SI_575_1-9_sorted.bam

#call haplotypes using *gasp* GATK's Haplotype Caller module
gatk HaplotypeCaller -R ${OUTDIR}/Sitalica_312_v2_1-9.fa -I ${OUTDIR}/SI_1716_A_1-9_sorted.bam -ERC GVCF -O ${OUTDIR}/SI_1716_A_1-9_gatk.vcf.gz

gatk HaplotypeCaller -R ${OUTDIR}/Sitalica_312_v2_1-9.fa -I ${OUTDIR}/SI_1716_G_1-9_sorted.bam -ERC GVCF -O ${OUTDIR}/SI_1716_G_1-9_gatk.vcf.gz

gatk HaplotypeCaller -R ${OUTDIR}/Sitalica_312_v2_1-9.fa -I ${OUTDIR}/SI_1716_Y_1-9_sorted.bam -ERC GVCF -O ${OUTDIR}/SI_1716_Y_1-9_gatk.vcf.gz

gatk HaplotypeCaller -R ${OUTDIR}/Sitalica_312_v2_1-9.fa -I ${OUTDIR}/SI_575_1-9_sorted.bam -ERC GVCF -O ${OUTDIR}/SI_575_1-9_gatk.vcf.gz

##merge GVCF files

gatk CombineGVCFs -R ${OUTDIR}/Sitalica_312_v2_1-9.fa --variant ${OUTDIR}/SI_1716_A_filtered_gatk.vcf.gz --variant ${OUTDIR}/SI_1716_G_filtered_gatk.vcf.gz --variant ${OUTDIR}/SI_1716_Y_filtered_gatk.vcf.gz --variant ${OUTDIR}/SI_575_filtered_gatk.vcf.gz -O ${OUTDIR}/SI_filtered_merged_gatk.vcf.gz

#gatk GenomicsDBImport -genomicsdb-workspace-path ${OUTDIR}/gatk_genomicsdb --variant ${OUTDIR}/SI_1716_A_1-9_gatk.vcf.gz --variant ${OUTDIR}/SI_1716_G_1-9_gatk.vcf.gz --variant ${OUTDIR}/SI_1716_Y_1-9_gatk.vcf.gz --variant ${OUTDIR}/SI_575_1-9_gatk.vcf.gz

##create merged file where all pools have been genotyped

gatk GenotypeGVCFs -R ${OUTDIR}/Sitalica_312_v2_1-9.fa --variant ${OUTDIR}/SI_filtered_merged_gatk.vcf.gz -O ${OUTDIR}/SI_filtered_snps.vcf

gatk VariantsToTable -V SI_raw_SNPs_gatk_haplo.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -GF GT -GF AD -O ${OUTDIR}/test_table.tsv
