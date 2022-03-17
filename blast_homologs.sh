#!/bin/bash

## ABOUT THIS SCRIPT ##
# Author: Ben Long
# Date : 03.17.2022
# Description: This script blasts three Setaria italica gene sequences against the maize, sorghum, and foxtail millet genomes. "So then I started blasting" -Frank
# Run Information: This script is run manually.

## SLURM PARAMETERS ##
#SBATCH --job-name=blast_homologs              # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=8		                            # Number of cores per task
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=24:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/mplab/ben/output/blastlog.%j			    # Standard output and error log
#SBATCH --mail-user=bjl31194@uga.edu                    # Where to send mail
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
SEQDIR="/work/mplab/ben/genomes"
OUTDIR="/work/mplab/ben/output"

#download genomes

curl http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz > ${SEQDIR}/Zm-B73-REFERENCE-NAM-5.0.fa.gz

curl http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz > ${SEQDIR}/Sorghum_bicolor_NCBIv3.fa.gz

curl http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/setaria_italica/dna/Setaria_italica.Setaria_italica_v2.0.dna.toplevel.fa.gz > ${SEQDIR}/Setaria_italica_v2.0.fa.gz

#modules
ml BLAST+/2.11.0-gompi-2019b
module list
​
#make blast databases
if [ ! -f ${SEQDIR}/Zm_nt_blastdb.ndb ]; then
	makeblastdb -dbtype nucl -in ${SEQDIR}/Zm-B73-REFERENCE-NAM-5.0.fa.gz -out ${SEQDIR}/Zm_nt_blastdb
fi

if [ ! -f ${SEQDIR}/Sb_nt_blastdb.ndb ]; then
	makeblastdb -dbtype nucl -in ${SEQDIR}/Sorghum_bicolor_NCBIv3.fa.gz -out ${SEQDIR}/Sb_nt_blastdb
fi
​
if [ ! -f ${SEQDIR}/Si_nt_blastdb.ndb ]; then
	makeblastdb -dbtype nucl -in ${SEQDIR}/Setaria_italica_v2.0.fa.gz -out ${SEQDIR}/Si_nt_blastdb
fi

#blast sequences to maize db
blastn -num_threads 8 -query ${SEQDIR}/Seita.2G371500.fa -db ${SEQDIR}/Zm_nt_blastdb -out ${OUTDIR}/Zm_Seita.2G371500_blastresults.tsv -outfmt 6 -max_target_seqs 10

blastn -num_threads 8 -query ${SEQDIR}/Seita.3G253100.fa -db ${SEQDIR}/Zm_nt_blastdb -out ${OUTDIR}/Zm_Seita.3G253100_blastresults.tsv -outfmt 6 -max_target_seqs 10

blastn -num_threads 8 -query ${SEQDIR}/Seita.3G359900.fa -db ${SEQDIR}/Zm_nt_blastdb -out ${OUTDIR}/Zm_Seita.3G359900_blastresults.tsv -outfmt 6 -max_target_seqs 10


# blast to sorghum db
blastn -num_threads 8 -query ${SEQDIR}/Seita.2G371500.fa -db ${SEQDIR}/Sb_nt_blastdb -out ${OUTDIR}/Sb_Seita.2G371500_blastresults.tsv -outfmt 6 -max_target_seqs 10

blastn -num_threads 8 -query ${SEQDIR}/Seita.3G253100.fa -db ${SEQDIR}/Sb_nt_blastdb -out ${OUTDIR}/Sb_Seita.3G253100_blastresults.tsv -outfmt 6 -max_target_seqs 10

blastn -num_threads 8 -query ${SEQDIR}/Seita.3G359900.fa -db ${SEQDIR}/Sb_nt_blastdb -out ${OUTDIR}/Sb_Seita.3G359900_blastresults.tsv -outfmt 6 -max_target_seqs 10


#blast to setaria db
blastn -num_threads 8 -query ${SEQDIR}/Seita.2G371500.fa -db ${SEQDIR}/Si_nt_blastdb -out ${OUTDIR}/Si_Seita.2G371500_blastresults.tsv -outfmt 6 -max_target_seqs 10

blastn -num_threads 8 -query ${SEQDIR}/Seita.3G253100.fa -db ${SEQDIR}/Si_nt_blastdb -out ${OUTDIR}/Si_Seita.3G253100_blastresults.tsv -outfmt 6 -max_target_seqs 10

blastn -num_threads 8 -query ${SEQDIR}/Seita.3G359900.fa -db ${SEQDIR}/Si_nt_blastdb -out ${OUTDIR}/Si_Seita.3G359900_blastresults.tsv -outfmt 6 -max_target_seqs 10
