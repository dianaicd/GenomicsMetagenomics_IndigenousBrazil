#!/bin/bash
#SBATCH --account amalaspi_virome
#SBATCH --job-name unmapped
#SBATCH --partition ax-normal
#SBATCH --time 4:00:15
#SBATCH --mem-per-cpu 6G
#SBATCH --nodes 1
#SBATCH --cpus-per-task 4
#SBATCH --error unmapped_bams/unmapped_job%j.txt
#SBATCH --output unmapped_bams/unmapped_job%j.txt

module add UHTS/Analysis/samtools/1.8
#Jan 23 2020.
#Chack flags, select unmapped reads from samples MN01701 and MN1943

#Check flags
samtools flagstat MN01701_libs/MN01701_merged_L1_L2.bam > unmapped_bams/MN01701_merged_L1_L2_flagstat
#samtools flagstat MN1943_libs/MN1943_L1.bam > unmapped_bams/MN1943_L1_flagstat

#Get unmapped reads
samtools view -b -f 4 MN01701_libs/MN01701_merged_L1_L2.bam > unmapped_bams/MN01701_unmapped.bam
#samtools view -b -f 4 MN1943_libs/MN1943_L1.bam > unmapped_bams/MN1943_unmapped.bam

#Check flags for unmapped reads
for f in unmapped_bams/MN01701_unmapped.bam; do #unmapped_bams/MN1943_unmapped.bam; do
	samtools flagstat ${f} > ${f%%_unmapped.bam}_unmapped_flagstat
done

#Get FASTA files
for f in unmapped_bams/MN01701_unmapped.bam; do #unmapped_bams/MN1943_unmapped.bam; do
	f_name=${f#unmapped_bams/}
	f_name=unmapped_fnas/$(echo ${f_name%%_unmapped.bam})_unmapped.fna
	samtools fasta ${f} > ${f_name}
done

