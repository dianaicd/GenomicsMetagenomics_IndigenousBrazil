#!/bin/bash
#SBATCH --account amalaspi_virome
#SBATCH --job-name unmapped
#SBATCH --partition ax-normal
#SBATCH --time 10:00:00
#SBATCH --mem-per-cpu 7G
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --error unmapped_bams/unmapped_%A_%j.txt
#SBATCH --output unmapped_bams/unmapped_%A_%j.txt
#SBATCH --array=0-23

module add UHTS/Analysis/samtools/1.8

#Jan 29 2020.
#Script to check flags, select unmapped reads, and create fasta files.
#Fix dir name. Array mode.

MNS=('MN00010' 'MN00013' 'MN00016' 'MN00019' 'MN00021' 'MN00022' 'MN00023' 'MN00039' 'MN0003' 'MN00045' 'MN00056' 'MN00064' 'MN00066' 'MN00067' 'MN00068' 'MN00069' 'MN0008' 'MN0009' 'MN00118' 'MN00119' 'MN00316' 'MN00346')

#Check flags
samtools flagstat ${MNS[$SLURM_ARRAY_TASK_ID]}.hg19_low_qual.bam > unmapped_bams/${MNS[$SLURM_ARRAY_TASK_ID]}_flagstat

#Get unmapped reads
#samtools view -b -f 4 ${MNS[$SLURM_ARRAY_TASK_ID]}.hg19_low_qual.bam > unmapped_bams/${MNS[$SLURM_ARRAY_TASK_ID]}_unmapped.bam

#Check flags for unmapped reads
#samtools flagstat unmapped_bams/${MNS[$SLURM_ARRAY_TASK_ID]}_unmapped.bam > ${MNS[$SLURM_ARRAY_TASK_ID]}_unmapped_flagstat

#Get FASTA files
#samtools fasta unmapped_bams/${MNS[$SLURM_ARRAY_TASK_ID]}_unmapped.bam > unmapped_fnas/${MNS[$SLURM_ARRAY_TASK_ID]}_unmapped.fna
