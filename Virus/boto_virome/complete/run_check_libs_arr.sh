#!/bin/bash

#SBATCH --account=amalaspi_virome
#SBATCH --partition=ax-normal
#SBATCH --error=libcheck_job_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time 00:45:00
#SBATCH --array=0-23

module add UHTS/Analysis/samtools/1.8

#Jan 17 2020
#Check which libraries are present in bam files.

#Jan 29 2020
#Array mode.

MNS=('MN00010' 'MN00013' 'MN00016' 'MN00019' 'MN00021' 'MN00022' 'MN00023' 'MN00039' 'MN0003' 'MN00045' 'MN00056' 'MN00064' 'MN00066' 'MN00067' 'MN00068' 'MN00069' 'MN0008' 'MN0009' 'MN00118' 'MN00119' 'MN00316' 'MN00346' 'MN01701' 'MN1943')

samtools view -H ${MNS[$SLURM_ARRAY_TASK_ID]}.hg19_low_qual.bam > ${MNS[$SLURM_ARRAY_TASK_ID]}_header.txt
grep "^@RG" ${MNS[$SLURM_ARRAY_TASK_ID]}_header.txt > ${MNS[$SLURM_ARRAY_TASK_ID]}_rg.txt
cut -f3 ${MNS[$SLURM_ARRAY_TASK_ID]}_rg.txt | sort | uniq > ${MNS[$SLURM_ARRAY_TASK_ID]}_libs.txt

