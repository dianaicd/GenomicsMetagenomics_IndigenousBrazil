#!/bin/bash

#SBATCH --account=amalaspi_virome
#SBATCH --partition=ax-normal
#SBATCH --error=split_%A_%a.err
#SBATCH --output=output_%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time 20:00:00
#SBATCH --array=0-23%7

module add UHTS/Analysis/picard-tools/2.18.11

#Jan 20 2020
#Split by libraries using picard tools.

#Jan 29 2020
#Array mode.

MNS=('MN00010' 'MN00013' 'MN00016' 'MN00019' 'MN00021' 'MN00022' 'MN00023' 'MN00039' 'MN0003' 'MN00045' 'MN00056' 'MN00064' 'MN00066' 'MN00067' 'MN00068' 'MN00069' 'MN0008' 'MN0009' 'MN00118' 'MN00119' 'MN00316' 'MN00346' 'MN01701' 'MN1943')

mkdir ${MNS[$SLURM_ARRAY_TASK_ID]}_libs
picard-tools SplitSamByLibrary I=${MNS[$SLURM_ARRAY_TASK_ID]}.hg19_low_qual.bam O=${MNS[$SLURM_ARRAY_TASK_ID]}_libs VALIDATION_STRINGENCY=LENIENT
