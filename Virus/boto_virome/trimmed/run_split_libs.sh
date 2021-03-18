#!/bin/bash

#SBATCH --account=amalaspi_virome
#SBATCH --partition=ax-normal
#SBATCH --error=%j.log
#SBATCH --output=output_%j_smallrun.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time 23:00:00

module add UHTS/Analysis/picard-tools/2.18.11

#Jan 20 2020
#Split by libraries using picard tools.

for i in *_low_qual.bam; do
	picard-tools SplitSamByLibrary I=${i} O=${i%%.hg19_low_qual.bam}_libs VALIDATION_STRINGENCY=LENIENT
done
