#!/bin/bash

#SBATCH --account=amalaspi_virome
#SBATCH --partition=ax-normal
#SBATCH --error=%j.log
#SBATCH --output=output_%j_smallrun.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time 01:00:00

module add UHTS/Analysis/samtools/1.8

#Jan 17 2020
#Check which libraries are present in bam files.

for i in *_low_qual.bam; do
	samtools view -H ${i} > ${i%%_low_qual.bam}_header.txt
	grep "^@RG" ${i%%_low_qual.bam}_header.txt > ${i%%_low_qual.bam}_rg.txt
	cut -f3 ${i%%_low_qual.bam}_rg.txt | sort | uniq > ${i%%_low_qual.bam}_libs.txt
done
