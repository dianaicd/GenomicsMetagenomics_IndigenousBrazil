#!/bin/bash
#SBATCH --account amalaspi_virome
#SBATCH --job-name unmapped
#SBATCH --partition ax-normal
#SBATCH --time 23:30:15
#SBATCH --mem-per-cpu 10G
#SBATCH --nodes 1
#SBATCH --cpus-per-task 4
#SBATCH --error unmapped_bams/unmapped_job%j.txt
#SBATCH --output unmapped_bams/unmapped_job%j.txt

module add UHTS/Analysis/samtools/1.8

#Check flags
for f in MN00010 MN00013 MN00016 MN00019 MN00021 MN00022 MN00023 MN00039 MN0003 MN00045 MN00056 MN00064 MN00066 MN00067 MN00068 MN00069 MN0008 MN0009 MN00118 MN00119 MN00316 MN00346; do
	samtools flagstat ${f}_low_qual.bam > unmapped_bams/${f}_flagstat
done

#wait

#Get unmapped reads
for f in MN00010 MN00013 MN00016 MN00019 MN00021 MN00022 MN00023 MN00039 MN0003 MN00045 MN00056 MN00064 MN00066 MN00067 MN00068 MN00069 MN0008 MN0009 MN00118 MN00119 MN00316 MN00346; do
	samtools view -b -f 4 ${f}.hg19_low_qual.bam > unmapped_bams/${f}_unmapped.bam
done

#wait

#Check flags for unmapped reads
for f in unmapped_bams/*_unmapped.bam; do
	samtools flagstat ${f} > ${f%%_unmapped.bam}_unmapped_flagstat
done

#wait

#Get FASTA files
for f in unmapped_bams/*_unmapped.bam; do
	f_name=${f#unmapped_bams/}
	f_name=unmapped_fnas/$(echo ${f_name%%_unmapped.bam})_unmapped.fna
	samtools fasta ${f} > ${f_name}
done

#wait
