#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J mapHuman[1-158]
#BSUB -o log/out.%I.txt
#BSUB -e log/err.%I.txt
#BSUB -R "rusage[mem=10000]"
#BSUB -M 10000000

#### Your shell commands below this line ####
date=2018_09_11
wd=~/scratch_monthly/Botocudos/BAM
cd $wd

# Reference sequence
ref=/archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa
echo $ref
if [ ! -e ma_mn.txt ]
then
  paste ma_ids.txt mn_ids.txt >ma_mn.txt
fi

#path to fastq file
fastq=(0 $(cat fastq.txt|  sed 's/,/ /g' fastq.txt ))

f=${fastq[$LSB_JOBINDEX]}
echo $f

# MA id
ma=$(echo $f | rev |cut -f 1 -d '/' |rev | cut -f 3 -d '_' )
echo $ma
# MN id
mn=$(grep $ma ma_mn.txt |cut -f 2)
echo $mn
# Date in which I downloaded the data
seq_date=$( echo $f | cut -f 7 -d '/' )
echo $seq_date
# Basename for bam file
base_bam=$( echo $f | rev | cut -f 1 -d '/' | rev |sed 's/.fastq.gz//' )
echo $base_bam
# Directory structure:
# date/MNxxxxx/seq_date/
if [ ! -d $date ]
then
  mkdir $date
fi
cd $date

if [ ! -d $mn ]
then
  mkdir $mn
fi
cd $mn

if [ ! -d $seq_date ]
then
  mkdir $seq_date
fi
cd $seq_date

mapping_aDNA.sh --fastq1 ${f} --base $base_bam \
  --ref $ref --noClean -p 9 -q 30 --save_unmapped --save_low_qual
