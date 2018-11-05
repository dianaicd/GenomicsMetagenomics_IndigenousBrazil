#!/bin/bash

date=2018_07_19
wd=~/scratch_monthly/Botocudos/BAM/$date

ind=($(cat ma_ids.txt))
f=($(cat fastq.txt))
#ref=~/archive/References/Human/hg19/chromosomes_1-22_XY_mito.fasta
ref=/archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa

cd $wd

for i in $(seq 0 22)
do

  mkdir ${ind[$i]}
  cd ${ind[$i]}

  echo ${ind[$i]}
  mapping_aDNA.sh --fastq1 ${f[$i]} --base ${ind[$i]}_Human \
    --ref $ref \
    --noClean \
    -p 9 --unmapped > out_Human_${ind[$i]} 2> err_Human_${ind[$i]} \ &
  
  cd $wd

done
