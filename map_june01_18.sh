#!/bin/bash

ind=($(cat ma_ids.txt))
f=($(cat fastq.txt))
ref=~/archive/References/HBV/Hepatitis_B_ref.fasta
dataset1=~/archive/References/HBV/HBV_ref_10.fasta
a=~/archive/References/HBV/AM282986.fasta
c=~/archive/References/HBV/AB117758.fasta
suf=~/archive/References/HBV/
id=(NC_003977 D23678 AM282986 AF241409 AB126581 AY090454 AB486012 AB117758 AB205192 AB056513 X69798)
genotype=(ReferenceHBV B A I D H J C E G F)

for i in $(seq 0 22)
do

#  mkdir ${ind[$i]}
  cd ${ind[$i]}
  echo ${ind[$i]}

  for g in 1 2 3 4 5 6 7 8 9 10
  do

    mapping_aDNA.sh --fastq1 ${f[$i]} --base ${ind[$i]}_${genotype[$g]}_Q0 \
      --ref ${suf}${id[$g]}.fasta --skip1 \
      -p 9 -q 0 > out_${genotype[$g]}_${ind[$i]}_Q0 2> err_${genotype[$g]}_${ind[$i]}_Q0 &
  done  
  
  cd ~/scratch_monthly/HBV/


done
