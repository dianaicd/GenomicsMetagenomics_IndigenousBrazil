wd=~/scratch_monthly/Botocudos/BAM/2018_09_11/
cd $wd
mn=($(cat mn_ids.txt))

for ind in ${mn[@]}
do
  echo $ind
  cd $ind
  ls */*calmd.bam > list.bam
  first=$(head -n1 list.bam)
  samtools view -H $first |grep -v "@RG" |grep -v "@PG" >header.txt
  for sam in $(cat list.bam)
  do
    samtools view -H $sam | grep "@RG"  >> header.txt
    samtools view -H $sam | grep "@PG"  >> header.txt
  done

  samtools merge -h header.txt -b list.bam -f test.bam --threads 10
  picard-tools MarkDuplicates I=test.bam O=${mn}.bam
  samtools index ${mn}.bam
  rm test.bam
  cd $wd
done
