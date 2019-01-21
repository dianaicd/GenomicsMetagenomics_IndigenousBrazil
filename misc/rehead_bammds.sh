# In order to run bammds, you need to make sure of two things:
# 1. The chromosomes in the panel are named in the same way as in the bam
# (e.g., all NN or all chrNN)
# 2. That your bam file contains only the chromosomes that are reported
# in the panel

# These two objectives can be achieved by reheading a  bam file:
ch=($(seq 1 22))
ind=$1
name=$2

if [ ! -d Reheaded ]
then
 mkdir Reheaded
fi

samtools view -H $ind.bam > Reheaded/$name.header

if [ ! -e Reheaded/$name.newheader.txt ]
then
  head -n1 Reheaded/$name.header > Reheaded/$name.newheader.txt
  for chr in "${ch[@]}"
  do
    grep -P "SN:chr$chr\t" Reheaded/$name.header |sed 's/chr//' >> Reheaded/$name.newheader.txt 
  done
  grep '@RG' Reheaded/$name.header >> Reheaded/$name.newheader.txt 
  grep '@PG' Reheaded/$name.header >> Reheaded/$name.newheader.txt
fi

rm Reheaded/$name.header

samtools reheader Reheaded/$name.newheader.txt ${ind}.bam > Reheaded/$name.bam
#samtools index Reheaded/$name.bam
