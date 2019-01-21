# In order to run bammds, you need to make sure of two things:
# 1. The chromosomes in the panel are named in the same way as in the bam
# (e.g., all NN or all chrNN)
# 2. That your bam file contains only the chromosomes that are reported
# in the panel

# These two objectives can be achieved by reheading a  bam file:
ch=($(seq 1 22))
ind=$1
name=$2


if [ ! -e keep_chrs.txt ]
then
  for chr in "${ch[@]}"
  do
    echo "SN:$chr"
  done >keep_chrs.txt
fi

if [ ! -d Reheaded ]
then
 mkdir Reheaded
fi

samtools view -H $ind.bam > Reheaded/$name.header
sed -i 's/SN:chr/SN:/' Reheaded/$name.header

grep '^@SQ' Reheaded/$name.header |grep -vf keep_chrs.txt  > Reheaded/$name.remove
grep -vf Reheaded/$name.remove Reheaded/$name.header > Reheaded/$name.new.header
rm Reheaded/$name.header Reheaded/$name.remove

samtools reheader Reheaded/$name.new.header ${ind}.bam > Reheaded/$name.bam
samtools index Reheaded/$name.bam
