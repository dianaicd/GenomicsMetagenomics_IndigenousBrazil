# In order to run bammds, you need to make sure of two things:
# 1. The chromosomes in the panel are named in the same way as in the bam
# (e.g., all NN or all chrNN)
# 2. That your bam file contains only the chromosomes that are reported
# in the panel

# These two objectives can be achieved by reheading a  bam file:
ch=($(seq 1 22) MT X Y)
ind=$1
name=$2

samtools view -H $ind.bam >$name.header

if [ ! -e keep_chrs.txt ]
then
  for chr in "${ch[@]}"
  do
    echo "SN:$chr"
  done >keep_chrs.txt
fi

grep '^@SQ' $name.header |grep -vf keep_chrs.txt  >$name.remove
grep -vf $name.remove $name.header > $name.new.header
rm $name.header $name.remove
samtools reheader $name.new.header $ind.bam >$name.reheaded.bam
samtools index $name.reheaded.bam
