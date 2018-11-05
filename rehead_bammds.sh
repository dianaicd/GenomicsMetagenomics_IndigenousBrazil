# In order to run bammds, you need to make sure of two things:
# 1. The chromosomes in the panel are named in the same way as in the bam
# (e.g., all NN or all chrNN)
# 2. That your bam file contains only the chromosomes that are reported
# in the panel

# These two objectives can be achieved by reheading a  bam file:
ch=($(seq 1 22))
ind=$1

samtools view -H $ind.bam >$ind.header

if [ ! -e keep_chrs.txt ]
then
  for chr in "${ch[@]}"
  do
    echo "SN:$chr"
  done >keep_chrs.txt
fi

grep '^@SQ' $ind.header |grep -vf keep_chrs.txt  >$ind.remove
grep -vf $ind.remove $ind.header > $ind.new.header
rm $ind.header $ind.remove
samtools reheader $ind.new.header $ind.bam >$ind.reheaded.bam
samtools index $ind.reheaded.bam
