

# Script to count reads (works on flagstat idxstats and bam)
bam=$1
chrs=( "$2" )
seq_retained_reads=$3
type=$4

hits=""
#"hits_raw"
hits_raw=0
if [ ! -e ${bam}_idxstats.txt ]
then
  samtools idxstats ${bam}.bam > ${bam}_idxstats.txt
fi

for chr in ${chrs[@]}
do
  x=$(grep -P $(echo "^$chr\t") ${bam}_idxstats.txt |cut -f 3)
  #hits_raw=$(samtools view  -c ${bam}.bam ${chrs[@]})
  hits_raw=$(echo $hits_raw + $x |bc -l )
done
#echo "hits_raw = $hits_raw"
hits=$(echo $hits$'\t'$hits_raw)

#[22] "hits_raw_frac"
hits_raw_frac=$(echo "$hits_raw / $seq_retained_reads" |bc -l)
#echo "hits_raw_frac = $hits_raw_frac"
hits=$(echo $hits$'\t'$hits_raw_frac)

#"hits_clonality"
hits_clonality=$(samtools view -f1024 -c ${bam}.bam ${chrs[@]})
#echo "hits_clonality = $hits_clonality"
hits=$(echo $hits$'\t'$hits_clonality)

#"hits_unique"
hits_unique=$(samtools view -cF1024 ${bam}.bam ${chrs[@]})
#hits_unique=$(echo $hits_raw - $hits_clonality |bc -l )
#echo "hits_unique = $hits_unique"
hits=$(echo $hits$'\t'${hits_unique})

#[25] "hits_unique_frac"
hits_unique_frac=$(echo "$hits_unique / $seq_retained_reads" | bc -l)
#echo "hits_unique_frac = $hits_unique_frac"
hits=$(echo $hits$'\t'${hits_unique_frac})

#"hits_coverage"
hits_coverage=0
bases=0
length=0
for chr in ${chrs[@]}
do
  if  grep -q "^$chr:" ${bam}_depth.txt
  then
    x=$(grep "^$chr:" ${bam}_depth.txt |cut -f2)
  else
    x=0
  fi

  bases=$(echo $bases + $x |bc -l)
  x=$(samtools view -H ${bam}.bam |grep -P $(echo "SN:${chr}\t") |cut -f 3 |sed 's/LN://')
  length=$(echo $length + $x |bc -l )

done

hits_coverage=$( echo "$bases / $length " |bc -l )
#echo "hits_coverage = $hits_coverage"
hits=$(echo $hits$'\t'${hits_coverage})

#"hits_length.mitochondrial."
hits_length=$(samtools view ${bam}.bam ${chrs[@]} |perl ~/data/Scripts/length.pl -o ${bam}_${type}_length.txt -type $type)
#echo "hits_length = $hits_length"
hits=$(echo $hits$'\t'${hits_length})

echo ${hits}
