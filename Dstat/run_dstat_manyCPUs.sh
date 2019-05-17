# 160 CPU are available
maxCPU=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
nCPU=$(echo "$maxCPU * 0.8 " |bc)
samplebam=$(ls bams/*bam |sed 's/@//'|head -n1)
totalbases=$(samtools view -H $samplebam |grep "^@SQ" |\
  sed 's/LN://' |awk 'BEGIN {FS = "\t" } ; {sum+=$3} END {print sum}')

samtools view -H $samplebam |grep "^@SQ" |\
  cut -f2,3 |sed 's/SN://;s/LN://' > chrom_size.txt
blockSize=$(echo "$totalbases / $nCPU" | bc)
Rscript print_regions.R $blockSize

for region in $(ls region_*.txt)
do
    sites=$(echo $region |sed 's/region_//; s/.txt//')

    angsd -doAbbababa2 1 -bam ordered_bam.filelist \
        -sizeFile sizeFile.size -doCounts 1 \
        -out dstat_Feb25_rmtrans_${sites} -anc Yoruba.fa \
        -rf $region -useLast 1 -minQ 20 -rmTrans 1 &
done
