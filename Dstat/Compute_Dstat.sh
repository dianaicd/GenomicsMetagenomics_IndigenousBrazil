#
# Script to compute the D-statistic in ANGSD
#

# Chromosomes on which we will compute the D-stat
CHR=($(seq 1 22))
# Basename for the output. E.g., "24Botocudos"
ind=MN0008
# Basename of the file with the paths to the panel and
# target's bam files
panel=fromVictor
# Fasta file of outgroup's genome
anc=Yoruba.fa

# Run angsd, per chromosome, removing transitions
if [ ! -d Raw ]
then
  mkdir Raw
fi

for chr in ${CHR[@]}
do 
  echo $chr
  for bam in BAM/*bam
  do
    if [ ! -d BAM/$chr ]
      then
      mkdir BAM/$chr
    fi
    name=$(basename $bam .bam )
    samtools view -b $bam $chr >BAM/$chr/${name}_${chr}.bam &
  done 
  
  angsd -doabbababa 1 \
  -out Raw/${ind}_${panel}_${chr}_rmtrans \
  -bam ${panel}.txt \
  -doCounts 1 -useLast 1 \
  -minQ 20 -minMapQ 30 -r ${chr} -nThreads 4 &
done
{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
wait
echo all jobs are done!

# Run angsd, per chromosome, removing transitions, 
# and correcting for error
if [ ! -d Corrected ]
then
  mkdir Corrected
fi

nThreads=8
for chr in ${CHR[@]}
do 
  echo $chr
  angsd -doAbbababa2 1 \
  -out Corrected/${ind}_${panel}_${chr}_rmtrans \
  -bam ${ind}_${panel}.txt \
  -doCounts 1 -anc $anc \
  -minQ 20 -minMapQ 30 -r ${chr} -nThreads $nThreads \
  -rmTrans 1 -useLast 0 \
  -sizeFile ${ind}_${panel}.size >out_${ind}_${panel}_${chr}_rmtrans \
  2>err_${ind}_${panel}_${chr}_rmtrans &
done

{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
  wait
  echo all jobs are done!
