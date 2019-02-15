panel=Maanasa_americas_reheaded_filtered
bamlist=24indLS.bam.list

CHR=($(seq 1 22))

for chr in ${CHR[@]}
do
    echo $chr
    angsd -out out_${chr} -doCounts 1 -dumpCounts 4 -bam $bamlist \
    -rf ${panel}_${chr}_sites.txt -minQ 20 -nThreads 8 -p 2 &
done
