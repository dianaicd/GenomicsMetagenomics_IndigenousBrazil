for ind in $(ls *bam) ;
do
    bam=$(basename $ind .bam@)
    echo "Working on $bam"
    echo "------------------------------ First... ------------------------"
    ~/install/angsd/angsd -i $bam.bam -r X: -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 20 -out $bam
    echo "------------------------------ Second... ------------------------"
    ~/install/angsd/misc/contamination -b 5000000 -c 154900000 -k 1 -m 0.05 -d 3 -e 20 -h HapMapCEU.gz -a ${bam}.icnts.gz > ${bam}_counts
    echo "------------------------------ Third... ------------------------"
    Rscript ContaEst.R counts=${bam}_counts freqs=HapMapCEU.gz maxsites=1000 nthr=60 outfile=${bam}_res
done