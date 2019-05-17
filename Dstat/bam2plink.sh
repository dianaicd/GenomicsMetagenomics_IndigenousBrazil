
# Script to:
#  - Sample a random base per defined site (output: plink)
#  - Merge multiple plinks to plink panel
#  - Sample a random allele for all the extended panel

panel=~/archive/Panels/fromVictor/Maanasa_0.01_minind
bamlist=bamlist

pathPanel=$(dirname $panel)
panel=$(basename $panel .bed)


ln -s ~/data/Git/Botocudos-scripts/MDS/sample_ped.pl ./
ln -s $pathPanel/$panel.{bed,bim,fam} ./

if [ ! -e ${panel}_regions ]
then
    cat $panel.bim |awk '{print $1"\t"$4-1"\t"$4"\t"$5"\t"$6}' \
    >${panel}_regions
fi

#-----------------------------------------------------------------------------#
# Sample a random allele from the bam files
while read line
    do
    echo $line
    ind=$(basename $line .bam)
    bam2plink.py bamfile=$line plinkpref=$panel \
    MinBQ=20 indname=$ind.sampled popname=$ind doCns=F \
    trim=5 MinMQ=30 >out_${ind}_bam2plink.txt \
    2> err_${ind}_bam2plink.txt &
done < $bamlist

{ sleep 5; echo waking up after 5 seconds; } &
  wait
  echo all jobs are done!

#-----------------------------------------------------------------------------#
# Merge plink binary files to make an extended panel
extendedPanel=$panel

while read line 
do
    ind=$(basename $line .bam)
    echo "-------------------- Merging $ind to $panel ---------------------"

    plink --bfile $extendedPanel --bmerge ${ind}.sampled_${ind}_${panel} \
      --make-bed --out ${panel}_${ind} --allow-no-sex
    rm $extendedPanel.*
    extendedPanel=${panel}_${ind}

done < $bamlist

for ext in bed bim fam log nosex
do
  mv $extendedPanel.$ext ${panel}_${bamlist}.$ext
done

#-----------------------------------------------------------------------------#
# Sample a random allele for the whole panel; this is useful to do an MDS

plink --recode --bfile ${panel}_${bamlist} --out ${panel}_${bamlist}

perl sample_ped.pl -ped ${panel}_${bamlist}.ped \
  -out ${panel}_${bamlist}_sampled.ped

cp ${panel}_${bamlist}.map ${panel}_${bamlist}_sampled.map

plink --distance square gz 'flat-missing' \
 --file ${panel}_${bamlist}_sampled \
 --out ${panel}_${bamlist}_sampled
