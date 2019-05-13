#
#
# First argument: panel's Directory
#     e.g. dir=/home/ubelix/iee/dcruz/data/Merged/500k
# Second argument: panel
#     e.g. panel=merged_noduplicates
# Third argument: K
#     e.g. k=2
#dir=/home/ubelix/iee/dcruz/data/Merged/500k
#panel=America_Oceania

SHORTOPTS="p:h:b:r"
LONGOPTS="panel: bamlist: homozygous: rmdamage:"
ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
retVal=$?
if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi
eval set -- "$ARGS"
while  [ $# -gt 0 ]; do
  case "$1" in
    -p|--panel)          panel=$2; shift;;
    -b|--bamlist)        bamlist=$2; shift;;
    -h|--homozygous)     homo=$2; shift;;
    -r|--rmdamage)       rmdamage=$2; shift;;
  esac 
  shift
done

# Directory where we can find the VCF panel
 dir=$(dirname $panel)
 type=$(echo $panel |rev |cut -f1 -d. |rev)
panel=$(basename $panel .$type)
# pass it as
# MAke only homozygous sites? (usually no)
# Path to bam files

# remove damaged sites from the panel? (usually yes)

echo "Parameters:
  dir=$dir
  panel=$panel
  name=$bamlist
  homo=$homo
  rmdamage=$rmdamage
"

#-------------------------------------------------------------------------
# VCF to .beagle
#
if [ ! -e ${panel}.beagle ]
then
  if [ $type == "bed" ] 
  then plink --recode vcf --bfile $panel --out $panel 
  fi
  
  if [ ! -e vcftogenolike.pl ] 
  then ln -s ~/data/Git/Botocudos-scripts/GenoLike/vcftogenolike.pl ./
  fi
  perl vcftogenolike.pl -i $panel.vcf  -rmdamage $rmdamage -o ${panel}.beagle
fi

echo "${panel}.beagle ready"

#--------------------------------------------------------------------------
# Sites file for ANGSD
#
if [ ! -e ${panel}_sites.txt ]
then
  CHRS=($(tail -n +2 $panel.beagle | cut -f1 -d'_' |sort |uniq))
  for chr in ${CHRS[@]}
  do 
    cut -f 1-3 ${panel}.beagle | sed 1d |\
    grep "^${chr}_" | sed 's/_/\t/' |sort -V > ${panel}_${chr}_sites.txt
    angsd sites index ${panel}_${chr}_sites.txt
    cut -f1 ${panel}_${chr}_sites.txt |sort -V |uniq >chr${chr}.txt
  done
  cut -f 1 ${panel}.beagle > ${panel}_ids.txt
fi

#------------------------------------------------------------------------------
# Compute genotype likelihoods for the samples
#
 n=$(cat ${panel}_*_sites.txt | wc -l|cut -f1 -d ' ')
if [ ! -e ${bamlist}_${n}_sites.beagle ]
then
  CHRS=($(tail -n +2 $panel.beagle | cut -f1 -d'_' |sort |uniq))
  for chr in ${CHRS[@]}
  do
    angsd -GL 1 -out ${bamlist}_${n}_${chr}_sites -doGlf 2 -doMajorMinor 3  \
      -bam $bamlist-minQ 35 -minmapQ 30 -trim 5\
      -sites ${panel}_${chr}_sites.txt -rf chr${chr}.txt  -checkbamheaders 0 \
       -nThreads 4  >out_gl_${chr}.txt 2>err_gl_${chr}.txt& #-minInd 1
  done 
  wait
  echo all jobs are done!
  cat ${bamlist}_${n}_*_sites.beagle.gz >${bamlist}_${n}_sites.beagle.gz
  gunzip ${bamlist}_${n}_sites.beagle.gz -c > ${bamlist}_${n}_tmp.beagle
  head -n1 ${bamlist}_${n}_tmp.beagle > header_${bamlist}_${n}.txt
  sed -i '/marker/d' ${bamlist}_${n}_tmp.beagle
  cat header_${bamlist}_${n}.txt ${bamlist}_${n}_tmp.beagle >  ${bamlist}_${n}_sites.beagle
  rm header_${bamlist}_${n}.txt ${bamlist}_${n}_tmp.beagle ${bamlist}_${n}_sites.beagle.gz 
  #sed -i 's/chr//' ${bamlist}_sites.beagle

fi

echo "${bamlist}_${n}_sites.beagle ready"

#------------------------------------------------------------------------------
# Merge genotype likelihoods (samples,panel)
#
if [ ! -e ${panel}_${bamlist}.beagle ]
then
  echo "perl ~/data/Scripts/merge_genos.pl -g2 $panel.beagle \
    -g1 ${bamlist}_${n}_sites.beagle -id ${dir}/GenoLike/ids.txt \
    -o ${panel}_${bamlist}.beagle -homozygous $homo"

  if [ ! -e merge_genos.pl ]
  then ln -s ~/data/Git/Botocudos-scripts/GenoLike/merge_genos.pl ./
  fi 
  perl merge_genos.pl -g2 $panel.beagle \
    -g1 ${bamlist}_${n}_sites.beagle -id ${panel}_ids.txt \
    -o ${panel}_${bamlist}.beagle -homozygous $homo

fi

echo "Panel and samples merged to ${panel}_${bamlist}.beagle"

#/home/ubelix/iee/dcruz/install/NGSadmix -likes ${panel}_8botocudos.beagle \
#-K $k  -o ${panel}_8botocudos_k$k


#if [ ! -e ${panel}_likelihoods.txt]
#then
#  grep like *log |sed 's/.*_k/k/' | sed 's/.log.*=/\t/' | sed "s/ after.*//" |\
#   sort -V  > ${panel}_likelihoods.txt
#fi
