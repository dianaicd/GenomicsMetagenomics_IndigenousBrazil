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

# Directory where we can find the VCF panel
dir=$1
# Name of the panel
panel=$2
# Name for the output
name=$3
# Selected individuals
# pass it as
# "$(echo ${selected[@]})"
selected=( "$4" ) #(MA2394 MA2387 MA2384 MA2400 MA2395 MA2402 MA2398 MA2399 MA2382 MA2392)
# MAke only homozygous sites? (usually no)
homo=$5
# Path to bam files
bam_path=$6

# remove damaged sites from the panel? (usually yes)
rmdamage=$7
#k=$4

echo "Parameters:
  dir=$dir
  panel=$panel
  name=$name
  selected=${selected[@]}
  homo=$homo
  bam_path=$bam_path
  rmdamage=$rmdamage
"

#-------------------------------------------------------------------------
# VCF to .beagle
#
if [ ! -e ${panel}.beagle ]
then
  perl ~/data/Scripts/vcftogenolike.pl -i ${dir}/$panel.vcf -rmdamage $rmdamage \
  -o ${panel}.beagle
fi

echo "${panel}.beagle ready"

#-------------------------------------------------------------------------
# List with bam files
#
if [ ! -e $name.bam.list ]
then
  touch $name.bam.list
  for i in ${selected[@]}
    do
    #
    # Notice the complete name of your bam file
    #
    ls ${bam_path}/${i}.bam >>$name.bam.list
  done
fi

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
if [ ! -e ${name}_${n}_sites.beagle ]
then
  CHRS=($(tail -n +2 $panel.beagle | cut -f1 -d'_' |sort |uniq))
  for chr in ${CHRS[@]}
  do
    angsd -GL 1 -out ${name}_${n}_${chr}_sites -doGlf 2 -doMajorMinor 3  \
      -bam $name.bam.list -minQ 35 -minmapQ 30 -trim 5\
      -sites ${panel}_${chr}_sites.txt -rf chr${chr}.txt -nThreads 4  & #-minInd 1
  done 
  { sleep 5; echo waking up after 5 seconds; } &
  { sleep 1; echo waking up after 1 second; } &
  wait
  echo all jobs are done!
  cat ${name}_${n}_*_sites.beagle.gz >${name}_${n}_sites.beagle.gz
  gunzip ${name}_${n}_sites.beagle.gz -c > ${name}_${n}_tmp.beagle
  head -n1 ${name}_${n}_tmp.beagle > header_${name}_${n}.txt
  sed -i '/marker/d' ${name}_${n}_tmp.beagle
  cat header_${name}_${n}.txt ${name}_${n}_tmp.beagle >  ${name}_${n}_sites.beagle
  rm header_${name}_${n}.txt ${name}_${n}_tmp.beagle ${name}_${n}_sites.beagle.gz 
  #sed -i 's/chr//' ${name}_sites.beagle

fi

echo "${name}_${n}_sites.beagle ready"

#------------------------------------------------------------------------------
# Merge genotype likelihoods (samples,panel)
#
if [ ! -e ${panel}_${name}.beagle ]
then
  echo "perl ~/data/Scripts/merge_genos.pl -g2 $panel.beagle \
    -g1 ${name}_${n}_sites.beagle -id ${dir}/GenoLike/ids.txt \
    -o ${panel}_${name}.beagle -homozygous $homo"

  perl ~/data/Scripts/merge_genos.pl -g2 $panel.beagle \
    -g1 ${name}_${n}_sites.beagle -id ${panel}_ids.txt \
    -o ${panel}_${name}.beagle -homozygous $homo

fi

echo "Panel and samples merged to ${panel}_${name}.beagle"

#/home/ubelix/iee/dcruz/install/NGSadmix -likes ${panel}_8botocudos.beagle \
#-K $k  -o ${panel}_8botocudos_k$k


#if [ ! -e ${panel}_likelihoods.txt]
#then
#  grep like *log |sed 's/.*_k/k/' | sed 's/.log.*=/\t/' | sed "s/ after.*//" |\
#   sort -V  > ${panel}_likelihoods.txt
#fi
