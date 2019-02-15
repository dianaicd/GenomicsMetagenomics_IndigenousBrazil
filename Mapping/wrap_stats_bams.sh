


# Script to wrap up summary statistics for a mapping

bam=$1
name=$2
organism=$3

if [ ! -e ${bam}.bam ]
then
  echo "${bam}.bam does not exist."
  exit 1
fi

if [ ! -e ${bam}.bam.bai ]
then
  samtools index ${bam}.bam
fi

libs=($(samtools view -H ${bam}.bam |grep "^@RG" |cut -f 3 | sed 's/LB://'))
nlibs=$(echo ${#libs[*]})

echo "The file $bam.bam contains $nlibs library(ies)."

for library in ${libs[@]}
do
  echo "Working on $library"
  IDs=($(samtools view -H ${bam}.bam |grep "^@RG" | grep $library |cut -f 2 | sed 's/ID://'))
  nIDs=$(echo ${#IDs[*]})
  # Verify that .settings files exist
  for file in "${IDs[@]}"
  do
    echo $file
    if [ ! -e $file.settings ]
      then
        echo "Could not find ${file}.settings"
          exit 1
    fi
  done

  #-----------------------------------------------------------------------------
  # Make a .settings file per library
  echo "Library $library was split into $nIDs file(s)."
  ~/data/Scripts/summary_settings.sh $library "$(echo ${IDs[@]})"
done

#-------------------------------------------------------------------------------
# Make a .settings file per sample
if [ ! -e ${bam}.settings ]
then
  ~/data/Scripts/summary_settings.sh $bam "$(echo ${libs[@]})"
fi

settings=$(~/data/Scripts/summary_settings.sh $bam "$(echo ${libs[@]})")
#-------------------------------------------------------------------------------
if [ ! -e ${bam}_depth.txt ]
then
  depth.py $bam.bam > ${bam}_depth.txt
fi

chrs=($(samtools view -H ${bam}.bam |grep "^@SQ" |sed 's/SN://' |cut -f 2))
seq_retained_reads=$(echo $settings |cut -f 4 -d ' ')
endogenous=$(~/data/Scripts/count_hits.sh $bam "$(echo ${chrs[@]})" ${seq_retained_reads} endo)

if [ $organism = "Human" ]
then
  mito_chr=(MT)
  nuclear_chr=( "${chrs[@]/$mito_chr}")


  nuclear=$(~/data/Scripts/count_hits.sh $bam "$(echo ${nuclear_chr[@]})" ${seq_retained_reads} nuclear)
  mitochondrial=$(~/data/Scripts/count_hits.sh $bam "$(echo ${mito_chr[@]})" ${seq_retained_reads} MT)

  echo "Target Sample Library lib_type seq_reads_se seq_trash_se seq_trash_se_frac \
seq_retained_reads seq_retained_nts seq_retained_length hits_raw_endogenous \
hits_raw_frac_endogenous hits_clonality_endogenous hits_unique_endogenous \
hits_unique_frac_endogenous hits_coverage_endogenous hits_length_endogenous \
hits_raw_mitochondrial hits_raw_frac_mitochondrial hits_clonality_mitochondrial \
hits_unique_mitochondrial hits_unique_frac_mitochondrial \
hits_coverage_mitochondrial hits_length_mitochondrial hits_raw_nuclear \
hits_raw_frac_nuclear hits_clonality_nuclear hits_unique_nuclear \
hits_unique_frac_nuclear hits_coverage_nuclear hits_length_nuclear"  >${name}.summary

  echo $bam $bam $bam SE $settings $endogenous $mitochondrial $nuclear >>${name}.summary
else
  echo "Target Sample Library lib_type seq_reads_se seq_trash_se seq_trash_se_frac \
seq_retained_reads seq_retained_nts seq_retained_length hits_raw \
hits_raw_frac hits_clonality hits_unique \
hits_unique_frac hits_coverage hits_length"  >${name}.summary

echo $bam $bam $bam SE $settings $endogenous >>${name}.summary

fi
