

#[1] "Target"
x=$1
name=$2
#cd $x
  # Name of mitochondrial chromosome
mito=MT

  #_Human_rmdup_realign_calmd.bam
  #target=$(basename $x _Human_rmdup_realign_calmd.bam)

#------------------------------------------------------------------------------#
# Verify that the inputs are present
if [ ! -e $x.bam ]
  then
    echo "Could not find ${x}.bam"
    exit 1
fi
if [ ! -e $x.settings ]
  then
    echo "Could not find ${x}.settings"
    exit 1
fi
if [ ! -e ${x}.bam.bai ]
    then
      samtools index ${x}.bam
fi

if [ ! -e $x.coverage ]
  then
    paleomix coverage ${x}.bam ${x}.coverage
fi

if [ ! -e ${x}_flagstat.txt ]
  then
    samtools flagstat ${x}.bam > ${x}_flagstat.txt
fi

if [ ! -e ${x}_idxstats.txt ]
  then
    samtools idxstats ${x}.bam > ${x}_idxstats.txt
fi

# For simplicity, we will keep Sample and Library the same as Target
#"Sample"
#"Library"

# We are assuming it is a single-end type (SE)
#[4] "lib_type"

all_chr=($(samtools view -H ${x}.bam |grep '@SQ' |cut -f 2 |sed 's/SN://'))


#------------------------------------------------------------------------------#
# Statistics refering to total reads and trimming
~/data/Scripts/summary_settings.sh
#------------------------------------------------------------------------------#
# Endogenous reads
# Endogenous=all contigs defined in the reference
# raw = before removing duplicates
#"hits_raw.endogenous."
# Notice that this final bam file includes duplicates, but they are marked
hits_raw_endogenous=$(grep mapped ${x}_flagstat.txt| cut -f 1 -d ' ')
echo "hits_raw_endogenous = $hits_raw_endogenous"
#"hits_raw_frac.endogenous."
hits_raw_frac_endogenous=$(echo "$hits_raw_endogenous / $seq_retained_reads" \
| bc -l)
echo "hits_raw_frac_endogenous = $hits_raw_frac_endogenous"
# clonality = clonal duplicates
#[13] "hits_clonality.endogenous."
hits_clonality_endogenous=$(grep duplicates ${x}_flagstat.txt| cut -f1 -d ' ')
echo "hits_clonality_endogenous = $hits_clonality_endogenous"
# unique = total reads excluding clonal duplicates
#"hits_unique.endogenous."
hits_unique_endogenous=$(grep -v '#' ${x}.coverage |head -n 2 | tail -n1 |cut -f 6)
echo "hits_unique_endogenous = $hits_unique_endogenous"
#"hits_unique_frac.endogenous."
hits_unique_frac_endogenous=$(echo \
"$hits_unique_endogenous / $seq_retained_reads" | bc -l)
echo "hits_unique_frac_endogenous = $hits_unique_frac_endogenous"
# Easy to parse from paleomix output
#[16] "hits_coverage.endogenous."
hits_coverage_endogenous=$(grep -v '#' ${x}.coverage |head -n 2 | tail -n1 |\
  cut -f 14)
echo "hits_coverage_endogenous = $hits_coverage_endogenous"
#"hits_length.endogenous."
hits_length_endogenous=$(samtools view ${x}.bam |perl ~/data/Scripts/length.pl -o ${x}_endo_length.txt -type endo)
echo "hits_length_endogenous = $hits_length_endogenous"
# Repeat for mitochondiral DNA

#------------------------------------------------------------------------------#
# Mitochondrial reads
#"hits_raw.mitochondrial."
hits_raw_mitochondrial=$(samtools view \
  -c ${x}.bam $mito)
echo "hits_raw_mitochondrial = $hits_raw_mitochondrial"
#[22] "hits_raw_frac.mitochondrial."
hits_raw_frac_mitochondrial=$(echo "$hits_raw_mitochondrial / \
$seq_retained_reads" |bc -l)
echo "hits_raw_frac_mitochondrial = $hits_raw_frac_mitochondrial"
#"hits_clonality.mitochondrial."
hits_clonality_mitochondrial=$(samtools view -f1024 \
-c ${x}.bam $mito)
echo "hits_clonality_mitochondrial = $hits_clonality_mitochondrial"
#"hits_unique.mitochondrial."
hits_unique_mitochondrial=$(grep -v '#' ${x}.coverage | grep -P '\*\t\*' |\
  grep $mito |cut -f 6)
echo "hits_unique_mitochondrial = $hits_unique_mitochondrial"
#[25] "hits_unique_frac.mitochondrial."
hits_unique_frac_mitochondrial=$(echo "$hits_unique_mitochondrial / \
$seq_retained_reads" | bc -l)
echo "hits_unique_frac_mitochondrial = $hits_unique_frac_mitochondrial"
#"hits_coverage.mitochondrial."
hits_coverage_mitochondrial=$(grep -v '#' ${x}.coverage | grep -P '\*\t\*' |\
  grep $mito |cut -f 14)
echo "hits_coverage_mitochondrial = $hits_coverage_mitochondrial"
#"hits_length.mitochondrial."
hits_length_mitochondrial=$(samtools view ${x}.bam $mito |perl ~/data/Scripts/length.pl -o ${x}_mito_length.txt -type MT)
echo "hits_length_mitochondrial = $hits_length_mitochondrial"


#------------------------------------------------------------------------------#
# Nuclear reads
# Nuclear = all contigs in the reference, excluding mtDNA
#[28] "hits_raw.nuclear."
nuc_chr=($(samtools view -H ${x}.bam \
 |grep "@SQ" |sed 's/.*SN://' |sed 's/\tLN:.*//' | grep -v $mito))
echo "nuc_chr = ${nuc_chr[@]}"
hits_raw_nuclear=$(samtools view -c ${x}.bam \
 $(echo ${nuc_chr[@]}))
echo "hits_raw_nuclear = $hits_raw_nuclear"
#"hits_raw_frac.nuclear."
hits_raw_frac_nuclear=$(echo "$hits_raw_nuclear / $seq_retained_reads" | bc -l)
echo "hits_raw_frac_nuclear = $hits_raw_frac_nuclear"
#"hits_clonality.nuclear."
hits_clonality_nuclear=$(samtools view -c -f1024 \
${x}.bam $(echo ${nuc_chr[@]})) #$(seq 1 22) X Y)
echo "hits_clonality_nuclear = $hits_clonality_nuclear"
#[31] "hits_unique.nuclear."
hits_unique_nuclear=$(($hits_raw_nuclear - $hits_clonality_nuclear))
echo "hits_unique_nuclear = $hits_unique_nuclear"
#"hits_unique_frac.nuclear."
hits_unique_frac_nuclear=$(echo "$hits_unique_nuclear / $seq_retained_reads" | bc -l)
echo "hits_unique_frac_nuclear = $hits_unique_frac_nuclear"
#"hits_coverage.nuclear."
bases=$(grep -v '#' ${x}.coverage |\
grep -P "\*\t\*\t\*" -  |cut -f 11 )
bases_mito=$(grep -v '#' ${x}.coverage |\
grep -P "\*\t\*\t${mito}" -  |cut -f 11)
bases=$(echo "${bases} - ${bases_mito} "|bc)
echo "bases = $bases"
bases=$(echo ${bases[@]} | sed 's/ /\+/g' - |bc)
echo "bases = $bases"
l_genome=$(grep -v '#' ${x}.coverage |\
grep -P "\*\t\*" - | grep -v $mito |cut -f 5 | grep -v Size)
echo "l_genome = $l_genome"
l_genome=$(echo ${l_genome[@]}| sed 's/ /\+/g' - |bc )
hits_coverage_nuclear=$(echo "$bases / $l_genome" |bc -l )
echo "l_genome = $l_genome"
#[34] "hits_length.nuclear."
hits_length_nuclear=$(samtools view ${x}.bam | perl ~/data/Scripts/length.pl -o ${x}_nuclear_length.txt -type nuclear)
echo "hits_length_nuclear = $hits_length_nuclear"


#------------------------------------------------------------------------------#
# Some ratios that might be of interest
#"ratio_reads.nuc.mito."
ratio_reads_nuc_mito=$(echo "$hits_unique_nuclear / \
$hits_unique_mitochondrial" | bc -l)
echo "ratio_reads_nuc_mito = $ratio_reads_nuc_mito"
#[19] "ratio_genome.mito.nuc."
l_mito=$(grep -v '#' ${x}.coverage |\
grep -P "\*\t\*" - | grep $mito |cut -f 5 )
echo "l_mito = $l_mito"
bases_mito=$(grep -v '#' ${x}.coverage |\
grep -P "\*\t\*" - | grep $mito |cut -f 11)
echo "bases_mito = $bases_mito"
ratio_genome_mito_nuc=$(echo "$bases_mito / $bases / ($l_mito) / (2 * $l_genome)" | bc -l)
echo "ratio_genome_mito_nuc = $ratio_genome_mito_nuc"
#"ratio_genome.nuc.mito."
ratio_genome_nuc_mito=$(echo "$bases / $bases_mito / (2 * $l_genome) / $l_mito" | bc -l)
echo "ratio_genome_nuc_mito = $ratio_genome_nuc_mito"


#------------------------------------------------------------------------------#
# Print to output file
echo "Target Sample Library lib_type seq_reads_se seq_trash_se seq_trash_se_frac \
seq_retained_reads seq_retained_nts seq_retained_length hits_raw_endogenous \
hits_raw_frac_endogenous hits_clonality_endogenous hits_unique_endogenous \
hits_unique_frac_endogenous hits_coverage_endogenous hits_length_endogenous \
ratio_reads_nuc_mito ratio_genome_mito_nuc ratio_genome_nuc_mito \
hits_raw_mitochondrial hits_raw_frac_mitochondrial hits_clonality_mitochondrial \
hits_unique_mitochondrial hits_unique_frac_mitochondrial \
hits_coverage_mitochondrial hits_length_mitochondrial hits_raw_nuclear \
hits_raw_frac_nuclear hits_clonality_nuclear hits_unique_nuclear \
hits_unique_frac_nuclear hits_coverage_nuclear hits_length_nuclear" | \
sed 's/ /\t/g' >$name.summary

echo "$x $x $x SE $seq_reads_se $seq_trash_se $seq_trash_se_frac \
$seq_retained_reads $seq_retained_nts $seq_retained_length $hits_raw_endogenous \
$hits_raw_frac_endogenous $hits_clonality_endogenous $hits_unique_endogenous \
$hits_unique_frac_endogenous $hits_coverage_endogenous $hits_length_endogenous \
$ratio_reads_nuc_mito $ratio_genome_mito_nuc $ratio_genome_nuc_mito \
$hits_raw_mitochondrial $hits_raw_frac_mitochondrial $hits_clonality_mitochondrial \
$hits_unique_mitochondrial $hits_unique_frac_mitochondrial \
$hits_coverage_mitochondrial $hits_length_mitochondrial $hits_raw_nuclear \
$hits_raw_frac_nuclear $hits_clonality_nuclear $hits_unique_nuclear \
$hits_unique_frac_nuclear $hits_coverage_nuclear $hits_length_nuclear" | \
sed 's/ /\t/g' >>$name.summary
