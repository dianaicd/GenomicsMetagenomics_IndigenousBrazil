# Script to run contamMix
# Author: Diana Cruz
# July 23rd, 2018

#-------------------------------------------------------------------------------
# Required modules:
# mafft
# R
# samtools
# I have no spceific reasons to load the following versions
module add SequenceAnalysis/MultipleSequenceAlignment/mafft/7.310
module add R/3.4.2
module add UHTS/Analysis/samtools/1.4

# bam [file.bam] is the name of a bam file. E.g., Quack.bam
# rmtrans [yes\no]indicates whether to ignore transitions when running contamMix
bam=$1
rmtrans=$2

ind=$(basename $bam .bam)
mito=MT

mkdir $mito
cd $mito

#-------------------------------------------------------------------------------
# The first step is to create a
# consensus sequence for the mtDNA of the sample
# (fasta format)
if [ ! -e ${ind}_${mito}.fa ]
then
  echo "Consensus..."
  samtools view -b ../$bam $mito > ${ind}_${mito}.bam
  samtools index ${ind}_${mito}.bam

  angsd -doFasta 2 -doCounts 1 -i ${ind}_${mito}.bam -out ${ind}_${mito}
  gunzip ${ind}_${mito}.fa.gz
  bwa index ${ind}_${mito}.fa
  samtools faidx ${ind}_${mito}.fa
  picard-tools CreateSequenceDictionary REFERENCE=${ind}_${mito}.fa OUTPUT=${ind}_${mito}.dict
fi
#-------------------------------------------------------------------------------
# Once we have the consensus, we should realign
# the reads to this fasta "reference"
if [ ! -e ${ind}_${mito}_consensus.bam ]
then
  echo "Realigning..."
  # _consensus.truncated so that it works with mapping_aDNA.sh
  bedtools bamtofastq -i ${ind}_${mito}.bam -fq ${ind}_${mito}_consensus.truncated
  gzip ${ind}_${mito}_consensus.truncated

  echo -e "ID\tData\tMAPQ\tLB\tPL\tSM" > mito.txt
  paste <(ls ${ind}_${mito}_consensus.truncated.gz | cut -d'/' -f3 | cut -d'.' -f1) \
  <(ls ${ind}_${mito}_consensus.truncated.gz) \
  <(ls ${ind}_${mito}_consensus.truncated.gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f-3) \
  <(ls ${ind}_${mito}_consensus.truncated.gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f1 | sed "s/$/_${mito}_consensus/") | \
  awk '{print $1,$2,"30",$3,"ILLUMINA",$4}' OFS='\t' >> mito.txt

  echo "------------------Input file for mapping_aDNA.sh------------------------"
  cat mito.txt
  mapping_aDNA.sh --ref ${ind}_${mito}.fa -i mito.txt --skip1 1
fi

#------------------------------------------------------------------------------
# Multiple alignment of the contaminants and the consensus
if [ ! -e 311_${ind[$i]}_aligned.fasta ]
then
  echo "–--------------------------Multiple alignment--------------------------"

  cat ~/archive/Contamination/311humans.fasta ${ind}_${mito}.fa > 311_${ind}.fasta
  mafft --auto 311_${ind}.fasta > 311_${ind}_aligned.fasta
fi

#-------------------------------------------------------------------------------
# Finally, run contamMix
if [ "$rmtrans" == "rmtrans" ]
then
  echo "-----------------Running contamMix, removing transitions---------------"
  ~/data/Scripts/contamMix/exec/estimate.R --samFn ${ind}_${mito}_consensus.bam  --malnFn 311_${ind}_aligned.fasta  --figure ${ind}_rmtrans_fig --nIter 100000 --alpha 0.1 --transverOnly --saveData ${ind}_rmtrans_data  > out_rmtrans_${ind}.txt
else
  echo "--------------Running contamMix, using all the mt reads----------------"
  ~/data/Scripts/contamMix/exec/estimate.R --samFn ${ind}_${mito}_consensus.bam  --malnFn 311_${ind}_aligned.fasta  --figure ${ind}_fig --nIter 100000 --alpha 0.1 --saveData ${ind}_data > out_${ind}.txt
fi
