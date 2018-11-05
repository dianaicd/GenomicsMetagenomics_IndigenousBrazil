#!/bin/bash

# You must specify a valid email address!
#SBATCH --mail-user=diana.cruz@iee.unibe.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="bammds"

# Runtime and memory
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=95G

# Partition
##SBATCH --partition=long

# For parallel jobs
##SBATCH --cpus-per-task=8
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


##SBATCH --output=log/out_%A_%a.pca
##SBATCH --error=log/err_%A_%a.pca
#SBATCH --output=out_8.mds
#SBATCH --error=err_8.mds
# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=0-3

#### Your shell commands below this line ####
module add UHTS/Analysis/vcftools/0.1.14
panel=America
#rm *ba*
#ln ~/data/Botocudos/BAM/2017_09_15/*nuc*realigned.bam ./
#ln ~/data/Botocudos/BAM/2017_09_15/*nuc*realigned.bai ./

#for bam in *bam
#do
#  samtools view -H $bam >header.txt
#  sed -i 's/SN:chr/SN:/' header.txt
#  samtools reheader header.txt $bam > Botocudo.$(basename $bam .bam).bam
#  samtools index Botocudo.$(basename $bam .bam).bam
#mv $bam Botocudo.$(basename $bam .hg19_nuc.realigned.bam).bam
#done
#panel=(America America_Oceania_nopapus_noaus America_Oceania_reheaded \
#merged_noduplicates_reheaded)
boto=(MA2391 MA2403 MA2386 MA2381 MA2397 MA2404 MA2388 MA2385 MA2396 \
 MA2401 MA2389 MA2383 MA2387 MA2395 MA2384 MA2400 \
 MA2394 MA2402 MA2398 MA2399 \
 MA2382 MA2392)
for panel in America America_Oceania_nopapus_noaus America_Oceania_reheaded \
merged_noduplicates_reheaded
do
  echo $panel
  #

#MA2402, MA2394, MA2393, MA2400, MA2384, MA2399, MA2382, MA2398, MA2392
  bammds Botocudo.${boto[21]}.bam Botocudo.${boto[20]}.bam  \
  Botocudo.${boto[19]}.bam Botocudo.${boto[18]}.bam \
  Botocudo.${boto[17]}.bam Botocudo.${boto[16]}.bam \
  Botocudo.${boto[15]}.bam Botocudo.${boto[14]}.bam \
  $panel.vcf --dim1 1 --dim2 2 --output ${panel}_8_d1_d2.csv

   bammds Botocudo.${boto[21]}.bam \
   Botocudo.${boto[20]}.bam  Botocudo.${boto[19]}.bam \
   Botocudo.${boto[18]}.bam Botocudo.${boto[17]}.bam \
   Botocudo.${boto[16]}.bam Botocudo.${boto[15]}.bam \
   Botocudo.${boto[14]}.bam $panel.vcf --dim1 2 --dim2 3 \
  --output ${panel}_8_d2_d3.csv
done
