#!/bin/bash

# You must specify a valid email address!
#SBATCH --mail-user=diana.cruz@iee.unibe.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="pca-all"

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

##SBATCH --output=out.pca
##SBATCH --error=err.pca
#SBATCH --output=log/out_%A_%a.pca
#SBATCH --error=log/err_%A_%a.pca

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
#SBATCH --array=0-22

#### Your shell commands below this line ####
module add EcologyEvolution/eigensoft/6.1.4

mybam=($(ls *bam))
old_panel=merged_noduplicates_pop.vcf
mypanel=merged_noduplicates
ref=/home/ubelix/iee/dcruz/data/References/Human/hg19/chromosomes_1-22_XY_Un.fasta

#bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' all_samples_merged.vcf.gz | bcftools view -m2 -M2 -v snps -o ${mypanel}.vcf

#plink --recode --vcf $old_panel --out $mypanel --geno 0.05

#plink --file $mypanel --make-bed --out $mypanel
#
#cut -f 1,4 ${mypanel}.map | sed 's/^/chr/' | sed 's/chr23/X/' | sed 's/chr24/Y/' > ${mypanel}.pos

#awk '{printf $1" "$2" "$3" "$4" "$5" "$6; for (i=7; i<=NF; i=i+2){x=int(.5 + rand()); printf (" "$(i+x)" "$(i+x));} print "" }' ${mypanel}.ped > ${mypanel}.rand.ped
#cp ${mypanel}.map ${mypanel}.rand.map

./make_evec_per_sample_Diana.sh ${mybam[$SLURM_ARRAY_TASK_ID]} $mypanel $ref
plink --pca --file $(basename ${mybam[$SLURM_ARRAY_TASK_ID]} .bam)\.filtSites_MERGED --out $(basename ${mybam[$SLURM_ARRAY_TASK_ID]} .hg19_nuc.realigned.bam)
