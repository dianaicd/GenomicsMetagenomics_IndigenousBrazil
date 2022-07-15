#!/bin/bash


# Don't forget to pass the correct --account 
# for each job! (${PI}_${project})
#SBATCH --account=amalaspi_americas

# You must specify a valid email address!
#SBATCH --mail-user=diana-ivette.cruz-davalos@unil.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="mask please"

# Runtime and memory
# default: 12h (normal, ax-normal), 5d (long, ax-long)
# limit: 24h (normal, ax-normal), 10d (long, ax-long)
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=2G

# Partition
# Options are:
# cpu, gpu
#SBATCH --partition=cpu

# For parallel jobs
#SBATCH --cpus-per-task=8
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


#SBATCH --output=logs/make_mask_%A-%a.out 
#SBATCH --error=logs/make_mask_%A-%a.err

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10
#SBATCH --array=1-22

# <!-- # merge and look hets
#-----------------------------------------------------------------#
# Mask description
# The chromosome (can be an arbitrary string, but has to be the same for all rows in a file).
# 1. chr
# 2.# The position on the chromosome.
# 3.# The number of called homozygous sites since the last segregating site, 
#   which includes the given location. This number must always be greater than zero 
#  and cannot be larger than the difference between the current position and the previous position.
# 4.# The ordered and phased alleles of the multiple haplotypes. 
#    If the phasing is unknown or only partially known, multiple phasings can be given, 
#    separated by a comma to indicate the different possibilities 
#   (see the second-last line in the example). Unknown alleles can be indicated by “?”, 
#   but they can also simply be left out and expressed by a reduced number of called sites 
#   in the line of the next heterozygous site. -->
#----------------------------------------------------------------#
#----------------------------------------------------------------#
# make masks per chr per ind
# mkdir msmc2_masks
all_inds=( KaritianaA KaritianaB MN0008 MN00056 MayaG MayaH SuruiA SuruiB )

# actually what they want is the positions of all genotypes that were called
chr=$SLURM_ARRAY_TASK_ID

for ind in ${all_inds[@]}
do
    echo $ind
    # for chr in {1..22}
    # do
    echo $chr
    vcf=Final_phased/$ind.chr$chr.phased.vcf.gz
    #vcf=vcf/chr$chr/$ind.chr$chr.vcf.gz
    mask=msmc2_masks/mask_$ind.chr$chr.bed.gz
    less $vcf | grep -v "#" | awk '{print $1,$2-1,$2}' |gzip > $mask &
    # done
done

wait
#----------------------------------------------------------------#
# merge inds
module load gcc/9.3.0 python/3.5.10


all_inds=( KaritianaA KaritianaB MN0008 MN00056 MayaG MayaH SuruiA SuruiB )

# for chr in {1..22}
# do
    mask_command=$( for ind in  ${all_inds[@]} ; do echo "--mask msmc2_masks/mask_$ind.chr$chr.bed.gz" ; done)
    all_vcf=$( for ind in  ${all_inds[@]} ; do echo "Final_phased/$ind.chr$chr.phased.vcf.gz"; done)

    ./generate_multihetsep.py \
    --chr $chr \
        $mask_command \
        --mask masks/chr$chr.bed.gz \
        $all_vcf \
        > input_files_msmc2/chr${chr}.txt

# done

#----------------------------------------------------------------#
#
# threads=8
# Karitiana="0,1,2,3"
# Botocudo="4,5,6,7"
# Maya="8,9,10,11"
# Surui="12,13,14,15"

# # run it once per pop
# msmc2 -t $threads -I $Karitiana -o Karitiana_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Botocudo -o Botocudo_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Maya -o Maya_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Surui -o Surui_4hap multihetsep/chr*.txt

# # cross pop
# msmc2 -t $threads -o cross_Karitiana_Botocudo -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 multihetsep/chr*.txt
# msmc2 -t $threads -o cross_Karitiana_Maya \
#     -I 0-8,0-9,0-10,0-11,1-8,1-9,1-10,1-11,2-8,2-9,2-10,2-11,3-8,3-9,3-10,3-11 multihetsep/chr*.txt
# msmc2 -t $threads -o cross_Karitiana_Surui \
# -I 0-12,0-13,0-14,0-15,1-12,1-13,1-14,1-15,2-12,2-13,2-14,2-15,3-12,3-13,3-14,3-15 multihetsep/chr*.txt
# msmc2 -t $threads -o cross_Botocudo_Maya \
# -I 4-8,4-9,4-10,4-11,5-8,5-9,5-10,5-11,6-8,6-9,6-10,6-11,7-8,7-9,7-10,7-11 multihetsep/chr*.txt
# msmc2 -t $threads -o cross_Botocudo_Surui º
# -I 4-12,4-13,4-14,4-15,5-12,5-13,5-14,5-15,6-12,6-13,6-14,6-15,7-12,7-13,7-14,7-15 multihetsep/chr*.txt
# msmc2 -t $threads -o cross_Maya_Surui -I 8-12,8-13,8-14,8-15,9-12,9-13,9-14,9-15,10-12,10-13,10-14,10-15,11-12,11-13,11-14,11-15 multihetsep/chr*.txt