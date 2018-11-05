#!/bin/bash

# You must specify a valid email address!
#SBATCH --mail-user=diana.cruz@iee.unibe.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="ngsadmix-all"

# Runtime and memory
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=60G

# Partition
##SBATCH --partition=long

# For parallel jobs
##SBATCH --cpus-per-task=8
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


#SBATCH --output=log/out_%A_%a.all
#SBATCH --error=log/err_%A_%a.all

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
#SBATCH --array=0-19%20

#### Your shell commands below this line ####
k=($(seq 2 20))
./workflow_genolike.sh /home/ubelix/iee/dcruz/data/Merged/500k \
 America_Oceania_homozygous ${k[$SLURM_ARRAY_TASK_ID]} 
