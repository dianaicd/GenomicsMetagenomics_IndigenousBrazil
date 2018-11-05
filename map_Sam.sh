#!/bin/bash

# You must specify a valid email address!
#SBATCH --mail-user=diana.cruz@iee.unibe.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="Mapping"

# Runtime and memory
#SBATCH --time=95:00:00
#SBATCH --mem-per-cpu=90G

# Partition
##SBATCH --partition=empi

# For parallel jobs
##SBATCH --cpus-per-task=8
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


#SBATCH --output=log/out_second_%A_%a.txt
#SBATCH --error=log/err_second_%A_%a.txt

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
#SBATCH --array=1,10,16,17%4

#### Your shell commands below this line ####
ind=($(cat ma_ids.txt))
f=($(cat fastq.txt))

mkdir ${ind[$SLURM_ARRAY_TASK_ID]}
cd ${ind[$SLURM_ARRAY_TASK_ID]}

~/install/mapping_aDNA.sh --fastq1 ${f[$SLURM_ARRAY_TASK_ID]} --base ${ind[$SLURM_ARRAY_TASK_ID]} --ref /home/ubelix/iee/dcruz/data/References/Human/hg19/chromosomes_1-22_XY_mito.fasta -p 9 --ubelix

