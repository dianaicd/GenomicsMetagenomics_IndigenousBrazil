#!/bin/bash

# Don't forget to pass the correct --account 
# for each job! (${PI}_${project})
#SBATCH --account=amalaspi_americas

# You must specify a valid email address!
#SBATCH --mail-user=diana-ivette.cruz-davalos@unil.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="msmc2"

# Runtime and memory
# default: 12h (normal, ax-normal), 5d (long, ax-long)
# limit: 24h (normal, ax-normal), 10d (long, ax-long)
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=12G

# Partition
# Options are:
# cpu, gpu
#SBATCH --partition=cpu

# For parallel jobs
#SBATCH --cpus-per-task=8
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


#SBATCH --output=log/msmc_single_%A-%a.out 
#SBATCH --error=log/msmc_single_%A-%a.err

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10
#SBATCH --array=0-3
#### Your shell commands below this line ####

threads=8
PopNames=(Karitiana Botocudo Maya Surui)
Karitiana="0,1,2,3"
Botocudo="4,5" # remove MN00056
Maya="8,9,10,11"
Surui="12,13,14,15"
PopIndices=($Karitiana $Botocudo $Maya $Surui)

# run it once per pop
# msmc2 -t $threads -I $Karitiana -o Karitiana_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Botocudo -o Botocudo_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Maya -o Maya_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Surui -o Surui_4hap multihetsep/chr*.txt
echo $SLURM_ARRAY_TASK_ID
echo ${PopIndices[$SLURM_ARRAY_TASK_ID]}
echo ${PopNames[$SLURM_ARRAY_TASK_ID]}_4hap

~/install/msmc2/build/release/msmc2 -t $threads -I ${PopIndices[$SLURM_ARRAY_TASK_ID]} -o ${PopNames[$SLURM_ARRAY_TASK_ID]}_4hap multihetsep/chr*.txt
