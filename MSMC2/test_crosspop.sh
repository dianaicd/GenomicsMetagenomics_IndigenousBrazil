#!/bin/bash


# Don't forget to pass the correct --account 
# for each job! (${PI}_${project})
#SBATCH --account=amalaspi_americas

# You must specify a valid email address!
#SBATCH --mail-user=diana-ivette.cruz-davalos@unil.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="crosspop"

# Runtime and memory
# default: 12h (normal, ax-normal), 5d (long, ax-long)
# limit: 24h (normal, ax-normal), 10d (long, ax-long)
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=100G

# Partition
# Options are:
# cpu, gpu
#SBATCH --partition=cpu

# For parallel jobs
#SBATCH --cpus-per-task=8
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


#SBATCH --output=logs/test_crosspop.out 
#SBATCH --error=logs/test_crosspop.err

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10
##SBATCH --array=1-22

./msmc2 -p 1*1+20*1+1*5 -t 8 -I 4-10,4-11,5-10,5-11 -o crosspop/rep2_Botocudo_MayaG_paramB_crosspop bootstrapped_2/bootstrap_multihetsep.chr*