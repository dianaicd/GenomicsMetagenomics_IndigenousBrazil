#!/bin/bash

# You must specify a valid email address!
#SBATCH --mail-user=<user>@<group>.unibe.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="Example Job"

# Runtime and memory
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G

# Partition
##SBATCH --partition=long

# For parallel jobs
##SBATCH --cpus-per-task=8
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


##SBATCH --output=/path/to/outfile
##SBATCH --error=/path/to/errfile

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10

#### Your shell commands below this line ####

