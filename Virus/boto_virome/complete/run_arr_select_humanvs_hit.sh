#!/bin/bash
#SBATCH --job-name human_virus_hits
#SBATCH --partition ax-normal
#SBATCH --time 1:00:00
#SBATCH --mem-per-cpu 1G
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --error human_virus_hits/hvs_hits_job%A_%a.err
#SBATCH --output human_virus_hits/hvs_hits_job%A_%a.out
#SBATCH --array=0-23%

# February 13 2020
# Run the script count_reads.py in all samples to get the best hits that are against a human virus.

# Python 3
# conda activate myenv_python3

# Path to python script
select_hvs_hits=/scratch/axiom/FAC/FBM/DBC/amalaspi/virome/yarizmen/virome/Boto_virome/python_scripts/select_hvirus_hits.py
# Path to file with human virus taxonomic IDs
human_virus_taxids=/scratch/axiom/FAC/FBM/DBC/amalaspi/virome/yarizmen/virome/Boto_virome/true_complete/DIAMOND_output/HSapiensVirus_taxids

# Samples
MNS=('MN00010' 'MN00013' 'MN00016' 'MN00019' 'MN00021' 'MN00022' 'MN00023' 'MN00039' 'MN0003' 'MN00045' 'MN00056' 'MN00064' 'MN00066' 'MN00067' 'MN00068' 'MN00069' 'MN0008' 'MN0009' 'MN00118' 'MN00119' 'MN00316' 'MN00346' 'MN01701' 'MN1943')

# Set input file
alg_file=${MNS[$SLURM_ARRAY_TASK_ID]}_alg.txt

output=human_virus_hits/${MNS[$SLURM_ARRAY_TASK_ID]}_hvs_hits.csv

# Run count_reads.py
python $select_hvs_hits $human_virus_taxids $alg_file $output

