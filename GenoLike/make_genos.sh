#!/usr/bin/bash

# Path to bam files:
# Directory where we can find the VCF panel
dir=/home/dcruzdva/archive/Panels/Merged/500k
# Name of the panel
panel=88ind_nodamage
# Name for the output
name=24ind
# Selected individuals
# pass it as
# "$(echo ${selected[@]})"
selected=($(cut -f4 -d',' /home/dcruzdva/archive/Botocudo/database_names.csv |grep -v 'MN$' |grep -v MA |grep -v Quack))
# MAke only homozygous sites? (usually no)
homo=no
# Path to bam files
# soft links to final bams
bam_path=/home/dcruzdva/archive/BAM/Botocudos/2018_10_26/link_final
# remove damaged sites from the panel? (usually yes)
rmdamage=yes

workflow_genolike.sh $dir $panel $name "$(echo ${selected[@]})" $homo $bam_path $rmdamage