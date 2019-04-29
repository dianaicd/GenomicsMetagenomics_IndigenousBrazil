#!/bin/sh

# Suggested pipeline in ATLAS

#ln -s ~/install/atlas/atlas ./
ln -s ~/data/Git/Botocudos-scripts/Heterozygosity/plot_theta_altas.R ./ 
ref=/archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa
conserved=~/archive/ConservedSites/UCNE/Conserved_Homo_Sapiens_Vertebrata_Human_69_all.bed
bam=MN0008.bam 

atlas_het()
{
    bam=$1
    sex=$2
    # Pipeline for BAM files with single-end data:
    #-------------------------------------------------------------------------#
    # Split the read groups according to length with splitRGbyLength, 
    # use output file in subsequent steps
    name=$(basename $bam .bam)
    atlas task=BAMDiagnostics bam=$bam logFile=${name}_BAMDiagnostics.log verbose 

    readGroups=($(cut -f1 ${name}_readLength.txt |sort |uniq \
        |grep -v RG |grep -v allReadGroups))

    for rg in ${readGroups[@]}
    do
        grep -P "^$rg\t" ${name}_readLength.txt | \
        sort -nrk 3 |head -n1|cut -f1,2 >>$name.lengths.txt
    done

    atlas task=splitRGbyLength bam=$bam readGroups=$name.lengths.txt \
        logFile=${name}_splitRGbyLength.log verbose

    #-------------------------------------------------------------------------#
    # Estimate the post-mortem damage patterns with estimatePMD
    echo ${readGroups[@]} > $name.tomerge.txt
    atlas task=estimatePMD bam=$bam fasta=$ref length=25 \
        logFile=${name}_estimatePMD.log verbose

    #-------------------------------------------------------------------------#
    # If the genome is male, use X recal to recalibrate the base quality scores.    
    if [ $sex eq "male" ]
    then
        atlas task=recal bam=$bam pmdFile=${name}_PMD_input_Empiric.txt \
        equalBaseFreq logFile=${name}_recal.log verbose 
    else
    # If the genome is female, or you have a group of samples that you want to
    # compare, use BQSR or recal based on invariant sites to recalibrate.
        atlas task=recal bam=$bam \
            pmdFile=${name}_PMD_input_Empiric.txt \
            window=$conserved logFile=${name}_recal.log equalBaseFreq verbose     
    fi
    #-----------------------------------------------------------------------------#
    # Use one of the variant discovery tools, one of the population genetics tools,
    #  or qualityTransformation. Take the base quality score recalibration and 
    # the PMD patterns into account by providing the corresponding files.
    # atlas task=estimateTheta bam=$bam pmdFile=${name}_PMD_input_Empiric.txt \
    # verbose recal=${name}_recalibrationEM.txt logFile=${name}_estimateTheta.log

    Rscript plot_theta_altas.R 1:22 $name

    atlas task=estimateTheta bam=$bam pmdFile=${name}_PMD_input_Empiric.txt \
        thetaGenomeWide minCoverage=2 recal=${name}_recalibrationEM.txt \
        logFile=${name}_estimateTheta.log verbose
}