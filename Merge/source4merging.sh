
gitpath=~/data/Git/Botocudos-scripts

# Why do we need different formats?
# MDS: PED/BED
# ngsAdmix: Genotype likelihoods
# admixtools: PED -> eigenstrat
# Heterozygosity: Counts

# Output types:
# VCF
# PED/BED
# Genotype likelihoods
#-----------------------------------------------------------------------------#
enough_jobs(){
    local maxJobs=$1
    local i=$2
    local launchedJobs=$(echo "($i + 1) % $maxJobs" |bc)
    if [ $launchedJobs = 0 ]
    then
        echo "Waiting... $i"
        wait
    fi
}
#-----------------------------------------------------------------------------#
mpileup(){
    SHORTOPTS="p:b:c:"
    LONGOPTS="panel: bamlist: chromosomes:"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)    local panel=$2; shift;;
            -b|--bamlist)  local bamlist=$2; shift;;
            -c|--chromosomes) local chromosomes=($(echo $2)) ; shift ;;
        esac 
        shift
    done

    for chr in ${chromosomes[@]}
    do     
        samtools mpileup -r $chr -Bl ${panel}_${chr}_sites.bed \
        -b $bamlist \
        -a -o ${bamlist}_${chr}.mpileup -Q20 \
        > out_mpileup_${chr}.txt 2> err_mpileup_${chr}.txt & 
    done
    wait
}
#-----------------------------------------------------------------------------#
prepare_sites(){
    SHORTOPTS="p:t:c:"
    LONGOPTS="panel: type: chromosomes:"
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "error in preparing sites: wtf arguments"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)         local panel=$2 ; shift ;;
            -t|--type)          local type=$2 ; shift ;;
            -c|--chromosomes)   local chromosomes=($(echo $2)) ; shift ;;
        esac 
        shift
    done

    # Get sites and reference, alternative alleles
    if [ $type == "vcf" ]
    then
        cut -f 1,2,4,5 $panel.vcf |grep -v '#' \
        |sed 's/\t/_/ ;s/A/0/g ; s/C/1/g; s/G/2/g ; s/T/3/g' >$panel.refalt
    elif [ $type == "bed" ]
    then 
        cut -f 1,4,5,6 $panel.bim |\
        sed 's/\t/_/ ;s/A/0/g ; s/C/1/g; s/G/2/g ; s/T/3/g' >$panel.refalt
    elif [ $type == "ped" ]
    then
        cut -f 1,4,5,6 $panel.bim |\
        sed 's/\t/_/ ;s/A/0/g ; s/C/1/g; s/G/2/g ; s/T/3/g' >$panel.refalt
    fi 
    
    chromosomes=($(cut -f1 $panel.refalt |cut -f1 -d_| sort -n|uniq)); echo ${chromosomes[@]}

    for chr in ${chromosomes[@]}
    do
        cut -f 1,2 ${panel}.refalt |grep -P "^$chr\_"| sed 's/_/\t/'|\
        awk '{print($1"\t"$2-1"\t"$2)}' >${panel}_${chr}_sites.bed &

        grep -P "^$chr\_" ${panel}.refalt > ${panel}_${chr}.refalt &
    done
    wait
}
#-----------------------------------------------------------------------------#
# Count alleles and sample one base at random
sample_mpileup(){
    local ped=""
    local allmutations=""
    SHORTOPTS="p:b:Fc:a"
    LONGOPTS="panel: bamlist: chromosomes: pedformat allmutations"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)         local panel=$2; shift;;
            -b|--bamlist)       local bamlist=$2; shift;;
            -c:--chromosomes)   local chromosomes=($(echo $2)); shift;;
            -F|--pedformat)     ped="--ped"; shift;;
            -a|--allmutations)  allmutations="--allmutations"; shift;;
        esac 
        shift
    done

    if [ ! -e count_and_sample.py ] ;
    then ln -s $gitpath/AlleleCounts/count_and_sample.py ./ ; fi  

    for chr in ${chromosomes[@]}
    do
        python count_and_sample.py --mpileup ${bamlist}_${chr}.mpileup \
            --counts ${bamlist}_${chr}.counts.gz \
            --sampled ${bamlist}_${chr}.sampled.gz \
            --refalt ${panel}_${chr}.refalt $ped $allmutations &
    done
    wait

    if [ -e $bamlist.counts.sampled.txt ];
     then rm $bamlist.counts.sampled.txt ; fi
    for chr in ${chromosomes[@]}
    do
        gunzip -c ${bamlist}_${chr}.sampled.gz >> $bamlist.counts.sampled.txt
    done
}
#-----------------------------------------------------------------------------#
# Use angsd to count alleles and sample one base at random
# the best strategy is to split calculations by individuals
# trim 5 bases by default
angsd_haploid(){
    local ped=""
    local trim=5
    local nThreads=2
    local rmtrans=1
    SHORTOPTS="p:b:Fc:"
    LONGOPTS="panel: bamlist: chromosomes: pedformat"
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)         local panel=$2; shift;;
            -b|--bamlist)       local bamlist=$2; shift;;
            -c:--chromosomes)   local chromosomes=($(echo $2)); shift;;
            -F|--pedformat)     ped="--ped"; shift;;
        esac 
        shift
    done

    # Prepare sites for angsd
    cut -f1,4-6 $panel.bim |sort -nk1,2 >$panel.sites 
    angsd sites index $panel.sites

    # Sample a read per site per individual
    nLines=$(wc -l $bamlist |cut -f1 -d ' ')

    if [ ! -d log ] ; then mkdir log ; fi 
    for i in $(seq 1 $nLines)
    do
        bamdir=$(dirname $(sed -n ${i}p $bamlist))
        ind=$(basename $(sed -n ${i}p $bamlist ) .bam)

        if [ ! -d $ind ] ; then mkdir $ind ; fi 

        angsd -i $bamdir/$ind.bam -doHaploCall 1 -sites $panel.sites \
            -minQ 20 -minMapQ 30 -trim $trim -doCounts 1 -out $ind/$ind \
            -howOften 5000 -nThreads $nThreads -maxMis 2 -doMajorMinor 3 \
            > log/angsd.$ind.out 2> log/angsd.$ind.err &
     
    done ; wait 

    # Merge files


    for chr in ${chromosomes[@]}
    do
        python count_and_sample.py --mpileup ${bamlist}_${chr}.mpileup \
            --counts ${bamlist}_${chr}.counts.gz \
            --sampled ${bamlist}_${chr}.sampled.gz \
            --refalt ${panel}_${chr}.refalt $ped &
    done
    wait

    if [ -e $bamlist.counts.sampled.txt ];
     then rm $bamlist.counts.sampled.txt ; fi
    for chr in ${chromosomes[@]}
    do
        gunzip -c ${bamlist}_${chr}.sampled.gz >> $bamlist.counts.sampled.txt
    done
}
#-----------------------------------------------------------------------------#
# Convert VCF to haploid counts
vcf_haploid(){
    if [ ! -e vcf_sample_random.py ]
    then 
        ln -s ~/data/Git/Botocudos-scripts/AlleleCounts/vcf_sample_random.py ./
    fi

    panel=$1

    python vcf_sample_random.py $panel.vcf $panel.counts $panel.counts.sampled.txt
}
#-----------------------------------------------------------------------------#
# Make genotype likelihoods
geno_like(){
    SHORTOPTS="p:h:b:r"
    LONGOPTS="panel: homozygous: bamlist: rmdamage:"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)         local panel=$2; shift;;
            -h|--homozygous)    local homo=$2; shift;;
            -b|--bamlist)       local bamlist=$2; shift;;
            -r|--rmdamage)      local rmdamage=$2; shift;;
        esac 
        shift
    done

    if [ ! -e workflow_genolike.sh ]
    then
        ln -s ~/data/Git/Botocudos-scripts/GenoLike/workflow_genolike.sh ./
    fi

    ./workflow_genolike.sh --panel $panel --homozygous $homo \
    --bamlist $bamlist --rmdamage $rmdamage
}
#-----------------------------------------------------------------------------#
# BED to tPED
bed2tped(){
    local panel=$1
    plink --recode transpose --bfile $panel --out $panel 
}
#-----------------------------------------------------------------------------#
# BED to VCF
BED2VCF(){
    local panel=$1
    panel=$(basename $panel .bed)
    plink --recode vcf --bfile $panel --out $panel
}
#-----------------------------------------------------------------------------#
# Call genotypes for data with genome coverage > 10x
genos_modern(){
    if [ ! -e mpileup_bcftools_Moreno2018.sh ] ; then
        ln -s ~/data/Git/Botocudos-scripts/Genotypes/mpileup_bcftools_Moreno2018.sh ./
    fi
    panel=$1
    bamfile=$2
    nThreads=$3

    ./mpileup_bcftools_Moreno2018.sh $panel $bamfile $nThreads
}
#-----------------------------------------------------------------------------#
make_bamlist(){
    SHORTOPTS="s:o:"
    LONGOPTS="samples: output:"
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "HEEEELP"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -s|--samples) local samples=$2; shift;;
            -o|--output) local output=$2; shift;;
        esac
        shift
    done

    local botoAll=~/Project/Botocudos/BAM
    local Posth=~/Project/Americas/Posth/BAM
    local SGDP=~/scratch_monthly/Simons/*
    local Yana=~/archive/BAM/Yana
    local Malaspinas=~/Project/Botocudos/Malaspinas2014/
    local MalTa=~/archive/BAM/MalTa
    local MaanasaAncient=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/ancient
    local MaanasaNewWorld=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/newworld 
    local MaanasaOldWorld=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/oldworld 
    local MaanasaAll=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/*
    local Lindo=~/archive/Panels/Lindo2018
    local Scheib=~/scratch_monthly/Botocudos/BAM/Scheib/final 
    local fromVictorAncient=~/Project/Americas/fromVictor/Ancient 
    local fromVictorModern=~/Project/Americas/fromVictor/Modern 
    local TeamAB=~/archive/BAM/Reich_A_B

    case "$samples" in 
        BotocudosAll)       ls $botoAll/*bam ;;     # Botocudos24, Botocudos22, BotocudosAll,
        Botocudos22)        ls $botoAll/*bam |grep -v MN01701 |grep -v MN1943 ;;       
        Posth)              ls $Posth/*bam  ;;      # Posth
        SGDP)               ls $SGDP/*bam ;;        # SGDP
        Yana)               ls $Yana/*bam ;;        # Yana
        Malaspinas2014)     ls $Malaspinas/*bam ;;  # Malaspinas2014: Polynesians
        MalTa)              ls $MalTa/*bam ;;       # MaanasaAncient
        MaanasaAncient)     ls $MaanasaAncient/*bam ;;# MaanasaNewWorld
        MaanasaNewWorld)    ls $MaanasaNewWorld/*bam ;; # MaanasaOldWorld
        MaanasaOldWorld)    ls $MaanasaOldWorld/*bam ;; # MaanasaAll
        MaanasaAll)         ls $MaanasaAll/*bam ;;      # MalTa
        Lindo)              ls $Lindo/*bam ;;           # Lindo
        Scheib)             ls $Scheib/*bam ;;          # Scheib
        fromVictorAncient)  ls $fromVictorAncient/*bam ;;# fromVictorModern
        fromVictorModern)   ls $fromVictorModern/*bam ;; # fromVictorAncient
        TeamAB)             ls $TeamAB/*bam ;; # Team A, Team B
        ALL)    # ALL
    esac > $output 
}
#-----------------------------------------------------------------------------#
# Merge BAM to BED
mergeBAM2BED(){
    SHORTOPTS="p:b:a"
    LONGOPTS="panel: bamlist: allmutations"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    local allmutations=""
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)     local panel=$2; shift;;
            -b|--bamlist)  local bamlist=$2; shift;;
            -a|--allmutations) local allmutations="--allmutations"; shift;;
        esac 
        shift
    done

    local pathPanel=$(dirname $panel)
    local panel=$(basename $panel .bed)

    if [ ! -e sampled_ped.pl ] ; then
        ln -s ~/data/Git/Botocudos-scripts/MDS/sample_ped.pl ./
    fi

    ln -s $pathPanel/$panel.{bed,bim,fam} ./


    prepare_sites --panel $panel --type bed 
    local chromosomes=($(cut -f1 $panel.refalt |cut -f1 -d_| sort -n|uniq)); echo ${chromosomes[@]}

    mpileup --panel $panel --bamlist $bamlist \
        --chromosomes "$(echo ${chromosomes[@]})"

    sample_mpileup --bamlist $bamlist --panel $panel --pedformat \
        --chromosomes "$(echo ${chromosomes[@]})" $allmutations
    
    bed2tped $panel 
    paste $panel.tped $bamlist.counts.sampled.txt -d ' ' > $panel.$bamlist.tped
    cp $panel.tfam $panel.$bamlist.tfam 
    
    while read line
    do
        name=$(basename $line .bam)
        echo "$name $name 0 0 0 -9" >> $panel.$bamlist.tfam 
    done < $bamlist

    plink --recode --make-bed --tfile $panel.$bamlist --out $panel.$bamlist 
    awk 'BEGIN{DOF="\t"} {print $1,$2,$3,$4,$5,1}' \
            $panel.$bamlist.fam >tmp.fam ; mv tmp.fam $panel.$bamlist.fam 
    # Sample a random allele for the whole panel, making it as if it were haploid;
    # this is useful to do an MDS
    plink --recode --bfile ${panel}.${bamlist} --out ${panel}.${bamlist}

    perl sample_ped.pl -ped ${panel}.${bamlist}.ped \
    -out ${panel}.${bamlist}.haploid.ped

    cp ${panel}.${bamlist}.map ${panel}.${bamlist}.haploid.map

    for mind in 05 $(seq 10 5 95)
    do
        plink --distance square gz 'flat-missing' \
            --file ${panel}.${bamlist}.haploid \
            --mind 0.$mind \
            --out ${panel}.${bamlist}.mind0.$mind.haploid &
    done ;wait
    # Recode for Fred's format
    plink --recode 12 --file ${panel}.${bamlist}.haploid \
        --out ${panel}.${bamlist}.haploid.fred
    
    plink --recode transpose --file ${panel}.${bamlist}.haploid.fred \
        --out ${panel}.${bamlist}.haploid.fred


    local mind=0.95
    plink --tfile ${panel}.${bamlist}.haploid.fred --make-bed \
    --mind $mind --out ${panel}.${bamlist}.haploid.fred.mind${mind}

    bed2tped $panel.${bamlist}.haploid.fred.mind${mind}

    lastField=$(head -n1 ${panel}.${bamlist}.haploid.fred.mind${mind}.tped|wc -w)
    columns=$(echo $(seq 5 2 $lastField) | sed 's/ /,/g')

    cut -f $columns -d ' ' ${panel}.${bamlist}.haploid.fred.mind${mind}.tped \
     >${panel}.${bamlist}.haploid.fred.mind${mind}

     plink --distance square gz 'flat-missing' \
         --bfile ${panel}.${bamlist}.haploid.fred.mind${mind} \
         --out ${panel}.${bamlist}.haploid.fred.mind${mind}

}

bam2plink(){

    if [ ! -e ${panel}_regions ]
    then
        cat $panel.bim |awk '{print $1"\t"$4-1"\t"$4"\t"$5"\t"$6}' \
        >${panel}_regions
    fi

    # Sample a random allele from the bam files
    while read line
        do
        echo $line
        ind=$(basename $line .bam)
        bam2plink.py bamfile=$line plinkpref=$panel \
        MinBQ=20 indname=$ind.sampled popname=$ind doCns=F \
        trim=5 MinMQ=30 >out_${ind}_bam2plink.txt \
        2> err_${ind}_bam2plink.txt &
    done < $bamlist

    wait
    echo "An allele was sampled for each file in $bamlist"

    # Merge plink binary files to make an extended panel
    extendedPanel=$panel

    while read line 
    do
        ind=$(basename $line .bam)
        echo "-------------------- Merging $ind to $panel ---------------------"

        plink --bfile $extendedPanel --bmerge ${ind}.sampled_${ind}_${panel} \
        --make-bed --out ${panel}_${ind} --allow-no-sex
        rm $extendedPanel.*
        extendedPanel=${panel}_${ind}

    done < $bamlist

    for ext in bed bim fam log nosex
    do
        mv $extendedPanel.$ext ${panel}_${bamlist}.$ext
    done
}

#-----------------------------------------------------------------------------#
# Make par file for convertf
make_par(){
    SHORTOPTS="p:c:"
    LONGOPTS="panel: conversion:"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)         local panel=$2; shift;;
            -c|--conversion)    local conversion=$2; shift;;
        esac 
        shift
    done
    
    if [ -e ${panel}_${conversion}.par ] ; then rm ${panel}_${conversion}.par ; fi

    if [ $conversion == "ped2eigenstrat" ] ; then
        argument=(genotype snpname indivname outputformat genotypeoutname \
                    snpoutname indivoutname familynames)

        options=($panel.ped $panel.bim $panel.ped EIGENSTRAT $panel.eigenstratgeno \
                $panel.snp $panel.ind YES)
        
        for i in $(seq 0 7)
        do 
          echo "${argument[$i]}: ${options[$i]}" >> ${panel}_${conversion}.par
        done

    elif [ $conversion == "eigenstrat2ped" ] ; then 
        argument=(genotype snpname indivname outputformat genotypeoutname \
                    snpoutname indivoutname)

        options=($panel.eigenstratgeno $panel.snp $panel.ind PED $panel.ped \
                $panel.pedsnp $panel.pedind)

        for i in $(seq 0 7)
        do 
          echo "${argument[$i]}: ${options[$i]}" >> ${panel}_${conversion}.par
        done
    fi 
}