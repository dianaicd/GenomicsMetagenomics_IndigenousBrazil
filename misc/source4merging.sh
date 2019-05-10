
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
    maxJobs=$1
    i=$2
    launchedJobs=$(echo "($i + 1) % $maxJobs" |bc)
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
            -p|--panel)    panel=$2; shift;;
            -b|--bamlist)  bamlist=$2; shift;;
            -c|--chromosomes) chromosomes=($(echo $2)) ; shift ;;
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
            -p|--panel)     panel=$2 ; shift ;;
            -t|--type)   type=$2 ; shift ;;
            -c|--chromosomes) chromosomes=($(echo $2)) ; shift ;;
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
        cut -f 1,2 ${panel}.refalt |grep -P "^$chr_"| sed 's/_/\t/'|\
        awk '{print($1"\t"$2-1"\t"$2)}' >${panel}_${chr}_sites.bed &

        grep -P "^$chr\_" ${panel}.refalt > ${panel}_${chr}.refalt &
    done
    wait
}
#-----------------------------------------------------------------------------#
# Count alleles and sample one base at random
sample_mpileup(){
    ped=""

    SHORTOPTS="p:b:Fc:"
    LONGOPTS="panel: bamlist: chromosomes: pedformat"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)    panel=$2; shift;;
            -b|--bamlist)  bamlist=$2; shift;;
            -c:--chromosomes) chromosomes=($(echo $2)); shift;;
            -F|--pedformat) ped="--ped"; shift;;
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
    SHORTOPTS="p:n:s:h:b:r"
    LONGOPTS="panel: name: selected: homozygous: bampath: rmdamage:"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)         panel=$2; shift;;
            -n|--name)          name=$2; shift;;
            -s|--selected)      selected=($(echo $2)); shift;;
            -h|--homozygous)    homo=$2; shift;;
            -b|--bampath)       bam_path=$2; shift;;
            -r|--rmdamage)      rmdamage=$2; shift;;
        esac 
        shift
    done

    if [ -e workflow_genolike.sh ]
    then
        ln -s ~/data/Git/Botocudos-scripts/GenoLike/workflow_genolike.sh ./
    fi

    ./workflow_genolike.sh $dir $panel $name "$(echo ${selected[@]})" $homo $bam_path $rmdamage
}
#-----------------------------------------------------------------------------#
# BED to tPED
bed2tped(){
    panel=$1
    plink --recode transpose --bfile $panel --out $panel 
}
#-----------------------------------------------------------------------------#
# BED to VCF
BED2VCF(){
    panel=$1
    panel=$(basename $panel .bed)
    plink --recode vcf --bfile $panel --out $panel
}
#-----------------------------------------------------------------------------#
# Call genotypes for data with genome coverage > 10x
genos_modern(){
    if [ ! -e mpileup_bcftools_Moreno2018.sh ] ; then
        ln -s ~/data/Git/Botocudos-scripts/Genotype/mpileup_bcftools_Moreno2018.sh ./
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
            -s|--samples) samples=$2; shift;;
            -o|--output) output=$2; shift;;
        esac
        shift
    done

    botoAll=~/archive/BAM/Botocudos/2018_10_26/link_final
    Posth=~/archive/BAM/Posth
    SGDP=~/scratch_monthly/Simons/*
    Yana=~/archive/BAM/Yana
    MalTa=~/archive/BAM/MalTa
    MaanasaAncient=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/ancient
    MaanasaNewWorld=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/newworld 
    MaanasaOldWorld=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/oldworld 
    MaanasaAll=~/archive/Panels/Raghavan2015/www.cbs.dtu.dk/suppl/NativeAmerican/data/alignments/*
    Lindo=~/archive/Panels/Lindo2018
    Scheib="pending"
    fromVictorAncient=~/Project/Americas/fromVictor/Ancient 
    fromVictorModern=~/Project/Americas/fromVictor/Modern 

    case "$samples" in 
        BotocudosAll)       ls $botoAll/*bam ;;     # Botocudos24, Botocudos22, BotocudosAll,
        Botocudos22)        ls $botoAll/*bam |grep -v MN01701 |grep -v MN1943 ;;       
        Posth)              ls $Posth/*bam  ;;      # Posth
        SGDP)               ls $SGDP/*bam ;;        # SGDP
        Yana)               ls $Yana/*bam ;;        # Yana
        MalTa)              ls $MalTa/*bam ;;       # MaanasaAncient
        MaanasaAncient)     ls $MaanasaAncient/*bam ;;# MaanasaNewWorld
        MaanasaNewWorld)    ls $MaanasaNewWorld/*bam ;; # MaanasaOldWorld
        MaanasaOldWorld)    ls $MaanasaOldWorld/*bam ;; # MaanasaAll
        MaanasaAll)         ls $MaanasaAll/*bam ;;      # MalTa
        Lindo)              ls $Lindo/*bam ;;           # Lindo
        Scheib)             ls $Scheib/*bam ;;          # Scheib
        fromVictorAncient)  ls $fromVictorAncient/*bam ;;# fromVictorModern
        fromVictorModern)   ls $fromVictorModern/*bam ;; # fromVictorAncient
        ALL)    # ALL
    esac > $output 
}
#-----------------------------------------------------------------------------#
# Merge BAM to BED
mergeBAM2BED(){
    SHORTOPTS="p:b:"
    LONGOPTS="panel: bamlist:"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)     panel=$2; shift;;
             -b|--bamlist)  bamlist=$2; shift;;
        esac 
        shift
    done

    pathPanel=$(dirname $panel)
    panel=$(basename $panel .bed)

    if [ ! -e sampled_ped.pl ] ; then
        ln -s ~/data/Git/Botocudos-scripts/MDS/sample_ped.pl ./
    fi

    ln -s $pathPanel/$panel.{bed,bim,fam} ./


    prepare_sites --panel $panel --type bed 
    chromosomes=($(cut -f1 $panel.refalt |cut -f1 -d_| sort -n|uniq)); echo ${chromosomes[@]}

    mpileup --panel $panel --bamlist $bamlist \
        --chromosomes "$(echo ${chromosomes[@]})"

    sample_mpileup --bamlist $bamlist --panel $panel --pedformat \
        --chromosomes "$(echo ${chromosomes[@]})"
    
    bed2tped $panel 
    paste $panel.tped $bamlist.counts.sampled.txt -d ' ' > $panel.$bamlist.tped
    cp $panel.tfam $panel.$bamlist.tfam 
    
    while read line
    do
        name=$(basename $line .bam)
        echo "$name $name 0 0 0 -9" >> $panel.$bamlist.tfam 
    done < $bamlist

    plink --recode --make-bed --tfile $panel.$bamlist --out $panel.$bamlist 

    # Sample a random allele for the whole panel, making it as if it were haploid;
    # this is useful to do an MDS
    plink --recode --bfile ${panel}.${bamlist} --out ${panel}.${bamlist}

    perl sample_ped.pl -ped ${panel}.${bamlist}.ped \
    -out ${panel}.${bamlist}.haploid.ped

    cp ${panel}.${bamlist}.map ${panel}.${bamlist}.haploid.map

    plink --distance square gz 'flat-missing' \
    --file ${panel}.${bamlist}.haploid \
    --out ${panel}.${bamlist}.haploid
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