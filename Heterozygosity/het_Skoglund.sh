#!/bin/sh
workdir=~/scratch_monthly/Botocudos/Het/2019_05_13
cd $workdir
gitpath=~/data/Git/Botocudos-scripts
ln -s $gitpath/AlleleCounts/count_and_sample.py ./
ln -s $gitpath/AlleleCounts/vcf_sample_random.py ./
ln -s $gitpath/MDS/counts2dist.py ./
ln -s $gitpath/het_pi_Skoglund.py ./
ln -s $gitpath/Heterozygosity/remove_nondiallelic.r ./

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

# Compute heterozygosity as Pontus did before washing his hands
#-----------------------------------------------------------------------------#
# Select segregating sites in Sub-Saharan African populations
vcfpath=~/scratch_monthly/Simons/VCF
#cd $vcfpath

nThreads=6
# Countries to select from
countries=(BotswanaOrNamibia Central Congo Gambia Namibia Nigeria SierraLeone South)

if [ -e $vcfpath/Africa/selected_Africans.txt ]
then
    rm $vcfpath/Africa/selected_Africans.txt
fi

for COUNTRY in ${countries[@]}
do
    ls $vcfpath/Africa/$COUNTRY*gz >>selected_Africans.txt
done

while read line
do
    name=$(basename $line .vcf.gz)
    echo $name
    bcftools view -m2 -M2 -Ou  -v snps  $line \
        -e  'GT="hom" || (REF ~ "C" & ALT ~ "T") 
        || (REF ~ "G" & ALT ~ "A") || 
        (REF ~ "T" & ALT ~ "C") || (REF ~ "A" & ALT ~ "G")
        || (ALT ~ "C" & ALT ~ "T") ||(ALT ~ "G"& ALT ~ "A")' \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        > $vcfpath/Africa/sites_het_${name}.txt &

done < $vcfpath/Africa/selected_Africans.txt

# Identify non-diallelic sites
# This is bullshit
allvcf=($(ls $vcfpath/*/*vcf.gz))

if [ ! -d $vcfpath/NonDiallelic ] ; then mkdir $vcfpath/NonDiallelic ; fi 
for vcf in ${allvcf[@]}
do
    name=$( basename $vcf .vcf.gz)
    echo $name $vcf 
    bcftools view -m3 -Ou -v snps $vcf \
    | bcftools query -f '%CHROM\t%POS\n' > \
    $vcfpath/NonDiallelic/nondiallelic_${name}.txt &
done

cat $vcfpath/NonDiallelic/* |sort -nk1,2|uniq > $vcfpath/nondiallelic.txt
 
#-----------------------------------------------------------------------------#
# Select one African individual e.g., Yoruba
# Sample alleles for Botocudos according to heterozygous sites in Yoruba
cd $workdir

yoruba=$(grep -i yoruba $vcfpath/Africa/selected_Africans.txt |head -n1)
# file with paths to 22 BAM files (Botocudos)
bamfile=Botocudos
panel=$(basename $yoruba .vcf.gz)

# Remove non-diallelic sites
Rscript remove_nondiallelic.r $vcfpath/nondiallelic.txt \
    $vcfpath/Africa/sites_het_${panel}.txt

chromosomes=($(seq 1 22) X Y)

# Get sites and reference, alternative alleles
ln -s $vcfpath/Africa/sites_het_${panel}.txt ./

chromInPanel=($(cut -f1 sites_het_${panel}.txt |sort |uniq))
for chr in ${chromInPanel[@]}
do
    if [[ " ${chromosomes[*]} " == *$chr* ]]
    then
        : #echo $chr is desired
    else
        echo will remove $chr
        sed -i "/$chr/d" sites_het_${panel}.txt
    fi
done

cat sites_het_${panel}.txt |\
sed 's/\t/_/ ; s/A/0/g ; s/C/1/g ; s/G/2/g ; s/T/3/g' >$panel.refalt

# Run mpileup
for chr in ${chromosomes[@]}
do     
    cut -f 1,2 sites_het_${panel}.txt |grep -P "^$chr\t"|\
        awk '{print($1"\t"$2-1"\t"$2)}' >${panel}_${chr}_sites.bed

    samtools mpileup -r $chr -Bl ${panel}_${chr}_sites.bed \
    -b Boto/$bamfile \
     -a -o Boto/$bamfile_${chr}.mpileup -Q20 \
      > Boto/out_mpileup_${chr}.txt 2> Boto/err_mpileup_${chr}.txt & 
done
  { sleep 5; echo waking up after 5 seconds; } &
  { sleep 1; echo waking up after 1 second; } &
  wait
  echo "Samtools mpileup done"

# Count alleles and sample one base at random
if [ -e Boto/$bamfile.counts.sampled.txt ]
then rm Boto/$bamfile.counts.sampled.txt 
fi

for chr in ${chromosomes[@]}
do
    python count_and_sample.py Boto/$chr.mpileup \
        Boto/$chr.counts.gz Boto/$chr.sampled.gz Boto/$panel.refalt 

    gunzip -c Boto/${chr}.sampled.gz >> Boto/$bamfile.counts.sampled.txt
done

python3.5 het_pi_Skoglund.py -c Boto/$bamfile.counts.sampled.txt \
 -s Boto/$panel.refalt -o Boto/$bamfile

#-----------------------------------------------------------------------------#
# Sample random alleles for BAM files from the Americas shared by Víctor
fromVictor=~/Project/Americas/fromVictor
allpops=($( cat $fromVictor/ind_pop_region.txt \
        |grep Americas | cut -f2 |sed 's/ Americas//'|sort |uniq))
yoruba=$(grep -i yoruba $vcfpath/Africa/selected_Africans.txt |head -n1)

panel=$(basename $yoruba .vcf.gz)

chromosomes=($(seq 1 22) X Y)

# Get sites and reference, alternative alleles
ln -s $vcfpath/sites_het_${panel}.txt ./

chromInPanel=($(cut -f1 sites_het_${panel}.txt |sort |uniq))
for chr in ${chromInPanel[@]}
do
    if [[ " ${chromosomes[*]} " == *$chr* ]]
    then
        : #echo $chr is desired
    else
        echo will remove $chr
        sed -i "/$chr/d" sites_het_${panel}.txt
    fi
done

for chr in ${chromosomes[@]}
do     
    cut -f 1,2 sites_het_${panel}.txt |grep -P "^$chr\t"|\
            awk '{print($1"\t"$2-1"\t"$2)}' >${panel}_${chr}_sites.bed
done

cat sites_het_${panel}.txt |\
    sed 's/\t/_/ ; s/A/0/g ; s/C/1/g ; s/G/2/g ; s/T/3/g' >$panel.refalt

for pop in ${allpops[@]}
do
    echo $pop
    mkdir $pop 
    samples=($(grep $pop ~/Project/Americas/fromVictor/ind_pop_region.txt |cut -f1 |sed 's/$/.bam/'))

    if [ -e $pop/$pop ] ; then rm $pop/$pop ; fi
    for s in ${samples[@]} ; do ls $fromVictor/{Ancient,Modern}/$s >>$pop/$pop ; done

    bamfile=$pop/$pop

    # Run mpileup
    for chr in ${chromosomes[@]}
    do
        samtools mpileup -r $chr -Bl ${panel}_${chr}_sites.bed -b $bamfile \
        -a -o $pop/$bamfile_${chr}.mpileup -Q20 \
        > $pop/out_mpileup_${chr}.txt 2> $pop/err_mpileup_${chr}.txt & 
    done
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait
    echo "Samtools mpileup done"

    # Count alleles and sample one base at random
    if [ -e $bamfile.counts.sampled.txt ] ; then rm $bamfile.counts.sampled.txt ; fi     

    for chr in ${chromosomes[@]}
    do

        python count_and_sample.py --mpileup $pop/$chr.mpileup \
            --counts $pop/$chr.counts.gz --sampled $pop/$chr.sampled.gz --refalt $panel.refalt 

        gunzip -c $pop/${chr}.sampled.gz >> $bamfile.counts.sampled.txt
    done

    python3.5 het_pi_Skoglund.py --counts $pop/$pop.counts.sampled.txt  \
    --sites $panel.refalt --output $pop/$pop --autosomes_only

done

#-----------------------------------------------------------------------------#
# Sample random alleles for all populations in SGDP (VCF)
nThreads=1
maxJobs=150
allpops=($(cut -f4 ~/Project/Simons/Simons_sample_pop_region_country.txt \
            |sort |uniq |grep -v Population))

i=0

if [ ! -d Merged ] ; then mkdir Merged ; fi
chromosomes=($(seq 1 22) X Y)
for pop in ${allpops[@]}
do
    echo $pop
    samples=($(ls ~/scratch_monthly/Simons/VCF/*/*${pop}*vcf.gz))
    for chr in ${chromosomes[@]}
    do
        launchedJobs=$(echo "($i + 1) % $maxJobs" |bc)
        if [ $launchedJobs = 0 ]
        then
            echo "Waiting... $i"
            wait
        fi
        if [ ! -d $pop ] ; then mkdir $pop ; fi
        bcftools merge -0R ${panel}_${chr}_sites.bed --force-samples \
         -Ob $yoruba ${samples[@]} >Merged/${panel}_${pop}_${chr}.bcf &
    done
done ; wait 

# Look for non-diallelic sites
if [ ! -d NonDiallelic ] ; then mkdir NonDiallelic ; fi
rm NonDiallelic/*
for pop in ${allpops[@]}
do
    echo $pop
    for chr in ${chromosomes[@]}
    do
        bcftools index -f Merged/${panel}_${pop}_${chr}.bcf
        bcftools view -m3 -v snps Merged/${panel}_${pop}_${chr}.bcf |
        bcftools query -f '%CHROM\t%POS\n' >> NonDiallelic/${panel}_${chr}.txt
    done

done ; wait 

if [ -e $panel.refalt ] ; then rm ${panel}.refalt ; fi
for chr in ${chromosomes[@]}
do
    cut -f 1,2 sites_het_${panel}.txt |grep -P "^$chr\t" >${panel}_${chr}.txt
    Rscript remove_nondiallelic.r  NonDiallelic/${panel}_${chr}.txt ${panel}_${chr}.txt
    cat ${panel}_${chr}.txt >> $panel.refalt
done ; wait


maxJobs=5
i=0
for pop in ${allpops[@]}
do    
    echo $pop
    if [ -e $pop/$pop.counts ] ; then rm $pop/$pop.counts ; fi
    samples=($(ls ~/scratch_monthly/Simons/VCF/*/*${pop}*vcf.gz))
    lastField=$( echo ${#samples[@]} +1|bc)

    for chr in ${chromosomes[@]}
    do 
        echo $chr
        bcftools view -R ${panel}_${chr}.txt Merged/${panel}_${pop}_${chr}.bcf -Ou| \
         bcftools query -f '[%GT\t]\n' | \
            cut -f2-$lastField |sed 's|0/1|0\t1|g ; s|1/1|1\t1|g ; s|0/0|0\t0|g' \
            >$pop/$pop.${chr}.counts &
    done
    enough_jobs $maxJobs $i
    i=$(echo $i + 1|bc )
done ; wait


maxJobs=130
i=0
for pop in ${allpops[@]}
do    
    echo $pop
    if [ -e $pop/$pop.counts.txt ] ; then rm $pop/$pop.counts.txt ; fi
    for chr in ${chromosomes[@]}
    do 
        cat  $pop/$pop.${chr}.counts   >>$pop/$pop.counts.txt 
    done
    ncol=$(head -n1 $pop/$pop.counts.txt |wc -c)
    if [ $ncol -ge 4 ]
    then
        python3.5 het_pi_Skoglund.py --counts $pop/$pop.counts.txt \
            --sites $panel.refalt --output $pop/$pop --random_seed 123 \
            --called_genotypes &
        sleep 5s
    
        enough_jobs $maxJobs $i
        i=$(echo $i + 1|bc )
    fi
done ; wait

#-----------------------------------------------------------------------------#
# Sample random alleles for all populations in Team A,B (BAM)
source ~/data/Git/Botocudos-scripts/misc/source4merging.sh
bamlist=TeamAB
make_bamlist -s $bamlist -o $bamlist
sed -i '/0.0/d' $bamlist
mpileup --panel $panel --bamlist $bamlist \
    --chromosomes "$(echo ${chromosomes[@]})"

sample_mpileup --bamlist $bamlist --panel $panel \
    --chromosomes "$(echo ${chromosomes[@]})"

indices=/home/dcruzdva/archive/BAM/Reich_A_B/README_Bteam

# What about sampling with replacement?
while read line 
do
    id=$(echo $line | sed 's|.*/|| ; s/-dedup.*//')
    index=$(grep $id $indices |cut -f1)
    pop=$(grep $id $indices |cut -f2)
    echo "$id corresponds to line $index and pop $pop"
    i=$(echo "($index-1)*2 + 1" |bc -l )
    j=$(echo "$i+1" |bc) 
    echo $i,$j
    #cut -f$i,$j -d ' ' $bamlist.counts.sampled.txt > $pop.counts.txt
    if [ -e $pop.counts.txt ] ; then rm $pop.counts.txt ; fi 
    for chr in ${chromosomes[@]}
    do 
        less ${bamlist}_${chr}.counts.gz | cut -f $i,$j -d ' ' >> $pop.counts.txt
    done

    paste -d ' ' $pop.counts.txt $pop.counts.txt > tmp.counts ; mv tmp.counts $pop.counts.txt
    python count_and_sample.py \
        --counts $pop.counts.txt \
            --sampled ${pop}.sampled.gz \
            --refalt ${panel}.refalt

    python3.5 het_pi_Skoglund.py --counts $pop.sampled.gz  \
        --sites $panel.refalt --output $pop &
done < $bamlist 

python3.5 het_pi_Skoglund.py --counts $pop.sampled.gz  \
    --sites $panel.refalt --output $pop 

# And what about sampling without replacement?
while read line 
do
    id=$(echo $line | sed 's|.*/|| ; s/-dedup.*//')
    index=$(grep $id $indices |cut -f1)
    pop=$(grep $id $indices |cut -f2)
    echo "$id corresponds to line $index and pop $pop"
    i=$(echo "($index-1)*2 + 1" |bc -l )
    j=$(echo "$i+1" |bc) 
    echo $i,$j
    #cut -f$i,$j -d ' ' $bamlist.counts.sampled.txt > $pop.counts.txt
    if [ -e $pop.counts.txt ] ; then rm $pop.counts.txt ; fi 
    for chr in ${chromosomes[@]}
    do 
        less ${bamlist}_${chr}.counts.gz | cut -f $i,$j -d ' ' >> $pop.counts.txt
    done
    
    python3.5 het_pi_Skoglund.py --counts $pop.counts.txt  \
        --sites $panel.refalt --output $pop 

done < $bamlist 