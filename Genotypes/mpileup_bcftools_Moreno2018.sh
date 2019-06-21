# Genotype calling for the panel

# Average depth >10X

# depth of coverage lower than 2 times the average depth 
# and greater than 1/3 the average depth
ref=/archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa

bamfile=$1
panel=$(basename $2)
nThreads=$3

chromosomes=($(cut -f 1 $panel.bim |sort |uniq ))

if [ ! -d log ] ; then mkdir log ; fi
if [ ! -d RAW ] ; then mkdir RAW ; fi
for chr in ${chromosomes[@]}
do
    grep -P "^$chr\t" ${panel}.bim |cut -f1,4 > ${panel}_${chr}_positions.txt

        bcftools mpileup -C 50 -b $bamfile -q30 -Q20 -a FMT/DP,SP \
        -f $ref --threads $nThreads -r $chr \
        -R ${panel}_${chr}_positions.txt -Ou | \
        bcftools annotate -c RPB | \
        bcftools call --threads $nThreads -c -V indels \
        -Ob -o RAW/${bamfile}_${chr}.bcf \
        >log/out.raw.${bamfile}.txt 2> log/err.raw.${bamfile}.txt &
done ; wait

if [ ! -d Filter1 ] ; then mkdir Filter1 ; fi   

if [ -e ${bamfile}_toconcat.txt ] ; then rm ${bamfile}_toconcat.txt ; fi
for chr in ${chromosomes[@]}
do
    bcftools filter --threads $nThreads -e \
        "(SP > 40) 

        || (GT='het' && 
                ( (DP4[0]+DP4[1])/SUM(DP4) <= 0.2 ||
                  (DP4[2]+DP4[3])/SUM(DP4) <= 0.2   )
            )
        || (REF='C' & ALT='T') || (REF=='T' & ALT == 'C')
        || (REF='G' & ALT='A') || (REF=='A' & ALT == 'G')
        || N_ALT > 1 || QUAL <=30 || VDB == 0 ||RPB==0" \
        -o Filter1/${bamfile}_${chr}_f1.bcf RAW/${bamfile}_${chr}.bcf &

    echo "Filter1/${bamfile}_${chr}_f1.bcf" >> ${bamfile}_toconcat.txt
done; wait

bcftools concat -f ${bamfile}_toconcat.txt -Ob -o ${bamfile}_raw.bcf --threads $nThreads

sample=$(bcftools query -f '[%SAMPLE]\n' ${bamfile}_raw.bcf \
            |head -n1)
        
avgdp=$(bcftools stats -s $sample ${bamfile}_raw.bcf \
     |grep -P "^PSC" |cut -f 10); echo "average depth: $avgdp "

mindp=$avgdp/3
maxdp=2*$avgdp
mindp=10
maxdp=30
bcftools filter --threads $nThreads -e \
    "(SUM(DP4)< $mindp| SUM(DP4) > $maxdp )" ${bamfile}_raw.bcf\
    -Ob -o ${bamfile}.bcf 
        
bcftools stats -s $sample ${bamfile}.bcf  |less -S

#-----------------------------------------------------------------------------#
bcftools view -Ov -o ${bamfile}.vcf ${bamfile}.bcf
plink --recode --vcf ${bamfile}.vcf --out ${bamfile} --make-bed --set-missing-var-ids @:#
plink --bfile $panel --bmerge $bamfile --out ${panel}_${bamfile}

# plink --recode vcf --real-ref-alleles --bfile h --keep ls.txt --out lagoa
# filter out:
# variants located within 5 bp of each other
# variants with a phred posterior probability lower than 30 / SHOULD BE FIELD QUAL < 30
# variants with a significant strand and/or end distance bias
# (p < 1e-4)
# heterozygous calls with an allelic balance lower than 0.2 (accept a/sum(A,a) >0.2)

# bcftools stats lagoa.vcf -s Sumidouro5.variant_Sumidouro5.variant |less -S
# Victorvcf:  4,941,388; het: 193,682 ; nNonRefHom: 75,816 ;    nRefHom: 4,671,890
# bcf:  5,003,052; het: 197,570 ;       nNonRefHom: 112,780 ; nRefHom: 4,670,326

#chr18
#bcf: 147332cd 

# Filter from 10x to --x
# accept >=10, <=floor(2*$avgdp)