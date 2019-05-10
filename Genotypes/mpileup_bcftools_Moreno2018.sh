# Genotype calling for the panel

# Average depth >10X

# depth of coverage lower than 2 times the average depth 
# and greater than 1/3 the average depth
ref=/archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa

bamfile=$1
panel=$(basename $2)
nThreads=$3

chromosomes=($(cut -f 1 $panel.bim |sort |uniq ))

for chr in ${chromosomes[@]}
do
    grep -P "^$chr\t" ${panel}.bim |cut -f1,4 > ${panel}_${chr}_positions.txt

        bcftools mpileup -C 50 -b $bamfile -Q20 -q30 -a DP,SP \
        -f $ref --threads $nThreads -R ${panel}_${chr}_positions.txt | \
        bcftools call --threads $nThreads -c -V indels | \
        bcftools filter --threads $nThreads -e \
        '( DP > (2*AVG(DP)) ) || (DP < (AVG(DP)/3) ) 
        || (SP > 40 ) 
        || (GT="het" & 
            (COUNT(GT="A")/COUNT(GT="R") < 0.4) ) 
        || (REF="C" & ALT="T") || (REF=="T" & ALT == "C")
        || (REF="G" & ALT="A") || (REF=="A" & ALT == "G")' \
        -Ob -o ${bamfile}_${chr}.bcf >out_${chr} 2>err_${chr}.txt &

    echo "${bamfile}_${chr}.bcf" >> ${bamfile}_toconcat.txt
done

bcftools concat -f ${bamfile}_toconcat.txt -Ob -o ${bamfile}.bcf --threads $nThreads
bcftools view -Ov -o ${bamfile}.vcf ${bamfile}.bcf
plink --recode --vcf ${bamfile}.vcf --out ${bamfile} --make-bed --set-missing-var-ids @:#
plink --bfile $panel --bmerge $bamfile --out ${panel}_${bamfile}
# filter out:
# variants located within 5 bp of each other
# variants with a phred posterior probability lower than 30
# variants with a significant strand and/or end distance bias
# (p < 1e-4)
# heterozygous calls with an allelic balance lower than 0.2 (accept a/sum(A,a) >0.2)

