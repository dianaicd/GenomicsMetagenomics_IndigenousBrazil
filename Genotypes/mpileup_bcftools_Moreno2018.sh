# Genotype calling for the panel

# Average depth >10X

# depth of coverage lower than 2 times the average depth 
# and greater than 1/3 the average depth
ref=/archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa
nThreads=12
chromosomes=($(seq 1 22))
panel=h

for chr in ${chromosomes[@]}
do

    bcftools mpileup -C 50 -b $bamfile -Q20 -q30 -a DP,SP \
    -f $ref --threads $nThreads -R ${panel}_${chr}_positions.txt | \
    bcftools call --threads $nThreads -v -c | \
    bcftools filter --threads $nThreads -e \
    '(DP > (2*AVG(DP))) || (DP < (AVG(DP)/3)) 
    || (SP < 1e-4) 
    || (GT="het" & 
        ( (COUNT(GT="R")/COUNT(GT="A")  
           | COUNT(GT="A")/COUNT(GT="R") < 0.2) ) )' \
     -Ob -o ${bamfile}_${chr}.bcf >out_${chr} 2>err_${chr}.txt &

done
# filter out:
# variants located within 5 bp of each other
# variants with a phred posterior probability lower than 30
# variants with a significant strand and/or end distance bias
# (p < 1e-4)
# heterozygous calls with an allelic balance lower than 0.2

