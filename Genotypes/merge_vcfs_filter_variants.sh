
samples=($(ls Filtered/*depth_filter.vcf.gz |cut -f2 -d/ |cut -f1 -d_)) ; echo "${samples[@]}"

cp Filtered/${samples[0]}_Nigeria_B_Yoruba-3_depth_filter.vcf out.vcf
bgzip out.vcf
tabix out.vcf.gz
# samples="${samples[@]/${samples[0]}}" ; echo "${samples[@]}"

for s in ${samples[@]}
do
    echo "Merging $s"
    bgzip Filtered/${s}_Nigeria_B_Yoruba-3_depth_filter.vcf
    tabix Filtered/${s}_Nigeria_B_Yoruba-3_depth_filter.vcf.gz

    bcftools merge \
    -Oz -o out.tmp.vcf.gz -m snps \
    out.vcf.gz Filtered/${s}_Nigeria_B_Yoruba-3_depth_filter.vcf.gz

    mv out.tmp.vcf.gz out.vcf.gz
    tabix out.vcf.gz
done


samples=($(ls Filtered/*depth_filter.vcf.gz |cut -f2 -d/ |cut -f1 -d_)) ; echo "${samples[@]}"

to_remove=(LS8 LS4 LS6 LS7)
for rem in ${to_remove[@]}
do
    samples=("${samples[@]/$rem}") ; echo "${samples[@]}"
done

bcftools view -s $(echo "${samples[@]}" |sed 's/ \+/,/g' ) -Ov -o out_diallelic.vcf -M2 out.vcf.gz 

geno=0.1
plink --recode \
    --make-bed --vcf out_diallelic.vcf \
    --out out_diallelic_geno${geno}_maf0.05 \
    --geno $geno --maf 0.05 --double-id

plink --homozyg --bfile out_diallelic_geno${geno}_maf0.05 --out out_diallelic_geno${geno}_maf0.05
