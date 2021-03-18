chr=1
yoruba_prefix=Yoruba/B_Yoruba_maskStrict_repeatsUCSC_dp_GQ30_chr
guarani=guarani_tupiniquim_masked.chr
out=merged/yoruba_guarani_tupiniquim_masked.chr

for chr in $(seq 1 22)
do 

    echo "----------------------- chr $chr --------------------"
#    plink --bfile guarani_tupiniquim_masked/guarani_tupiniquim_masked \
#        --chr $chr --out $guarani$chr --snps-only 'just-acgt' --make-bed

#    awk 'BEGIN {OFS="\t"} {print $1,$4-1,$4}' \
#        $guarani$chr.bim  > positions.chr$chr.bed

#    tabix ${yoruba_prefix}$chr.bcf
#    bcftools view -r $chr -R positions.chr$chr.bed \
#        -Ov -e 'INDEL=1' -o Yoruba.chr$chr.vcf \
#        ${yoruba_prefix}$chr.bcf

    plink --make-bed --vcf Yoruba.chr$chr.vcf --out Yoruba.chr$chr

    plink --set-missing-var-ids @:# \
        --bfile Yoruba.chr$chr --out Yoruba.chr$chr --make-bed

    awk '{print $1":"$4,$2}' $guarani$chr.bim > update.chr$chr.txt
    plink --update-name update.chr$chr.txt \
        --bfile Yoruba.chr$chr --out Yoruba.chr$chr --make-bed

    plink --bfile $guarani$chr \
        --bmerge Yoruba.chr$chr --make-bed \
        --out $out$chr

    if [ -e $out$chr-merge.missnp ]
    then
        echo "============= removing triallelic for chr $chr ==========="
        for file in  $guarani$chr Yoruba.chr$chr
        do
        plink --exclude $out$chr-merge.missnp \
            --bfile $file --make-bed --out $file.diallelic
        done

        plink --bfile $guarani$chr.diallelic \
            --bmerge Yoruba.chr$chr.diallelic --make-bed \
            --out $out$chr

        cat $out$chr-merge.missnp >> triallelic.txt
    fi
done

plink --merge-list to_merge.txt --make-bed --out ${out}All