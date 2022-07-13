# match variants to populations
# coordinates SuruiA_decoded.Summary.txt 
# which ones have data: SuruiA_decoded.All_posterior_probs.txt
localrules: merge_plink
# Rscript variants2bed.R SuruiA_decoded.Summary.txt SuruiA_decoded.All_posterior_probs.txt Australasian aus
ref = "/work/FAC/FBM/DBC/amalaspi/popgen/shared_ressources/reference_human/hs.build37.1/hs.build37.1.fa"
chromosomes = [str(i) for i in range(1, 23)]
rule all:
    input:
        ascertained = expand( "ascertained_SuruiA/{ind2}.tped", ind2 = ["Mixe", "French", "WCD02", "Han"]),
        target_ind =  "plink/SuruiA_subset.tped",
        match_rate = expand("matches/SuruiA_{ind2}.match", ind2 = ["Mixe", "French", "WCD02", "Han"])

rule variants_to_bed:
    input:
        summary = "{ind}_decoded.Summary.txt",
        all = "{ind}_decoded.All_posterior_probs.txt"
    output:
        vars_bed = "{ind}_introgressed.bed",
        vars_block = "{ind}_introgressed.txt"
    params:
        state = "Australasian",
        prefix = "{ind}_introgressed"
    shell:
        """
        Rscript variants2bed.R {input.summary} {input.all} {params.state} {params.prefix}
        """

∫
# for chr in $(seq 1 22) ; do 
# bcftools view --types snps -Ob -o \
# SuruiA_subset_chr$chr.bcf \
# -R ../aus.bed SuruiA_maskStrict_repeatsUCSC_dp_GQ30_chr$chr.bcf ; done


rule extract_introgressed:
    input:
        bcf = "bcf/{ind}_maskStrict_repeatsUCSC_dp_GQ30_chr{chr}.bcf",
        bed = "{ind}_introgressed.bed"
    output:
        bcf = "bcf_introgressed/{ind}_subset_chr{chr}.bcf"
    params:
    shell:
        """
        bcftools view --types snps \
        -Ob -o  {output.bcf} \
            -R {input.bed} {input.bcf}
        """


#for chr in $(seq 1 22) ; do plink --recode --bcf SuruiA_subset_chr$chr.bcf --out SuruiA_subset_chr$chr   ; cp SuruiA_subset_chr$chr.map old_chr$chr.map ; awk '{print $1,$1"_"$4,$3,$4}' old_chr$chr.map > SuruiA_subset_chr$chr.map ; done

rule recode_plink:
    input:
        bcf = "bcf_introgressed/{ind}_subset_chr{chr}.bcf"
    output:
        ped = "plink/{ind}_subset_chr{chr}.ped",
        map = "plink/{ind}_subset_chr{chr}.map"
    params:
        prefix = "plink/{ind}_subset_chr{chr}"
    shell:
        """
        mkdir -p plink
        plink --recode --bcf {input.bcf} --out {params.prefix} 
        cp {output.map} {output.map}_old ; awk '{{print $1,$1"_"$4,$3,$4}}' {output.map}_old > {output.map}
        """

# plink --merge-list to_merge --out SuruiA_subset --recode transpose
# ref=~/shared_ressources/reference_human/hs.build37.1/hs.build37.1.fa

rule merge_plink:
    input:
        ped = lambda wildcards: expand("plink/{ind}_subset_chr{chr}.ped", ind = wildcards.ind, chr = chromosomes),
        map = lambda wildcards: expand("plink/{ind}_subset_chr{chr}.map", ind = wildcards.ind, chr = chromosomes)
    output:
        "plink/{ind}_subset.tped",
    params:
        prefix = "plink/{ind}_subset"
    shell:
        """
        for f in {input.ped} ; do echo $f ; done  > peds
        for f in {input.map} ; do echo $f ; done  > maps
        paste peds maps > to_merge

        plink --merge-list to_merge --out {params.prefix} --recode transpose
        """

rule call_ascertained:
    input:
        bed = "{ind1}_introgressed.bed",
        ref = ref,
        bam = "/users/dcruzdav/popgen/Panels/fromVictor/Modern/{ind2}.bam"
    output:
        bcf = "ascertained_{ind1}/{ind2}.bcf",
        tped =  "ascertained_{ind1}/{ind2}.tped",
    params:
        prefix = "ascertained_{ind1}/{ind2}"
    resources:
        runtime = 4 * 60,
        mem = 2 * 1024
    shell:
        """
        bcftools mpileup -R {input.bed} -f {input.ref} \
            -Ob {input.bam} | bcftools call -c |bcftools view -V indels -Ob -o {output.bcf}
        mkdir -p ascertained 
        plink --recode transpose --bcf {output.bcf} --out {params.prefix}
        """

rule fasta_index:
    input:
        "{file}.fa"
    output:
        "{file}.fa.fai"
    shell:
        """
        samtools faidx {input}
        """


rule query_ancestral:
    input:
        bed = "{ind}_introgressed.bed",
        fasta = "homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{chr}.fa",
        index = "homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{chr}.fa.fai",
    output:
        bed = temp("anc_{ind}_chr_{chr}.bed"),
        alleles = "ancestral/{ind}/alleles_chr{chr}.txt"
    params:
    shell:
        """
        module load gcc/9.3.0 bedtools2/2.29.2
        newChr=$(cut -f1 {input.index}) 
        grep -P "^{wildcards.chr}\t" {input.bed} | sed "s/{wildcards.chr}/$newChr/" > {output.bed}
        bedtools getfasta -fi {input.fasta} -bed {output.bed} -fo {output.alleles} -tab
        """

rule merge_ancestral:
    input:
        expand("ancestral/{ind}/alleles_chr{chr}.txt", chr = chromosomes, ind = "{ind}")
    output:
        "ancestral/{ind}/alleles.txt"
    params:
    shell:
        """
        cat {input} > {output}
        """

rule calc_match_rate:
    input:
        tped_obs = "plink/{ind1}_subset.tped",
        vars_block = "{ind1}_introgressed.txt",
        tped_source =  "ascertained_{ind1}/{ind2}.tped",
        ancestral = "ancestral/{ind1}/alleles.txt"
    output:
        match_rate = "matches/{ind1}_{ind2}.match"
    params:
    shell:
        """
        Rscript calc_match_rate.R {input.tped_obs} {input.vars_block} {input.tped_source} {input.ancestral} {output.match_rate}
        """

# for chr in $(seq 1 22) ; do vcftools --gzvcf ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --freq --chr $chr --out ~/popgen/Botocudos/findAustralasian/2022_03_01/freqs_1000G/chr$chr --positions ~/popgen/Botocudos/findAustralasian/2022_03_01/for_vcftools ; done
        
# for chr in $(seq 1 22)
# do 
#     bcftools query -f '%CHROM\t%POS\n' bcf_introgressed/SuruiA_random_chr$chr.bcf
# done > random_for_vcftools

# module load gcc/9.3.0 vcftools/0.1.14
# for chr in $(seq 1 22) 
# do 
# vcftools --gzvcf ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
# --freq --chr $chr --out ~/popgen/Botocudos/findAustralasian/2022_03_01/freqs_1000G/random/chr$chr \
# --positions ~/popgen/Botocudos/findAustralasian/2022_03_01/random_for_vcftools 
# done