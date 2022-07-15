genomes = [
    "A_Dai",
    "A_Dinka",
    "A_French",
    "A_Han",
    "A_Karitiana",
    "A_Mandenka",
    "A_Mbuti",
    "A_Papuan",
    "A_San",
    "A_Sardinian",
    "Aymara",
    "A_Yoruba",
    "B_Australian1",
    "B_Australian2",
    "B_Dai",
    "B_Dinka",
    "B_French",
    "B_Han",
    "B_Karitiana",
    "B_Mandenka",
    "B_Mbuti",
    "B_Mixe",
    "B_Papuan",
    "B_San",
    "B_Sardinian",
    "B_Yoruba",
    "Huichol",
    "Karitiana_HGDP_TeamAB",
    "MayaG",
    "MayaH",
    "MN0008_L3U_trim2",
    "Quechua",
    "SpCave",
    "SuruiA",
    "SuruiB"
]

chromosomes = [str(i) for i in range(1,23)]

rule all:
    input:
        variant = expand(
            "variant/chr{chr}/{sample}_chr{chr}.vcf.gz",
            chr = chromosomes,
            sample = genomes
        ),
        masks = expand(
             "called_masks/{sample}/{sample}_chr{chr}.bed",
             chr = chromosomes,
             sample = genomes
        )


rule bcf_to_vcf:
    input:
        "Filtered/chr{chr}/{sample}_maskStrict_repeatsUCSC_dp_GQ30_chr{chr}.bcf"
    output:
        temp("VCF/chr{chr}/{sample}_chr{chr}.vcf")
    resources:
        runtime = 60 * 4
    shell:
        """
        bcftools view -Ov -o {output} {input}
        """

rule subtract:
    input:
        bed = "genome_bed/chr{chr}.bed",
        vcf = "VCF/chr{chr}/{sample}_chr{chr}.vcf"
    output:
        bed = "called_masks/{sample}/{sample}_chr{chr}.bed"
    resources:
        mem = 1024 * 128
    shell:
        """
        module load gcc/9.3.0 bedtools2/2.29.2
        bedtools subtract -a {input.bed} -b {input.vcf} > {output.bed}
        """


rule extract_variants:
    input:
        "Filtered/chr{chr}/{sample}_maskStrict_repeatsUCSC_dp_GQ30_chr{chr}.bcf"
    output:
        "variant/chr{chr}/{sample}_chr{chr}.vcf.gz"
    resources:
        runtime = 60 * 4
    shell:
        """
        bcftools view -m2 -M2 -Oz -o {output} {input} 
        """

#bedtools subtract -a genome_bed/chr22.bed -b test.vcf > output.bed