# merge VCFs

configfile: "config_merge_vcfs.yaml"
chromosomes = [str(i) for i in range(1, 23)]
inds_pop=config["individuals_to_merge"]
populations = list(inds_pop.keys())
groups = config["populations_to_merge"].keys()
groups_populations = config["populations_to_merge"]
print(f"groups_populations: f{groups_populations}")
ind_vcf={ind:vcf for population in populations for ind,vcf in inds_pop[population].items()}
print(f"ind_vcf: f{ind_vcf}")

wildcard_constraints:
    chr = "|".join(chromosomes)

rule all:
    input:
        grouped_vcf = [f"groups/{group}.vcf.{ext}" for group in groups for ext in ["gz", "gz.tbi"]],
        poplist = [f"poplists/{population}.txt" for population in populations]

rule reheader:
    input:
        bcf=lambda wildcards: ind_vcf[wildcards.ind]
        #"/users/dcruzdav/popgen/Botocudos/Genotypes/2020_11_13/Filtered/chr{chr}/{ind}_maskStrict_repeatsUCSC_dp_GQ30_chr{chr}.bcf"
    output:
        header=temp("reheaded/{ind}.chr{chr}.header"),
        bcf=temp("reheaded/{ind}.chr{chr}.bcf")
    resources:
        runtime = 60 * 4,
        mem = 1024 * 2
    shell:
        """
        echo {wildcards.ind} > {output.header}
        bcftools  reheader -s {output.header} -o {output.bcf} {input.bcf}
        """
    
rule index_bcf:
    input:
        "{file}.bcf"
    output:
        "{file}.bcf.csi"
    resources:
        runtime = 60 * 4,
        mem = 1024 * 2
    shell:
        """
        tabix -f {input}
        """

rule index_Vcf:
    input:
        "{file}.vcf.gz"
    output:
        "{file}.vcf.gz.tbi"
    resources:
        runtime = 60 * 4,
        mem = 1024 * 2
    shell:
        """
        tabix -f {input}
        """

rule make_pop_list:
    output:
        "poplists/{population}.txt"
    params:
        individuals = lambda wildcards: list(inds_pop[wildcards.population].keys())
    shell: 
        """
        for ind in {params.individuals}
        do
            echo "$ind\t{wildcards.population}"
        done > {output}
        """

rule fix_ancestral_allele:
    input:
        bcf="{file}.bcf",
        chimp_fa=config["ancestral"]
    output:
        temp("{file}.anc.bcf")
    resources:
        runtime = 60 * 4,
        mem = 1024 * 2
    shell:
        '''
        AC=$(echo "$(bcftools query -l {input.bcf} | wc -l) * 2\"  |bc)
        bcftools annotate -Ob -x INFO,^FORMAT/GT {input.bcf} | \
            bcftools norm -Ob -c s -f {input.chimp_fa} - | \
            bcftools view -Ob -e REF=\\"N\\" -a  -M2 -o {output}
            
        '''

# https://www.biostars.org/p/61267/
# https://github.com/jessstapley/Set_ancestral_allele_vcf
# https://github.com/CMPG/originsEarlyFarmers/blob/main/dataFilteringandAssembling

rule subset_neutral:
    input:
        bcf="{file}.anc.bcf",
        index="{file}.anc.bcf.csi",
        bed="neutral.bed"
    output:
        bcf="{file}.neutral.bcf",
    resources:
    shell:
        """
        bcftools view -R {input.bed} -Ob -o {output.bcf} {input.bcf}
        """

# Merge all individuals of a single population on the same chromosome
rule merge_same_chr:
    input:
        bcf=lambda wildcards: [
            f"reheaded/{ind}.chr{wildcards.chr}.neutral.bcf" 
            for ind in inds_pop[wildcards.population]
            ],
        bcf_csi=lambda wildcards: [
            f"reheaded/{ind}.chr{wildcards.chr}.neutral.bcf.csi" 
            for ind in inds_pop[wildcards.population]
            ]
    resources:
        runtime = 60 * 4,
        mem = 1024 * 2
    params:
        n_inds = lambda wildcards: len(inds_pop[wildcards.population])
    output:
        temp("populations/{population}/chr{chr}.bcf")
    shell:
        """
        if [ {params.n_inds} -eq 1 ]
        then
            cp {input.bcf} {output}
        else
            bcftools merge --merge snps -Ob -o {output} {input.bcf} 
        fi
        """

# Merge all chromosomes of a single population
rule merge_chrs_pop:
    input:
        bcf = lambda wildcards: [
            f"populations/{wildcards.population}/chr{chr}.bcf"
            for chr in chromosomes
        ],
        bcf_csi = lambda wildcards: [
            f"populations/{wildcards.population}/chr{chr}.bcf.csi"
            for chr in chromosomes
        ]

    output:
        vcf = "populations/{population}.vcf.gz"
    resources:
        runtime = 60 * 4,
        mem = 1024 * 2
    shell:
        """
        bcftools concat -Oz -o {output.vcf} {input.bcf}
        """

# # # Merge populations
rule merge_pops:
    input:
        vcf=lambda wildcards:[
            f"populations/{population}.vcf.gz"
            for population in groups_populations[wildcards.group]
        ],
        vcf_tbi=lambda wildcards:[
            f"populations/{population}.vcf.gz.tbi"
            for population in groups_populations[wildcards.group]
        ]
    resources:
        runtime = 60 * 4,
        mem = 1024 * 2
    params: 
        nChr = lambda wildcards: sum([2 for population in groups_populations[wildcards.group] for ind in inds_pop[population]])
    output:
        vcf="groups/{group}.vcf.gz"
    shell:
        """

        bcftools merge --merge snps -Ob  {input.vcf} | bcftools view -M2 -i 'AN={params.nChr}' -Oz -o {output}
        """