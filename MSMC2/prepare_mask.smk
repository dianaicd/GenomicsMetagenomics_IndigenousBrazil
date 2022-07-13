configfile: "multiple_purposes.yaml"

#----------------------------------------------------------------#
# extract variables
inds = config["msmc2"]["inds"]
populations = config["msmc2"]["populations"]
combinations = config["msmc2"]["combinations"] if "combinations" in config["msmc2"] else {}

bootstrap_dir = config["msmc2"]["bootstrap_prefix"]
n_bootstraps = config["msmc2"]["n_bootstraps"]
chromosomes = [str(i) for i in range(1, 23)]
input_dir = config["msmc2"]["input_msmc2_dir"]
vcf = config["msmc2"]["vcf"]
#----------------------------------------------------------------#
# masks have to be prepared
# per ind
# per pop
# per comparison 
wildcard_constraints:
    ind = "|".join(inds),
    pop1 = "|".join(populations.keys()),
    pop2 = "|".join(populations.keys())


rule all:
    input:
        ind_multihetsep = [f"{input_dir}/{ind}.chr{chr}.txt" for ind in inds for chr in chromosomes],
        pop_multihetsep = [f"{input_dir}/pop.{p}.chr{chr}.txt" for p in populations.keys() for chr in chromosomes],
        comb_multihetsep = [f"{input_dir}/comb.{comb}.chr{chr}.txt" for comb in combinations for chr in chromosomes],
        boot_ind_multihetsep = [f"{bootstrap_dir}/{ind}_{rep}/bootstrap_multihetsep.chr{chr}.txt" 
                                for ind in inds for chr in chromosomes 
                                for rep in range(1, n_bootstraps+1)],
        boot_pop_multihetsep = [f"{bootstrap_dir}/pop.{p}_{rep}/bootstrap_multihetsep.chr{chr}.txt" 
                                for p in populations.keys() 
                                for chr in chromosomes 
                                for rep in range(1, n_bootstraps+1)],
        boot_comb_multihetsep = [f"{bootstrap_dir}/comb.{comb}_{rep}/bootstrap_multihetsep.chr{chr}.txt" 
                                for comb in combinations 
                                for chr in chromosomes 
                                for rep in range(1, n_bootstraps+1)]
        #multihetsep = expand("input_files_msmc2/chr{chr}.txt", chr = [str(i) for i in range(1,23)]),
        #bootstap = [f"{bootstrap_dir}_{i}/bootstrap_multihetsep.chr{chr}.txt" for chr in chromosomes for i in range(1, n_bootstraps+1)]


rule make_ind_mask:
    input:
        vcf=vcf
        #vcf="Final_phased/{ind}.chr{chr}.phased.vcf.gz"
    output:
        mask="msmc2_masks/mask_{ind}.chr{chr}.bed.gz"
    resources:
        runtime = 60 * 2,
        mem = 1024 * 2
    threads:
        4
    shell:
        """
        gunzip -c {input.vcf} | grep -v "#" | awk '{{print $1,$2-1,$2}}' |gzip > {output.mask} 
        """
#----------------------------------------------------------------#
# multihetsep
rule multi_ind:
    input:
        vcf=vcf,
        mask="msmc2_masks/mask_{ind}.chr{chr}.bed.gz"
    output:
        mask= "{input_dir}/{ind}.chr{chr}.txt"
    resources:
        runtime = 60 * 2,
        mem = 1024 * 2
    threads:
        4
    params:
        mappability_mask = "mappable_msmc-im/masks/um75-hs37d5.bed.gz"
    shell:
        """
        module load gcc/9.3.0 python/3.5.10
        ./generate_multihetsep.py \
            --chr {wildcards.chr} \
            --mask {input.mask} \
            --mask {params.mappability_mask} \
            {input.vcf} \
            > {output}
        """

rule multi_pop:
    input:
        all_vcf = lambda wildcards: [vcf.format(ind=ind, chr=wildcards.chr) for ind in populations[wildcards.p] ],
        mask=lambda wildcards:[f"msmc2_masks/mask_{ind}.chr{wildcards.chr}.bed.gz" for ind in populations[wildcards.p] ]
    output:
        mask= "{input_dir}/pop.{p}.chr{chr}.txt"
    resources:
        runtime = 60 * 12,
        mem = 1024 * 2
    threads:
        4
    params:
        mask_command=lambda wildcards:" ".join([f"--mask msmc2_masks/mask_{ind}.chr{wildcards.chr}.bed.gz" for ind in populations[wildcards.p]]),
        mappability_mask = "mappable_msmc-im/masks/um75-hs37d5.bed.gz"
    shell:
        """
        module load gcc/9.3.0 python/3.5.10
        ./generate_multihetsep.py \
            --chr {wildcards.chr} \
            {params.mask_command} \
            --mask {params.mappability_mask} \
            {input.all_vcf} \
            > {output}
        """

rule multi_combi:
    input:
        all_vcf = lambda wildcards: [vcf.format(ind=ind, chr=wildcards.chr) for ind in populations[wildcards.pop1] + populations[wildcards.pop2] ],
        mask1=lambda wildcards:[f"msmc2_masks/mask_{ind}.chr{wildcards.chr}.bed.gz" for ind in populations[wildcards.pop1] ],
        mask2=lambda wildcards:[f"msmc2_masks/mask_{ind}.chr{wildcards.chr}.bed.gz" for ind in populations[wildcards.pop2] ]
    output:
        mask= "{input_dir}/comb.{pop1}_{pop2}.chr{chr}.txt"
    resources:
        runtime = 60 * 12,
        mem = 1024 * 2
    threads:
        4
    params:
        mask_command=lambda wildcards:" ".join([f"--mask msmc2_masks/mask_{ind}.chr{wildcards.chr}.bed.gz" for ind in populations[wildcards.pop1] + populations[wildcards.pop2] ]),
        mappability_mask = "mappable_msmc-im/masks/um75-hs37d5.bed.gz"
    shell:
        """
        module load gcc/9.3.0 python/3.5.10
        ./generate_multihetsep.py \
            --chr {wildcards.chr} \
            {params.mask_command} \
            --mask {params.mappability_mask} \
            {input.all_vcf} \
            > {output}
        """


#----------------------------------------------------------------#

# rule make_all_mask:
#     input:
#         all_vcf=lambda wildcards:[f"Final_phased/{ind}.chr{wildcards.chr}.phased.vcf.gz" for ind in inds],
#         mask=lambda wildcards:[f"msmc2_masks/mask_{ind}.chr{wildcards.chr}.bed.gz" for ind in inds]
#     output:
#         "input_files_msmc2/chr{chr}.txt"
#     params:
#         mask_command=lambda wildcards:" ".join([f"--mask msmc2_masks/mask_{ind}.chr{wildcards.chr}.bed.gz" for ind in inds]),
#         #mappability_mask = lambda wildcards: f"--mask masks/chr{wildcards.chr}.bed.gz"
#         mappability_mask = "mappable_msmc-im/masks/um75-hs37d5.bed.gz"
#     resources:
#         mem = 1024 * 8,
#         runtime = 60 * 24
#     shell:
#         """
#         module load gcc/9.3.0 python/3.5.10
#         ./generate_multihetsep.py \
#             --chr {wildcards.chr} \
#             {params.mask_command} \
#             {params.mappability_mask} \
#             {input.all_vcf} \
#             > {output}
#         """


rule bootstrap:
    input:
        lambda wildcards: [f"{input_dir}/{wildcards.file}.chr{chr}.txt" for chr in chromosomes]
    output:
        expand("{bootstrap_dir}/{file}_{rep}/bootstrap_multihetsep.chr{chr}.txt" ,
         chr=chromosomes, rep=[rep for rep in range(1, n_bootstraps+1)], file = "{file}", bootstrap_dir = bootstrap_dir)
    params:
        dir = lambda wildcards: f"{bootstrap_dir}/{wildcards.file}",
        n_bootstraps = n_bootstraps,
        input =  lambda wildcards: f"{input_dir}/{wildcards.file}.chr*.txt"
    resources:
        mem = 1024 * 8,
        runtime = 60 * 24
    shell:
        """
        module load gcc/9.3.0 python/3.5.10
        ./multihetsep_bootstrap.py -n {params.n_bootstraps} \
          {params.dir}  {params.input} 
        """