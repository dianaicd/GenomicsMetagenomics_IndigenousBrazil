



samples_age = {
    "KaritianaA": 0,
    "KaritianaB": 0,
    "MayaG": 0,
    "MayaH": 0,
    "MN0008": 200,
    "SuruiA": 0,
    "SuruiB": 0
}

all_samples = list( samples_age.keys() )
my_combinations = [f"{all_samples[i]}_{all_samples[j]}" for i in range(0, len(all_samples)) for j in range(0, len(all_samples))]
mutation_types = {"all":"SGDP_mutages/SGDP_half_ne_fixed", "transv": "SGDP_mutages/SGDP_half_ne_fixed_transv"}
rule all:
    input:
        bam_col = expand(
            "bam_output_{mut}/{comb}.coal", comb = my_combinations, mut = ["all", "transv"]
        )

rule precompute_colate_bam:
    input:
        target_bam = "bams/{sample}.bam",
    output:
        "bam_precomp_{mut}/{sample}.colate.in"
    params:
        ref_genome = "hs37",
        out = "bam_precomp_{mut}/{sample}",
        mut = lambda wildcards: mutation_types[wildcards.mut]
        #mut = "SGDP_mutages/SGDP_half_ne_fixed_transv"
    resources:
        runtime = 60 * 8
    shell:
        """
        ./Colate \
            --mode make_tmp \
            --mut {params.mut} \
            --target_bam {input.target_bam} \
            --ref_genome {params.ref_genome} \
            --chr chr.txt \
            -o {params.out}
        """

rule colate:
    input:
        precomp_target = "{type}_precomp_{mut}/{target}.colate.in",
        precomp_ref = "{type}_precomp_{mut}/{ref}.colate.in"
    output:
        "{type}_output_{mut}/{target}_{ref}.coal"
    params:
        #mut = "SGDP_mutages/SGDP_half_ne_fixed_transv",
        mut = lambda wildcards: mutation_types[wildcards.mut],
        bootstrap = 100,
        years_gen = 28,
        target_age = lambda wildcards: samples_age[wildcards.target],
        ref_age = lambda wildcards: samples_age[wildcards.ref],
        out = "{type}_output_{mut}/{target}_{ref}"
    resources:
        runtime = 60 * 8
    shell:
        """

        #mut="example_fixed" #name of .mut files obtained from step 1 (or downloaded)
        bins="3,7,0.2" #epochs in log10 years (format: start,end,stepsize)
        ./Colate \
            --mode mut \
            --mut {params.mut} \
            --target_tmp {input.precomp_target} \
            --reference_tmp {input.precomp_ref}\
            --bins ${{bins}} \
            --chr chr.txt \
            --num_bootstraps {params.bootstrap} \
            --target_age {params.target_age} \
            --reference_age {params.ref_age} \
            --years_per_gen {params.years_gen} \
            -o {params.out}
        """

#./Colate --mode mut --mut ${mut} --target_bam bams/MN0008_L3U.hg19_trim2.bam --reference_bam bams/KaritianaA.bam --ref_genome hs37   --bins ${bins}  --chr chr.txt  --target_age 200 --reference_age 28 -o example_output