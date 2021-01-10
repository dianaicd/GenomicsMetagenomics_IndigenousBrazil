configfile: "multiple_purposes.yaml"
include: "parse_resources.smk"
#include: "genotype_calling_all.smk"
#=============================================================================#
# Some defaults
default_roh = {"homozyg-window-kb ": 5000,
        "homozyg-window-snp": 50,
        "homozyg-window-het": 1,
        "homozyg-window-missing": 5,
        "homozyg-window-threshold": 0.05,
        "homozyg-snp": 100,
        "homozyg-kb": 1000,
        "homozyg-density": 50, 
        "homozyg-gap": 1000
        }

map_filter_name = {"window_kb":"homozyg-window-kb ",
        "window_snp": "homozyg-window-snp",
        "hets_window": "homozyg-window-het",
        "window_missing": "homozyg-window-missing",
        "window_threshold": "homozyg-window-threshold",
        "segment_snp": "homozyg-snp",
        "segment_kb":"homozyg-kb",
        "segment_density": "homozyg-density", 
        "homozyg_gap": "homozyg-gap"}

def build_filter(par_name):
    myFilter = ""
    if "filters" in config['roh_plink'].keys():
        filters = list(config["roh_plink"]["filters"].keys())
        if par_name in filters:
            myFilter = "--" + map_filter_name[par_name] + " " + str(config["roh_plink"]["filters"][par_name])
        # else:
            #myFilter = "--" + map_filter_name[par_name] + " " + str(default_roh[map_filter_name[par_name]])
            
    return(myFilter)

ref_genome = config["ref_genome"]

with open(ref_genome + ".fai", 'r') as index:
    chromosomes = [line.split()[0] for line in index.readlines()]
    chromosomes = [str(x) for x in range(1,23)]

bamlists = [l for l in list(config["roh_plink"]["bamlists"].keys())]

#=============================================================================#

rule all:
    input:
        roh = expand("ROH/{bamfile}.chr{chr}.phased_{rmTrans}_roh.hom",
                     bamfile = bamlists, chr = chromosomes, rmTrans = ["all", "rmTrans", "1240K"])


rule bcf_to_vcf:
    input:
        bcf = "bcf/{file}.bcf"
    output:
        vcf = temp("vcf/{file}_all.vcf.gz")
    log:
        'logs/bcf_to_vcf_{file}.txt'
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("bcf2vcf_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("bcf4vcf_time", attempt, 4)
    shell:
        """
        bcftools view -Oz -o {output.vcf} {input.bcf}
        """

rule rmTrans_vcf:
    input:
        bcf = "bcf/{file}.bcf"
    output:
        vcf = temp("vcf/{file}_rmTrans.vcf.gz")
    threads:
        4
    log:
        'logs/rmTrans_vcf_{file}.vcf'
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("vcf_rmTrans_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("vcf_rmTrans_time", attempt, 4)
    shell:
        """
        bcftools filter --exclude \
            "(REF=='C' & ALT='T') || (REF=='T' & ALT == 'C') || 
            (REF=='G' & ALT='A') || (REF=='A' & ALT == 'G')" \
            --threads {threads} \
            -Oz -o {output.vcf} {input.bcf}
        """


rule vcf_to_plink:
    input:
        vcf = "vcf/{file}.vcf.gz"
    output:
        bed = "bed/{file}.bed"
    log:
        'logs/vcf_to_plink_{file}.txt'
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("vcf2ed_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("vcf2bed_time", attempt, 4)
    params:
        out = "bed/{file}"
    shell:
        """
        plink --recode --geno 0.95 --make-bed --vcf {input.vcf} --out {params.out} --const-fid
        """

rule extract_1240K:
    input:
        bed = "bed/{file}.chr{chr}.phased_all.bed",
        to_extract = '1240K/to_extract_chr{chr}.txt'
    output:
        temp_bim = temp('bed/{file}.chr{chr}.phased_all_temp.bim'),
        bed = 'bed/{file}.chr{chr}.phased_1240K.bed'
    log:
        'logs/extract_1240K_{file}_{chr}.txt'
    resources:
    params:
        bed_in = 'bed/{file}.chr{chr}.phased_all',
        bed_out = 'bed/{file}.chr{chr}.phased_1240K',
        bim = 'bed/{file}.chr{chr}.phased_all.bim',
        map = 'bed/{file}.chr{chr}.phased_all.map'
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} {{print $1,$1"_"$4,$3,$4,$5,$6}}' {params.bim} > {output.temp_bim}
        cp {output.temp_bim} {params.bim}
        plink --make-bed --bfile {params.bed_in} --extract {input.to_extract} --out {params.bed_out}
        """
        
rule roh_plink:
    input:
        bed = "bed/{file}.bed"
    output:
        roh = "ROH/{file}_roh.hom"
    params:
        bed = "bed/{file}",
        roh = "ROH/{file}_roh",
        window_snp = build_filter("window_snp"),
        hets_window = build_filter("hets_window"),
        window_missing = build_filter("window_missing"),
        window_threshold = build_filter("window_threshold"),
        segment_density = build_filter("segment_density"),
        segment_snp = build_filter("segment_snp"),
        segment_kb = build_filter("segment_kb"),
        homozyg_gap = build_filter("homozyg_gap")
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("roh_plink_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("roh_plink_time", attempt, 4)
    log:
        'logs/roh_plink_{file}.txt'
    shell:
        """
        plink --homozyg \
        --bfile {params.bed} \
        --out {params.roh} \
        {params.window_snp} \
        {params.hets_window} \
        {params.window_missing} \
        {params.window_threshold} \
        {params.segment_density} \
        {params.segment_snp} \
        {params.segment_kb} \
        {params.homozyg_gap}
        """