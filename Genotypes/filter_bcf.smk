# filter bcfs

configfile: "genotype_calling.yaml"
include: "parse_resources.smk"

import pandas as pd
from fractions import Fraction
#-----------------------------------------------------------------------------#
# note that the filters will be applied in the specified order 
filters = list(config["geno_calls"]["filters"].keys())
bamlists = list(config["geno_calls"]["individuals"].keys())
chromosomes = [str(i) for i in range(1, 23)]
#-----------------------------------------------------------------------------#
def get_filter_value(f, entry, default=""):
    if f in filters:
        value = config["geno_calls"]["filters"][f][entry]
    else:
        value = default
    
    return value

mask_label = get_filter_value("mask", "label")
min_dp = str(get_filter_value("depth", "min"))
max_dp = str(get_filter_value("depth", "max"))
gq = str(get_filter_value("gq", "min"))
repeats_label = get_filter_value("repeats", "label")

sufixes = {
    "mask": "mask{mask_label}".format(mask_label=mask_label),
    "depth": "dp",
    "gq": "GQ{gq}".format(gq=gq), 
    "repeats": "repeats{repeats_label}".format(repeats_label=repeats_label)
    }

filters_order = config["geno_calls"]["filters"]["order"]

def build_sufix(filters, config):
    sufix_file = "_".join([sufixes[f] for f in filters_order])

    return sufix_file

def get_input_file(current_filter, basename, chr):

    filters_order = config["geno_calls"]["filters"]["order"]
    index_filter = filters_order.index(current_filter)

    if index_filter:
        sufix = "" #"_" + sufixes[filters_order[index_filter - 1]]
        prefix = f"Filtered/chr{chr}/"
    else:
        sufix = ""
        prefix = "Raw/"
    
    bcf_name = prefix + basename + sufix + f"_chr{chr}.bcf"
    
    return bcf_name

def get_DoC(genomecov_path, chromosomes = [str(i) for i in range(1,23)]):
    genomecov = pd.read_csv(
        genomecov_path, 
        delimiter="\t", 
        header=None, 
        names=["chr", "depth", "counts", "length", "fraction"],
        dtype = {
            "chr":"str", 
            "depth":"int64", 
            "counts":"int64", 
            "length":"int64", 
            "fraction":"float64"
        }
        )

    genomecov = genomecov[genomecov["chr"].isin(chromosomes)]
    bases_sequenced = sum(genomecov["counts"]*genomecov["depth"])
    length_chromosomes = sum(genomecov["length"].unique())
    DoC = bases_sequenced / length_chromosomes
    return(DoC)
  
def get_genomecov_file(file):
    filters_order = config["geno_calls"]["filters"]["order"]
    index_filter = filters_order.index("depth")

    if index_filter:
        for i in range(0, index_filter+1):
            file = file.replace("_"+sufixes[filters_order[i]], "")
    
    genomecov_name = "genomecov/" + file + ".genomecov"
    
    return genomecov_name

def get_depth_limit(boundary = "min", genomecov_path = ""):
    is_absolute = get_filter_value("depth", "absolute", False)
    if is_absolute:
        depth = get_filter_value("depth", boundary)
    else:
        DoC = get_DoC(genomecov_path)
        fraction = Fraction(get_filter_value("depth", boundary))
        depth =  fraction * DoC
    return depth
#-----------------------------------------------------------------------------#

sufix = build_sufix(filters, config)
print(sufix)

rule all:
    input:
        expand(
            "Filtered/chr{chr}/{bamlist}_" + sufix + "_chr{chr}.bcf", 
            chr = chromosomes, 
            bamlist = bamlists
            )


rule mask:
    input:
        bcf = lambda wildcards: get_input_file("mask", wildcards.file, wildcards.chr), 
        annotation = "masks/chr{chr}.bed.gz",
        header = "masks/header.txt"
    output:
        bcf = temp("Filtered/chr{chr}/{file}_mask"+mask_label+"_chr{chr}.bcf")
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mask_mem", attempt, 8),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mask_time", attempt, 1)
    log:
        "logs/filter_{file}_chr{chr}.log"
    shell:
        """        
        bcftools annotate -a {input.annotation} \
            -m MASK=strict \
            -h {input.header} \
            -c CHROM,POS,FROM,TO \
             {input.bcf} | \
            bcftools view -i 'MASK=1' \
                -Ob > {output.bcf} 
        """ 


rule cp_genomecov:
    input:
        # bcf = lambda wildcards: [get_input_file("depth", wildcards.file, chr)
        #             for chr in chromosomes],
        genomecov = lambda wildcards: get_genomecov_file(wildcards.file)
    output:
        genomecov = "DoC/{file}.genomecov"
    log:
        "logs/filter_{file}.log"
    shell:
        """
        mkdir -p DoC 
        cp {input.genomecov} {output.genomecov} 
        """

rule depth:
    input:
        bcf = lambda wildcards: get_input_file("depth", wildcards.file, wildcards.chr),
        genomecov = "DoC/{file}.genomecov"
    output:
        bcf = temp("Filtered/chr{chr}/{file}_dp_chr{chr}.bcf")
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("depth_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("depth_time", attempt, 2)
    params:
        min_depth = lambda wildcards: get_depth_limit(
            boundary="min", 
            genomecov_path=f"DoC/{wildcards.file}.genomecov"
            ),
        max_depth = lambda wildcards: get_depth_limit(
            boundary="max", 
            genomecov_path=f"DoC/{wildcards.file}.genomecov"
            )
    log:
        "logs/filter_{file}_chr{chr}.log"
    threads: 1
    shell:
        """
        bcftools filter \
            --threads {threads} \
            -Ob \
            --exclude "(SUM(DP4)< {params.min_depth} | SUM(DP4) > {params.max_depth} )" \
            {input.bcf} \
            > {output.bcf} \
        """

rule genotype_quality:
    input:
        bcf = lambda wildcards: get_input_file("gq", wildcards.file, wildcards.chr), 
    output:
        bcf = "Filtered/chr{chr}/{file}_"+f"GQ{gq}"+"_chr{chr}.bcf"
    log:
        "logs/filter_{file}_chr{chr}.log"
    params:
        gq = get_filter_value("gq", "min")
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("genoqual_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("genoqual_time", attempt, 2)
    shell:
        """
        bcftools filter \
            --threads {threads} \
                --exclude "QUAL <= {params.gq} "\
            -Ob {input.bcf} \
                > {output.bcf} 
        """

rule repeats:
    input:
        bcf = lambda wildcards: get_input_file("repeats", wildcards.file, wildcards.chr), 
        annotation = "repeats/chr{chr}.bed.gz",
        header = "repeats/header.txt"
    output:
        bcf = temp("Filtered/chr{chr}/{file}_"+f"repeats{repeats_label}"+"_chr{chr}.bcf")
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("repeats_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("repeats_time", attempt, 2)
    log:
        "logs/filter_{file}_chr{chr}.log"
    shell:
        """
        bcftools annotate -a {input.annotation} \
            -m REPEATS=repeats \
            -h {input.header} \
            -c CHROM,POS,FROM,TO \
             {input.bcf} | \
            bcftools filter --exclude REPEATS=1 \
                -Ob > {output.bcf} 
        """