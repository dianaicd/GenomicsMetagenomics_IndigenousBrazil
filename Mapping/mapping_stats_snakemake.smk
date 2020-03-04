configname="samples.yaml"
configfile: configname
mito=config["Mitochondrial"]
import re
import yaml

# Verify that all necessary files are present:
# BAM:
# 1 final per sample
# 1 rmdup per library
# 1 including dup per library
# Stats from AdapterRemoval:
# 1 per ID
def get_index_column(input, columnName, sep = " "):
    with open(input, 'r') as file:
        header = file.readline().replace("\n", "")
        columnIndex = [index for index,element in enumerate(header.split(sep)) if(element == columnName)][0]
    return(columnIndex)

def expand_input(depth = "ID", SM = False, LB = False, ID = False):
    # returns prefix for an input
    myPrefix = {}
    if depth == "SM":
        [myPrefix.update({sm+"/"+sm : 1}) for sm in config["Samples"].keys()]
    elif depth == "LB":
        if SM:
            [myPrefix.update({sm+"/"+lb+"/library_rmdup/"+lb : 1}) 
                for sm in config["Samples"].keys()
                for lb in list(config["Samples"][sm]["Libraries"].keys()) if sm == SM]
        else:
            [myPrefix.update({sm+"/"+lb+"/library_rmdup/"+lb : 1}) 
                for sm in config["Samples"].keys()
                for lb in list(config["Samples"][sm]["Libraries"].keys())]
    elif depth == "ID":
        if SM:
            if LB:
                [myPrefix.update({sm+"/"+lb+"/fastq_bam_filter/"+id : 1})
                for sm in config["Samples"].keys()
                    if sm == SM 
                for lb in list(config["Samples"][sm]["Libraries"].keys())
                    if lb == LB
                for id in config["Samples"][sm]["Libraries"][lb]["IDs"].split(" ")]
            else:
                [myPrefix.update({sm+"/"+lb+"/fastq_bam_filter/"+id : 1})
                for sm in config["Samples"].keys()
                    if sm == SM 
                for lb in list(config["Samples"][sm]["Libraries"].keys())
                for id in config["Samples"][sm]["Libraries"][lb]["IDs"].split(" ")]
        else:
            [myPrefix.update({sm+"/"+lb+"/fastq_bam_filter/"+id : 1})
                for sm in config["Samples"].keys()
                for lb in list(config["Samples"][sm]["Libraries"].keys())
                for id in config["Samples"][sm]["Libraries"][lb]["IDs"].split(" ")]

    return(list(myPrefix.keys()))

IDs = expand_input(depth = "ID")
LBs = expand_input(depth = "LB")
SMs = expand_input(depth = "SM")

wildcard_constraints:
    file =   "(" + "|".join([lb for lb in LBs]) + ")",
    sample =  "(" + "|".join([sm.split("/")[-1] for sm in SMs]) + ")"
# print(SMs)
# print(len(SMs))
# print(LBs)
# print(len(LBs))
# print(IDs)
# print(len(IDs))

rule all:
    input:
        summary = expand("{prefix}.summary", prefix = SMs)

rule index_bam:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    shell:
        "samtools index {input}"
rule flagstat:
    input:
        "{file}.bam"
    output:
        "{file}_flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"

rule idxstats:
    input:
        bam = "{sample}/{library}/library_rmdup/{library}.bam",
        bai = "{sample}/{library}/library_rmdup/{library}.bam.bai"
    output:
        "{sample}/{library}/library_rmdup/{library}_idxstats.txt"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule get_chrs:
    input:
        "{file}.bam"
    run:
        chromosomes = shell('samtools view -H {input} |grep "^@SQ" |sed "s/SN://" |cut -f 2')
        print(chromosomes)
        print(type(chromosomes))

rule bedtools:
    input:
        "{sample}/{library}/library_rmdup/{library}.bam"
    output:
        "{sample}/{library}/library_rmdup/{library}.genomecov"
    shell:
        "bedtools genomecov -ibam {input} > {output}"

rule length:
    input:
        nuc="{sample}/{library}/library_rmdup/{library}.bam",
        bai="{sample}/{library}/library_rmdup/{library}.bam.bai"
    output:
        nuc="{sample}/{library}/library_rmdup/{library}.length",
        mito="{sample}/{library}/library_rmdup/{library}.mito.length"
    shell:
        "cat {input.nuc} |"
        "python ~/data/Git/Botocudos-scripts/DataQuality/read_length.py -o {output.nuc} ;"
        "samtools view -b {input.nuc} {mito} |"
        "python ~/data/Git/Botocudos-scripts/DataQuality/read_length.py -o {output.mito} "

def genomecov_input(wildcards):
    return([output + ".genomecov" for output in expand_input(depth = "LB", SM = wildcards.sample)])
def rmdup_input(wildcards):
    return([output + ".stats" for output in expand_input(depth = "LB", SM = wildcards.sample)])
def idxstats_input(wildcards):
    return([output + "_idxstats.txt" for output in expand_input(depth = "LB", SM = wildcards.sample)])
def idxstats_input(wildcards):
    return([output + "_idxstats.txt" for output in expand_input(depth = "LB", SM = wildcards.sample)])
def nuclen_input(wildcards):
    return([output + ".length" for output in expand_input(depth = "LB", SM = wildcards.sample)])
def mitolen_input(wildcards):
    return([output + ".mito.length" for output in expand_input(depth = "LB", SM = wildcards.sample)])

rule summary:
    input:
        bam="{sample}/{sample}.bam",
        #adapter = "{sample}/{library}/fastq_files_trim/{library}.settings", #theirChildrenSettings,
        rmdup = rmdup_input,
        genomecov= genomecov_input,
        idxstats = idxstats_input,
        nucLength= nuclen_input,
        mitoLength= mitolen_input
    output:
        "{sample}/{sample}.summary"
    shell:
        """
        python ~/data/Git/Botocudos-scripts/DataQuality/summary_from_yaml.py \
        --sample {wildcards.sample} --output {output} --mitochondrial {mito} --config_yaml {configname}
        """