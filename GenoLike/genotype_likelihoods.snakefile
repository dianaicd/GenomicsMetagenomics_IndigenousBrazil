# Snakefile to compute genotype likelihoods for
# a set of bam files,
# convert a panel (bed|ped|vcf) to genotype likelihoods
# and merge likelihoods for (panel, bam) files.
configfile: "multiple_purposes.yaml"
import os, glob
include: "parse_resources.smk"
localrules: make_bamlist, merge_chr, merge_genos, plink_to_vcf, vcf_to_beagle

bamlists = list(config["geno_like"]["bamlists"].keys())
panels = list(config["geno_like"]["panels"].keys()) 


chromosomes = [str(x) for x in range(1,23)] 

# wildcard_constraints:
#     extension = "(bed|vcf|ped)",
#     chr =   "|".join(chromosomes)  ,
#     panel =  "(" + "|".join([p for p in panels])+ ")(?!_" + "|_".join([b for b in bamlists]) + ")" ,
#     bamlist = ")".join(["(?<!"+p for p in panels])  + "_)" + "|".join([b for b in bamlists])

def expand_path(wildcards):
    myDict = config["geno_like"]["bamlists"][wildcards.bamlist]["paths"]
    paths = [myDict[group] for group in myDict.keys()]
    #full_paths = [os.path.expanduser(p) for p in paths]
    bams = []
    for p in paths:
        myCommand = "shopt -s extglob ; ls " + p 
        files = os.popen(myCommand).read().split('\n')
        files.remove('')
        [bams.append(individual_path) for individual_path in files]
    #bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

localrules: all,make_bamlist,plink_to_vcf,ped_to_vcf,vcf_to_beagle,angsd_sites,merge_chr,merge_genos
rule all:
    input:
        merged = expand("{panel}/{panel}_{bamlist}.beagle", panel = panels, bamlist = bamlists)

rule make_bamlist:
    input:
        expand_path 
    output:
        "{bamlist}.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("make_bamlist_mem", attempt, 1),
        runtime=lambda wildcards, attempt: get_runtime_alloc("make_bamlist_time", attempt, 1)

    log:
        "logs/{bamlist}_make_bamlist.log"
    run:
        with open(output[0], 'w') as file:
            for line in input:
                file.write(line+"\n")


rule plink_to_vcf:
    input:
        panel = lambda wildcards: config["geno_like"]["panels"][wildcards.panel]["path"]#"{panel}/{panel}.bed"
    output:
        "{panel}/{panel}.vcf"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("plink_to_vcf_mem", attempt, 1),
        runtime=lambda wildcards, attempt: get_runtime_alloc("plink_to_vcf_time", attempt, 1)
    params:
        out = "{panel}/{panel}",
        basename = lambda wildcards: config["geno_like"]["panels"][wildcards.panel]["path"].replace(".bed", "")
        # basename = "{panel}/{panel}"
    shell:
        """
        plink --recode vcf --bfile {params.basename} --out {params.out} 
        
        """

# rule ped_to_vcf:
#     input:
#         panel = "{panel}/{panel}.ped"
#     output:
#         "{panel}/{panel}.vcf"
#     resources:
#         memory=lambda wildcards, attempt: get_memory_alloc("plink_to_vcf_mem", attempt, 1),
#         runtime=lambda wildcards, attempt: get_runtime_alloc("plink_to_vcf_time", attempt, 1)
#     params:
#         basename = "{panel}/{panel}",
#         out = "{panel}/{panel}",
#     shell:
#         """
#         plink --recode vcf --file {params.basename} --out {params.out} 
        
#         """

rule vcf_to_beagle:
    input:
        panel = "{panel}/{panel}.vcf"
    output:
        panel = "{panel}/{panel}.beagle",
        ids = "{panel}/{panel}_ids.txt",
        names = "{panel}/{panel}_names.txt"
    params:
        rmdamage = lambda wildcards: config["geno_like"]["panels"][wildcards.panel]["rmdamage"]
    log:
        "logs/{panel}.log"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("panel_to_beagle_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("panel_to_beagle_time", attempt, 6)
    shell:
        """
        if [ ! -e vcftogenolike.pl ] 
            then ln -s ~/data/Git/Botocudos-scripts/GenoLike/vcftogenolike.pl ./
        fi

        perl vcftogenolike.pl -i {input.panel}  \
            -rmdamage {params.rmdamage} \
            -o {output.panel} 2>{log}

        cut -f 1 {output.panel} > {output.ids}
        bcftools query -l {input.panel} > {output.names}
        """

rule angsd_sites:
    input:
        panel = "{panel}/{panel}.beagle",
    output:
        sites = "{panel}/chr{chr}_sites.txt",
        rf = "{panel}/chr{chr}.txt",
        # idx = "{panel}/chr{chr}_sites.txt.idx",
        # binary = "{panel}/chr{chr}_sites.txt.bin"
    params:
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("angsd_sites_mem", attempt, 1),
        runtime = lambda wildcards, attempt: get_runtime_alloc("angsd_sites_time", attempt, 1)
    shell:
        """
        set +e
        cut -f 1-3 {input.panel}| sed 1d |\
            grep "^{wildcards.chr}_" | sed 's/_/\t/' |sort -V > {output.sites} 
        cut -f1 {output.sites} |sort -V |uniq > {output.rf}
        angsd sites index {output.sites}
        exitcode=$? 
        if [ $exitcode -eq 1 ] ; then exit 1; else exit 0 ; fi
        """

rule genotype_likelihoods:
    input:
        sites = "{panel}/chr{chr}_sites.txt",
        rf = "{panel}/chr{chr}.txt",
        bamlist = "{bamlist}.txt"
    output:
        # idx = "{panel}/chr{chr}_sites.txt.idx",
        # binary = "{panel}/chr{chr}_sites.txt.bin",
        gl = temp("{panel}/{bamlist}_chr{chr}.beagle.gz"),
        arg = temp("{panel}/{bamlist}_chr{chr}.arg")
    params:
        trim = lambda wildcards: config["geno_like"]["bamlists"][wildcards.bamlist]["trim"] if "trim" in config["geno_like"]["bamlists"][wildcards.bamlist].keys() else "",
        minQ = lambda wildcards: config["geno_like"]["bamlists"][wildcards.bamlist]["minQ"] if "minQ" in config["geno_like"]["bamlists"][wildcards.bamlist].keys() else 20,
        minmapQ = lambda wildcards: config["geno_like"]["bamlists"][wildcards.bamlist]["minmapQ"] if "minmapQ" in config["geno_like"]["bamlists"][wildcards.bamlist].keys() else 30
    log:
        "logs/{bamlist}_{panel}_chr{chr}_genolike.log"
    threads:
        config["geno_like"]["nThreads"]
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("genos_chr_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("genos_chr_time", attempt, 24)
    shell:
        """ 
        sleep 65
        touch {wildcards.panel}/chr{wildcards.chr}_sites.txt.{{idx,bin}}

        angsd -GL 1 \
            -out {wildcards.panel}/{wildcards.bamlist}_chr{wildcards.chr} \
            -doGlf 2 -doMajorMinor 3  \
            -bam {input.bamlist} \
            -minQ {params.minQ} \
            -minmapQ {params.minmapQ} \
            {params.trim} \
            -sites {input.sites} \
            -rf {input.rf}  \
            -checkbamheaders 0 \
            -nThreads {threads} &>{log}
        """

rule merge_chr:
    input:
        beagle_gz = lambda wildcards: ["{panel}/{bamlist}_chr{chr}.beagle.gz".format(
                                             bamlist = wildcards.bamlist, chr = chr, panel = wildcards.panel)
                                             for chr in chromosomes]
    output:
        beagle = "{panel}/{bamlist}.beagle",
        header = temp("{panel}/{bamlist}_header.txt"),
        beagle_temp = temp("{panel}/{bamlist}_tmp.beagle")
    log:
        "logs/{bamlist}_{panel}_genolike.log"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merge_chr_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merge_chr_time", attempt, 1)
    shell:
        """
        cat {input.beagle_gz} > {output.beagle_temp}.gz        2>{log}
        gunzip {output.beagle_temp}.gz            2>>{log}
        head -n1 {output.beagle_temp} > {output.header}     2>>{log}
        sed -i '/marker/d' {output.beagle_temp}             2>>{log}
        cat {output.header} {output.beagle_temp} > {output.beagle}  2>>{log}
        """

rule merge_genos:
    input:
        panel = "{panel}/{panel}.beagle",
        bamlist = "{panel}/{bamlist}.beagle",
        ids = "{panel}/{panel}_ids.txt"
    output:
        merged = "{panel}/{panel}_{bamlist}.beagle"
    params:
        homozygous = lambda wildcards: config["geno_like"]["panels"][wildcards.panel]["homozygous"]
    log:
        "logs/{panel}_{bamlist}_genolike.log"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merge_genos_mem", attempt, 32),
        runtime = lambda wildcards, attempt: get_runtime_alloc("merge_genos_time", attempt, 1)
    shell:
        """
        if [ ! -e merge_genos.pl ]
        then ln -s ~/data/Git/Botocudos-scripts/GenoLike/merge_genos.pl ./
        fi 
        perl merge_genos.pl -g1 {input.panel}\
            -g2 {input.bamlist} -id {input.ids} \
            -o {output.merged} -homozygous {params.homozygous}       2>>{log}
        """
