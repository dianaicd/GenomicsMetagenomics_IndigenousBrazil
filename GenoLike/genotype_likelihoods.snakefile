# Snakefile to compute genotype likelihoods for
# a set of bam files,
# convert a panel (bed|ped|vcf) to genotype likelihoods
# and merge likelihoods for (panel, bam) files.
configfile: "genolike.yaml"
import os, glob

bamlists = list(config["bamlists"].keys())
panels = list(config["panels"].keys()) 
git_path = "/home/dcruzdva/data/Git/Botocudos-scripts/"

chromosomes = [str(x) for x in range(1,23)] 

wildcard_constraints:
    extension = "(bed|vcf|ped)",
    Chr =   "|".join(chromosomes)  ,
    panel =  "(" + "|".join([p for p in panels])+ ")(?!_" + "|_".join([b for b in bamlists]) + ")" ,
    bamlist = ")".join(["(?<!"+p for p in panels])  + "_)" + "|".join([b for b in bamlists])

rule all:
    input:
        merged = expand("{panel}/{bamlist}_{panel}.beagle", panel = panels, bamlist = bamlists)

def expand_path(wildcards):
    paths = [list(config["bamlists"][l]["paths"].values()) for l in bamlists][0]
    #full_paths = [os.path.expanduser(p) for p in paths]
    bams = []
    for p in paths:
        myCommand = "ls " + p 
        files = os.popen(myCommand).read().split('\n')
        files.remove('')
        [bams.append(individual_path) for individual_path in files]
    #bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

rule make_bamlist:
    input:
        expand_path 
    output:
        "{bamlist}.txt"
    log:
        "logs/{bamlist}_make_bamlist.log"
    run:
        with open(output[0], 'w') as file:
            for line in input:
                file.write(line+"\n")

rule panel_to_beagle:
    input:
        panel = lambda wildcards: config["panels"][wildcards.panel]["path"]
    output:
        panel = "{panel}/{panel}.beagle",
        ids = "{panel}/{panel}_ids.txt",
        names = "{panel}/{panel}_names.txt"
    params:
        rmdamage = lambda wildcards: config["panels"][wildcards.panel]["GenoLike"]["rmdamage"]
    log:
        "logs/{panel}.log"
    shell:
        """
            extension=$(echo {input.panel}|rev |cut -f1 -d. |rev) ;
            if [ $extension == "bed" ] 
            then plink --recode vcf --bfile {wildcards.panel} --out {wildcards.panel} 
            elif [ $extension == "ped" ]
            then plink --recode vcf-fid --file {wildcards.panel} --out {wildcards.panel}
            fi

            if [ ! -e vcftogenolike.pl ] 
            then ln -s ~/data/Git/Botocudos-scripts/GenoLike/vcftogenolike.pl ./
            fi
            perl vcftogenolike.pl -i {wildcards.panel}.vcf  -rmdamage {params.rmdamage} -o {output.panel} 2>{log}
            cut -f 1 {output.panel} > {output.ids}
            bcftools query -l {wildcards.panel}.vcf > {output.names}
        """

rule genos_chr:
    input:
        panel = "{panel}/{panel}.beagle",
        bamlist = "{bamlist}.txt"
    output:
        sites = "{panel}/{bamlist}_{chr}_sites.txt",
        idx = "{panel}/{bamlist}_{chr}_sites.txt.idx",
        binary = "{panel}/{bamlist}_{chr}_sites.txt.bin",
        rf = "{panel}/chr{chr}_{bamlist}.txt",
        gl = "{panel}/{bamlist}_{chr}.beagle.gz"
    params:
        trim = lambda wildcards: config["bamlists"][wildcards.bamlist]["trim"],
        minQ = lambda wildcards: config["bamlists"][wildcards.bamlist]["minQ"],
        minmapQ = lambda wildcards: config["bamlists"][wildcards.bamlist]["minmapQ"]
    log:
        "logs/{bamlist}_{panel}_{chr}_genolike.log"
    threads:
        config["nThreads"]
    shell:
        """ 
        cut -f 1-3 {input.panel}| sed 1d |\
            grep "^{wildcards.chr}_" | sed 's/_/\t/' |sort -V > {output.sites} 
        cut -f1 {output.sites} |sort -V |uniq > {output.rf}
        angsd sites index {output.sites}
        angsd -GL 1 -out {wildcards.panel}/{wildcards.bamlist}_{wildcards.chr} -doGlf 2 -doMajorMinor 3  \
            -bam {input.bamlist} -minQ {params.minQ} -minmapQ {params.minmapQ} {params.trim}\
            -sites {output.sites} -rf {output.rf}  -checkbamheaders 0 \
            -nThreads {threads} &>{log}
        """

rule merge_chr:
    input:
        beagle_gz = lambda wildcards: expand("{panel}/{bamlist}_{chr}.beagle.gz",
                                             bamlist = bamlists, chr = chromosomes, panel = wildcards.panel)
    output:
        beagle = "{panel}/{bamlist}.beagle",
        header = temp("{panel}/{bamlist}_header.txt"),
        beagle_temp = temp("{panel}/{bamlist}_tmp.beagle")
    log:
        "logs/{bamlist}_{panel}_genolike.log"
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
        merged = "{panel}/{bamlist}_{panel}.beagle"
    params:
        homozygous = lambda wildcards: config["panels"][wildcards.panel]["GenoLike"]["homozygous"]
    log:
        "logs/{bamlist}_{panel}_genolike.log"
    shell:
        """
        if [ ! -e merge_genos.pl ]
        then ln -s ~/data/Git/Botocudos-scripts/GenoLike/merge_genos.pl ./
        fi 
        perl merge_genos.pl -g2 {input.panel}\
            -g1 {input.bamlist} -id {input.ids} \
            -o {output.merged} -homozygous {params.homozygous}       2>>{log}
        """