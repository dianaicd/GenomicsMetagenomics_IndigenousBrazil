# Snakemake for Error analysis in ANGSD
myGit = "/users/dcruzdav/data/Git/Botocudos-scripts/"
config_name = "multiple_purposes.yaml"
configfile: config_name
import glob

def param_default(myPar, myDict, myDef):
    if myPar in myDict.keys():
        param = myDict[myPar]
    else:
        param = myDef
    return param

def fasta_name(outgroup_name, outgroup_file, index = False):
    format_file = outgroup_file.split(".")[-1]
    if format_file == "fa":
        name = outgroup_file
    elif format_file == "bam":
        name = f"{outgroup_name}_baseQ{minQ}_mapQ{mapQ}.fa"
    if index:
        name = f"{name}.fai"

    return name

outgroup_name = [key for key,value in config["error_angsd"]["Outgroup_file"].items()][0]
outgroup_file = config["error_angsd"]["Outgroup_file"][outgroup_name]


perfect_name = [key for key,value in config["error_angsd"]["Perfect_file"].items()][0]
perfect_file = config["error_angsd"]["Perfect_file"][perfect_name]

rscript_plot_path =  param_default("Rscript_err_angsd",config["error_angsd"], 
                                    "/software/UHTS/Analysis/ANGSD/0.931/R/estError.R")

threads = config["error_angsd"]["ancError_threads"] if "ancError_threads" in config["error_angsd"].keys() else 1

bamlists = list(config["error_angsd"]["bamlists"].keys())
wildcard_constraints:
    bamlist = "(" + "|".join([b for b in bamlists]) + ")"

minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30)
mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30)

my_files = {
    perfect_name:config["error_angsd"]["Perfect_file"][perfect_name],
    outgroup_name:config["error_angsd"]["Outgroup_file"][outgroup_name]
}

#
# Get chromosome size and break it in blocks
chromosomes = [i for i in range(1,23)]
blockSize = config["error_angsd"]["blockSize"] if "blockSize" in config["error_angsd"] else "1e7"

with open(config["ref_genome"] + ".fai", 'r') as index:
    chr_size = {}
    for line in index.readlines():
        chr,size = line.split()[0:2]
        chr_size[chr] = size
    lower_chr = {}
    upper_chr = {}
    for chr in chromosomes:
        chr = str(chr)

        lower_chr[chr] = [i for i in range(1, int(chr_size[chr]), int(float(blockSize)))]
        upper_chr[chr] = [i - 1 for i in lower_chr[chr][1:]]
        upper_chr[chr].append(chr_size[chr])

# This magic snakefile will make a list with all the bamfiles
# to run, plus sublists with bamfiles per population
# as specified in the config file
myDict = config["error_angsd"]["bamlists"]
include: "make_bamlist.smk"


rule all:
    input:
        bamlist = expand("{bamlist}/{bamlist}.txt", bamlist = bamlists),
        error = expand("{bamlist}/{bamlist}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancError",
                        bamlist = bamlists, 
                        perfect = perfect_name,
                        outgroup = outgroup_name
                        ),
        error_txt = expand("{bamlist}/{bamlist}_perfect.{perfect}_outgroup.{outgroup}_error.txt", 
                        bamlist = bamlists,
                        perfect = perfect_name,
                        outgroup = outgroup_name
                        )


rule do_fasta:
    input:
        bam = lambda wildcards: my_files[wildcards.file]
    output:
        fasta = "{file}_baseQ{minQ}_mapQ{mapQ}.fa.gz".format(minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30),
                                                            mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
                                                            file = "{file}")
    log:
        "logs/do_fasta_{file}_baseQ{minQ}_mapQ{mapQ}.log".format(minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30),
                                                            mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
                                                            file = "{file}")
    params:
        minQ =  param_default("minQ_perfect", config["error_angsd"], myDef = 30),
        minMapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
        prefix_out = "{file}_baseQ{minQ}_mapQ{mapQ}".format(minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30),
                                                            mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
                                                            file = "{file}")
    resources:
        runtime=8*60
    shell:
        """
        angsd -doFasta 2 -i {input.bam} -out {params.prefix_out} \
        -minQ {params.minQ} -minMapQ {params.minMapQ} -doCounts 1 >{log}
        """

rule gunzip_fasta:
    input:
        "{file}.fa.gz"
    output:
        "{file}.fa"
    log:
        "logs/gunzip_fasta_{file}.log"
    shell:
        """
        gunzip {input} &>{log}
        """

rule index_fasta:
    input:
        fasta = "{file}.fa"
    output:
        index = "{file}.fa.fai"
    log:
        "logs/index_fasta_{file}.log"
    shell:
        """
        samtools faidx {input.fasta} &>{log}
        """


rule do_AncError:
    input:
        outgroup = lambda wildcards: fasta_name(outgroup_name = wildcards.outgroup, outgroup_file = outgroup_file),
        outgroup_idx = lambda wildcards: fasta_name(outgroup_name = wildcards.outgroup, outgroup_file = outgroup_file, index = True),
        perfect = lambda wildcards: f"{wildcards.perfect}_baseQ{minQ}_mapQ{mapQ}.fa",
        perfect_idx = lambda wildcards: f"{wildcards.perfect}_baseQ{minQ}_mapQ{mapQ}.fa.fai",
        bam_group = "{group}/{group}.txt"
    output:
        error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr_chr{chr}_{start}_{end}.ancError"
    threads:
        threads
    params:
        minQ = config["error_angsd"]["BaseQuality"] if "BaseQuality" in config["error_angsd"].keys() else 30,
        minMapQ = config["error_angsd"]["MapQuality"] if "MapQuality" in config["error_angsd"].keys() else 30,
        basename = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr_chr{chr}_{start}_{end}".format(
            group = "{group}", perfect = "{perfect}", outgroup = "{outgroup}",
            chr = "{chr}", start = "{start}", end = "{end}"
            )
    resources:
        mem = 1*512,
        runtime = 1*10
    log:
        "logs/do_AncError_{group}_{perfect}_{outgroup}_chr{chr}_{start}_{end}.log"
    shell:
        """
        angsd -doAncError 1 -anc {input.outgroup} \
        -ref {input.perfect} \
        -out {params.basename} \
        -bam {input.bam_group} \
        -nThreads {threads} \
        -minQ {params.minQ} -minMapQ {params.minMapQ} \
        -r    {wildcards.chr}:{wildcards.start}-{wildcards.end} &>{log}
        """

def expand_errors(group, perfect, outgroup, chr):
    outputs = [
        "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr_chr{chr}_{start}_{end}.ancError".format(
        chr = chr,
        start = lower_chr[str(chr)][i], 
        end = upper_chr[str(chr)][i],
        group = group,
        outgroup = outgroup,
        perfect = perfect
        ) 
        for i in range(0, len(lower_chr[str(chr)]))
        ]

    return(outputs)


rule merge_ancErr_by_chr:
    input:
        bam_group = "{group}/{group}.txt",
        anc_error = lambda wildcards: expand_errors(
                wildcards.group,
                wildcards.perfect,
                wildcards.outgroup,
                wildcards.chr
                )
    output:
        anc_error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr_chr{chr}.ancError"
    wildcard_constraints:
        chr = "(" + "|".join([str(i) for i in chromosomes]) + ")"
    run:
        import numpy as np
        
        with open(input.bam_group, "r") as file:
            n_ind = len(file.readlines())

        final_matrix = np.zeros((n_ind, 125))
        for err in input.anc_error:
            current_matrix = np.loadtxt(err)
            final_matrix = final_matrix + current_matrix

        np.savetxt(output.anc_error, final_matrix,
        header = f"Chr:\t{wildcards.chr}", delimiter = "\t")


rule final_error_chr:
    input:
        anc_error = expand(
                    "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr_chr{chr}.ancError",
                    group = "{group}",
                    outgroup = "{outgroup}",
                    perfect = "{perfect}",
                    chr = chromosomes
                )
    output:
        anc_error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancErrorChr"
    shell:
        """
        for file in {input}
        do
            cat $file | sed 's/#Chr/Chr/'
        done > {output} 
        """

rule merge_ancErr_by_Ind:
    input:
        bam_group = "{group}/{group}.txt",
        anc_error = expand(
            "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr_chr{chr}.ancError",
            group = "{group}",
            perfect = "{perfect}",
            outgroup = "{outgroup}",
            chr = chromosomes
            )
    output:
        anc_error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancError"
    run:
        import numpy as np
        
        with open(input.bam_group, "r") as file:
            n_ind = len(file.readlines())

        final_matrix = np.zeros((n_ind, 125))
        for err in input.anc_error:
            current_matrix = np.loadtxt(err)
            final_matrix = final_matrix + current_matrix

        np.savetxt(output.anc_error, final_matrix, delimiter = "\t")


rule plot_error:
    input:
        bam_group = "{group}/{group}.txt",
        error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancError",
        error_chr = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancErrorChr"
    output:
        names = "{group}/{group}_{perfect}_{outgroup}_names.txt",
        err_txt = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_error.txt"
    params:
        name = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_error".format(group = "{group}", perfect = "{perfect}", outgroup = "{outgroup}")
    shell:
        """
        cat {input.bam_group} | awk 'BEGIN {{FS="/"}} {{print $NF}}' | sed 's/.hg19.bam//' > {output.names}
        Rscript {rscript_plot_path} file={input.error} doPng=T out="{params.name}" \
        indNames={output.names}
        """