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
    return(param)


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

my_files = {
    perfect_name:config["error_angsd"]["Perfect_file"][perfect_name],
    outgroup_name:config["error_angsd"]["Outgroup_file"][outgroup_name]
}
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
        runtime=4*60
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
        outgroup = lambda wildcards: "{outgroup}_baseQ{minQ}_mapQ{mapQ}.fa".format(
            minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30),
            mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
            outgroup = wildcards.outgroup
            ),
        outgroup_idx = lambda wildcards: "{outgroup}_baseQ{minQ}_mapQ{mapQ}.fa.fai".format(
            minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30),
            mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
            outgroup = wildcards.outgroup
            ),
        perfect = lambda wildcards: "{perfect}_baseQ{minQ}_mapQ{mapQ}.fa".format(
            minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30),
            mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
            perfect = wildcards.perfect
            ),
        perfect_idx = lambda wildcards: "{perfect}_baseQ{minQ}_mapQ{mapQ}.fa.fai".format(
            minQ = param_default("minQ_perfect", config["error_angsd"], myDef = 30),
            mapQ = param_default("minMapQ_perfect", config["error_angsd"], myDef = 30),
            perfect = wildcards.perfect
            ),
        bam_group = "{group}/{group}.txt"
    output:
        error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancError"
    threads:
        threads
    params:
        minQ = config["error_angsd"]["BaseQuality"] if "BaseQuality" in config["error_angsd"].keys() else 30,
        minMapQ = config["error_angsd"]["MapQuality"] if "MapQuality" in config["error_angsd"].keys() else 30,
        basename = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr".format(group = "{group}", perfect = "{perfect}", outgroup = "{outgroup}")
    resources:
        mem = 1*1024,
        runtime = 4*60
    log:
        "logs/do_AncError_{group}_{perfect}_{outgroup}.log"
    shell:
        """
        angsd -doAncError 1 -anc {input.outgroup} \
        -ref {input.perfect} \
        -out {params.basename} \
        -bam {input.bam_group} \
        -nThreads {threads} \
        -minQ {params.minQ} -minMapQ {params.minMapQ} &>{log}
        """

rule plot_error:
    input:
        bam_group = "{group}/{group}.txt",
        error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancError"
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