# Snakemake for Error analysis in ANGSD
myGit = "/users/dcruzdav/data/Git/Botocudos-scripts/"
config_name = "multiple_purposes.yaml"
configfile: config_name

# Default files from angsd
outgroup = config["Outgroup_file"] if "Outgroup_file" in config.keys() else "hg19ancNoChr.fa"
perfect = config["Perfect_file"] if "Perfect_file" in config.keys() else "NA12778.fa"
threads = config["ancError_threads"] if "ancError_threads" in config.keys() else 1

bamlists = list(config["bamlists"].keys())
wildcard_constraints:
    bamlist = "(" + "|".join([b for b in bamlists]) + ")"

# This magic snakefile will make a list with all the bamfiles
# to run, plus sublists with bamfiles per population
# as specified in the config file
include: "make_bamlist.smk"


# rule all:
#     input:
#         bamlist = expand("{bamlist}.txt", bamlist = bamlists),
#         bamlist_group =  ["Bams_group/{bamlist}/{group}.txt".format(group = group, bamlist = bamlist) 
#                             for bamlist in list(config["bamlists"].keys()) 
#                             for group in config["bamlists"][bamlist]["paths"]],
#         error = ["{bamlist}/{group}/{group}_ancErr".format(bamlist = bamlist, group = group) 
#                 for bamlist in list(config["bamlists"].keys())
#                 for group in config["bamlists"][bamlist]["paths"]]


rule do_fasta:
    input:
        bam = "{file}.bam",
    output:
        fasta = "{file}.fa.gz",
    shell:
        """
        angsd -doFasta 2 -i {input.bam} -out {wildcards.file} -doCounts 1
        """
rule gunzip_fasta:
    input:
        "{file}.fa.gz"
    output:
        "{file}.fa"
    shell:
        """
        gunzip {input}
        """

rule index_fasta:
    input:
        fasta = "{file}.fa"
    output:
        index = "{file}.fa.fai"
    shell:
        """
        samtools faidx {input.fasta}
        """

rule do_AncError:
    input:
        outgroup = outgroup.split(".")[0] + ".fa",
        outgroup_idx = outgroup.split(".")[0] + ".fa.fai",
        perfect = perfect.split(".")[0] + ".fa",
        perfect_idx = perfect.split(".")[0] + ".fa.fai",
        bam_group = "Bams_group/{group}.txt"
    output:
        error = "{group}/{group}_ancErr.ancError"
    threads:
        threads
    shell:
        """
        angsd -doAncError 1 -anc {input.outgroup} \
        -ref {input.perfect} \
        -out {wildcards.group}/{wildcards.group}_ancErr \
        -bam {input.bam_group} \
        -nThreads {threads}
        """