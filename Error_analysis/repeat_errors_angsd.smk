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

# Default files from angsd
outgroup_file = param_default("Outgroup_file", config["error_angsd"], "hg19ancNoChr.fa")

perfect_file = param_default("Perfect_file", config["error_angsd"], "NA12778.fa")
rscript_plot_path =  param_default("Rscript_err_angsd",config["error_angsd"], 
                                    "/software/UHTS/Analysis/ANGSD/0.931/R/estError.R")

threads = config["ancError_threads"] if "ancError_threads" in config.keys() else 1

bamlists = list(config["error_angsd"]["bamlists"].keys())
wildcard_constraints:
    bamlist = "(" + "|".join([b for b in bamlists]) + ")"

my_names = {perfect_file: perfect_file.split("/")[-1].replace(".bam", "").replace(".fa", ""),
            outgroup_file: outgroup_file.split("/")[-1].replace(".fa", "")}
my_files = { perfect_file.split("/")[-1].replace(".bam", "").replace(".fa", ""): perfect_file.replace(".bam", "").replace(".fa", ""),
             outgroup_file.split("/")[-1].replace(".fa", ""): outgroup_file.replace(".fa", "")
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
                        perfect = my_names[perfect_file], 
                        outgroup = my_names[outgroup_file]),
        error_txt = expand("{bamlist}/{bamlist}_perfect.{perfect}_outgroup.{outgroup}_error_rep{rep}.txt", 
                        bamlist = bamlists,
                        perfect = my_names[perfect_file], 
                        outgroup = my_names[outgroup_file],
                        rep = [i for i in range(1, 101)])

rule make_namelist:
    input:
        bam_group = "{group}/{group}.txt"
    output:
        names = "{group}/{group}_{perfect}_{outgroup}_names.txt"
    shell:
        """
        cat {input.bam_group} | awk 'BEGIN {{FS="/"}} {{print $NF}}' | sed 's/.hg19.bam//' > {output.names}
        """

rule estimate_error:
    input:
        names = "{group}/{group}_{perfect}_{outgroup}_names.txt",
        error = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_ancErr.ancError"
    output:
        err_txt = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_error_rep{rep}.txt"
    params:
        name = "{group}/{group}_perfect.{perfect}_outgroup.{outgroup}_error_rep{rep}".format(group = "{group}", 
                                                                                            perfect = "{perfect}",
                                                                                            outgroup = "{outgroup}", 
                                                                                            rep = "{rep}")
    shell:
        """
        Rscript {rscript_plot_path} file={input.error} doPng=T out="{params.name}" \
        indNames={input.names}
        """