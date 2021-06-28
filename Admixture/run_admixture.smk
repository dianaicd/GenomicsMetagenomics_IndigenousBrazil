# Snakefile version to merge datasets
configfile: "multiple_purposes.yaml"
import os, glob, math, subprocess
import numpy as np

include: "parse_resources.smk"

if "panels" not in globals():
    panels = list(config["admixture"]["panels"].keys())
    panel_dict = config["admixture"]["panels"]

min_K = {panel:panel_dict[panel]["min_K"] for panel in panels}
max_K = {panel:panel_dict[panel]["max_K"]+1 for panel in panels}
replicates = {panel:panel_dict[panel]["replicates"] for panel in panels}

def num_clusters(panel):
    K = panel_dict[panel]["max_K"]
    min_K = panel_dict[panel]["min_K"] if "min_K" in panel_dict[panel].keys() else 2
    nClusters = [x for x in range(min_K, K+1)]
    return(nClusters)

def num_replicates(panel):
        rep = panel_dict[panel]["replicates"]
        nRep = [x for x in range(1, rep+1)]
        return(nRep)

def expand_log(panel, k):
    myLogs = [log for log in expand("{k}/{panel}_k{k}_{rep}.log", panel = panel, k = k, rep = num_replicates(panel) )]
    return(myLogs)


#-------------------------------------------------------------------------------------------------#
# Rules start here
localrules: link_file,all_ngsadmix,parse_likelihood,find_best_run,merge_commands

rule all_ngsadmix:
    input:
        qopt = ["{k}/{panel}_k{k}.Q".format( panel = panel, k = str(k) ) for panel in panels for k in num_clusters( panel )],
        fopt = ["{k}/{panel}_k{k}.P".format( panel = panel, k = str(k) ) for panel in panels for k in num_clusters( panel ) ],
        commands =  expand("{panel}/best_runs_scp_{panel}.txt", panel = panels)

rule link_file:
    input:
        bed = lambda wildcards: panel_dict[wildcards.panel]["path"]
    output:
        bed = temp("{k}/{panel}_k{k}_{rep}.bed"),
        bim = temp("{k}/{panel}_k{k}_{rep}.bim"),
        fam = temp("{k}/{panel}_k{k}_{rep}.fam")
    shell:
        """
        bim=$(echo {input.bed} | sed 's/.bed/.bim/')
        fam=$(echo {input.bed} | sed 's/.bed/.fam/')
        ln -s $(realpath {input.bed}) {output.bed}
        ln -s $(realpath $bim) {output.bim}
        ln -s $(realpath $fam) {output.fam}
        """

rule admixture:
    input:
        bed = "{k}/{panel}_k{k}_{rep}.bed",
        bim = "{k}/{panel}_k{k}_{rep}.bim",
        fam = "{k}/{panel}_k{k}_{rep}.fam"
    output:
        Q = temp("{k}/{panel}_k{k}_{rep}.Q"),
        P = temp("{k}/{panel}_k{k}_{rep}.P"),
        log = temp("{k}/{panel}_k{k}_{rep}.log")
    threads:
        4
    resources:
        runtime = lambda wildcards, attempt: get_runtime_alloc( "admixture_time", 
                                                                attempt, 
                                                                12
                                                                ),
        memory = lambda wildcards, attempt: get_memory_alloc("admixture_mem", attempt, 8)
    shell:
        """
        admixture {input.bed} {wildcards.k} -j{threads} &> {output.log}
        mv {wildcards.panel}_k{wildcards.k}_{wildcards.rep}.{wildcards.k}.Q {output.Q}
        mv {wildcards.panel}_k{wildcards.k}_{wildcards.rep}.{wildcards.k}.P {output.P}
        """

rule parse_likelihood:
    input:
        log = lambda wildcards: expand("{k}/{panel}_k{k}_{rep}.log",
                                        k = "{k}",
                                        panel = "{panel}",
                                        rep = range(replicates[wildcards.panel])
                                        )
    output:
        likelihood = "{panel}/likelihoods_k{k}_{panel}.txt"
    wildcard_constraints:
        k = "[\d]+"
    log:
        "logs/{k}_{panel}.txt"
    resources:
        runtime = 15,
        memory = 512
    shell:
        """
        myInput=({input.log})
        for file in ${{myInput[@]}} ; do
            likelihood=$(grep Loglikelihood $file |tail -n1| cut -f2 -d " ")
            rep=$(echo $file |sed 's/.*_//;s/.log//')
            echo "$rep $likelihood" 
        done > {output.likelihood} 2>{log}

        """

rule find_best_run:
    input:
        likelihoods_path = "{panel}/likelihoods_k{k}_{panel}.txt",
        runs = lambda wildcards: expand("{k}/{panel}_k{k}_{rep}.{ext}", 
                                        rep = range(replicates[wildcards.panel]), 
                                        k = "{k}", 
                                        panel = "{panel}",
                                        ext = ["Q", "P"])
                                
    output:
        commands_scp = "{panel}/best_runs_scp_{panel}_k{k}.txt",
        best_run_qopt = "{k}/{panel}_k{k}.Q",
        best_run_fopt = "{k}/{panel}_k{k}.P"
    log:
        "logs/{k}_{panel}.txt"
    resources:
        runtime = 15,
        memory = 512
    run:
        with open( output.commands_scp, "w") as scp_commands:
            likelihoods = np.loadtxt( input.likelihoods_path )
            best_replicate = np.argmax( likelihoods[:, 1]) #+ 1
            k = input.likelihoods_path.split("/")[1].split("_")[1].replace("k", "")
            best_run = "{k}/{panel}_k{k}_{rep}.Q".format(k = k, 
                                                            panel = wildcards.panel,
                                                            rep = str( best_replicate )
                                                            )

            command = "scp Wally:" + os.getcwd() + "/" + best_run.replace(f"_{best_replicate}.Q", ".Q" ) + " ./\n"
            scp_commands.write(command)

            new_prefix = best_run.replace( f"_{best_replicate}.Q", ".Q" )

            for sufix in [".Q", ".P"]:
                command = "cp {old_file} {new_file} ".format( old_file = best_run.replace(".Q", sufix),
                                                            new_file = new_prefix.replace(".Q", sufix))
                print( command )
                os.system( command )


rule merge_commands:
    input:
        commands =  lambda wildcards: expand("{panel}/best_runs_scp_{panel}_k{k}.txt", 
                                                panel = wildcards.panel,
                                                k = [ k for k in range( min_K[ wildcards.panel ], 
                                                                        max_K[ wildcards.panel ] 
                                                                        )
                                                    ]
                                                )
    output:
        commands =  "{panel}/best_runs_scp_{panel}.txt"
    log:
        "logs/{panel}.txt"
    resources:
        runtime = 60,
        memory = 256
    shell:
        """
        cat {input.commands} | sort |uniq > {output.commands}
        """

