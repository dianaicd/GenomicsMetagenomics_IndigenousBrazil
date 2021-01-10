# Snakefile version to merge datasets
configfile: "multiple_purposes.yaml"
import os, glob, math, subprocess
import numpy as np

include: "parse_resources.smk"

if "panels" not in globals():
    panels = list(config["ngs_admix"]["panels"].keys())
    panel_dict = config["ngs_admix"]["panels"]

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

rule all_ngsadmix:
    input:
        qopt = ["{k}/{panel}_k{k}.qopt".format( panel = panel, k = str(k) ) for panel in panels for k in num_clusters( panel )],
        fopt = ["{k}/{panel}_k{k}.fopt.gz".format( panel = panel, k = str(k) ) for panel in panels for k in num_clusters( panel ) ],
        commands =  expand("{panel}/best_runs_scp_{panel}.txt", panel = panels)

rule NGSadmix:
    input:
        beagle = lambda wildcards: panel_dict[wildcards.panel]["path"]
    output:
        qopt = "{k}/{panel}_k{k}_{rep}.qopt",
        log = "{k}/{panel}_k{k}_{rep}.log"
    log:
        "logs/{k}_{panel}_{rep}.log"
    threads:
        8
    resources:
        runtime = lambda wildcards, attempt: get_runtime_alloc( "NGSAdmix_time", 
                                                                attempt, 
                                                                12
                                                                ),
        memory = lambda wildcards, attempt: get_memory_alloc("NGSAdmix_mem", attempt, 8)
    shell:
        """
        NGSadmix -likes {input.beagle} -K {wildcards.k} -P {threads} \
            -o {wildcards.k}/{wildcards.panel}_k{wildcards.k}_{wildcards.rep} \
            -printInfo 1 \
             &> {log}
        """

rule parse_likelihood:
    input:
        log = lambda wildcards: expand_log(wildcards.panel, wildcards.k)
    output:
        likelihood = "{panel}/likelihoods_k{k}_{panel}.txt"
    wildcard_constraints:
        k = "[\d]+"
    log:
        "logs/parse_likelihood_{panel}_{k}.txt"
    shell:
        """
        myInput=({input.log})
        for file in ${{myInput[@]}} ; do
            likelihood=$(grep like $file | sed 's/.*like=/\\t/; s/ after.*//')
            rep=$(echo $file |sed 's/.*_//;s/.log//')
            echo "$rep$likelihood" 
        done > {output.likelihood} 2>{log}

        """

rule find_best_run:
    input:
        likelihoods_path = lambda wildcards: ["{panel}/likelihoods_k{k}_{panel}.txt".format(
                                                panel = wildcards.panel, 
                                                k = k)
                                                for k in range( min_K[ wildcards.panel ], max_K[ wildcards.panel ] ) ]
                                
    output:
        commands_scp = "{panel}/best_runs_scp_{panel}_k{k}.txt",
        best_run_qopt = "{k}/{panel}_k{k}.qopt",
        best_run_fopt = "{k}/{panel}_k{k}.fopt.gz",
        best_filter = "{k}/{panel}_k{k}.filter"
    log:
        "logs/find_likelihood_{panel}_{k}.txt"
    run:
        with open( output.commands_scp, "w") as scp_commands:
            for likelihoods_path in input.likelihoods_path:
                likelihoods = np.loadtxt( likelihoods_path )
                best_replicate = np.argmax( likelihoods[:, 1]) + 1
                k = likelihoods_path.split("/")[1].split("_")[1].replace("k", "")
                best_run = "{k}/{panel}_k{k}_{rep}.qopt".format(k = k, 
                                                                panel = wildcards.panel,
                                                                rep = str( best_replicate )
                                                                )

                command = "scp Axiom:" + os.getcwd() + "/" + best_run.replace( "_" + str(best_replicate ), "" ) + " ./\n"
                scp_commands.write(command)

                new_prefix = best_run.replace( f"_{best_replicate}.qopt", ".qopt" )

                for sufix in [".qopt", ".fopt.gz", ".filter"]:
                    command = "cp {old_file} {new_file} ".format( old_file = best_run.replace(".qopt", sufix),
                                                                new_file = new_prefix.replace(".qopt", sufix))
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
    shell:
        """
        cat {input.commands} > {output.commands}
        """

