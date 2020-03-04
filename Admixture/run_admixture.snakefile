# Snakefile version to merge datasets
configfile: "multiple_purposes.yaml"
import os,glob
# bamlists = list(config[""]["bamlists"].keys())
# print(bamlists)  #config["samples"].values()
panels = list(config["ngs_admix"]["panels"].keys())
# print(panels)

def num_clusters(panel):
    K = config["ngs_admix"]["panels"][panel]["max_K"]
    min_K = config["ngs_admix"]["panels"][panel]["min_K"] if "min_K" in config["ngs_admix"]["panels"][panel].keys() else 2
    nClusters = [x for x in range(min_K, K+1)]
    return(nClusters)

def num_replicates(panel):
        rep = config["ngs_admix"]["panels"][panel]["replicates"]
        nRep = [x for x in range(1, rep+1)]
        return(nRep)

def expand_results(panels):
    myResults = [group for p in panels for group in expand("{k}/{panel}_k{k}_{rep}.qopt", panel = p, k = num_clusters(p), rep = num_replicates(p) )  for item in group]
    return(myResults)

def expand_log(panel, k):
    myLogs = [log for log in expand("{k}/{panel}_k{k}_{rep}.log", panel = panel, k = k, rep = num_replicates(panel) )]
    return(myLogs)

def give_likelihoods(panel):
    K = num_clusters(panel)
    myFiles = expand("likelihoods_k{k}_{panel}.txt",
    k = K, panel = panel)
    return(myFiles)

rule all:
    input:
        # qopt = expand("{k}/{panel}_k{k}_{rep}.qopt", k = num_clusters(), 
                        # panel = panels, rep = nRep)
        qopt = expand_results(panels),
        likelihoods = [file for panel in panels for file in give_likelihoods(panel)]

rule NGSadmix:
    input:
        beagle = lambda wildcards: config["ngs_admix"]["panels"][wildcards.panel]["path"]
    output:
        qopt = "{k}/{panel}_k{k}_{rep}.qopt"
    log:
        "logs/{k}_{panel}_{rep}.log"
    threads:
        8
    shell:
        """
        NGSadmix -likes {input.beagle} -K {wildcards.k} -P {threads} \
            -o {wildcards.k}/{wildcards.panel}_k{wildcards.k}_{wildcards.rep} \
             &> {log}
        """

rule parse_likelihood:
    input:
        log = lambda wildcards: expand_log(wildcards.panel, wildcards.k)#"{k}/{panel}_k{k}_{rep}.log"
    output:
        likelihood = "likelihoods_k{k}_{panel}.txt"
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


def search_likelihood_files(panel):
    nClusters = num_clusters(panel)
    myFiles = expand("likelihoods_k{k}_{panel}.txt", 
                k = nClusters, panel = panel)
    return(myFiles)

rule find_likelihood:
    input:
        lambda wildcards: search_likelihood_files(wildcards.panel)
    output:
        commands = "best_runs_scp_{panel}.txt"
    log:
        "logs/find_likelihood_{panel}.txt"
    shell:
        """
        Rscript ~/data/Git/Botocudos-scripts/Admixture/find_likelihood.R likelihoods \
            {wildcards.panel} {nRep} NGSadmix > {output.commands} 2>{log}
        """