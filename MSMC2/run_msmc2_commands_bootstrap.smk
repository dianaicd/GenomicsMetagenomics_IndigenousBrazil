configfile: "multiple_purposes.yaml"
localrules: combineCrossCoal

threads = config["msmc2"]["threads"]

inds = config["msmc2"]["inds"]
INDICES = {f"{inds[i-1]}": "{x},{y}".format(x=2*i-2,y=2*i-1) for i in range(1, len(inds)+1)}
pops = config["msmc2"]["populations"] if "populations" in config["msmc2"] else {}
INDICES.update({p: ",".join(INDICES[ind] for ind in pops[p]) for p in pops.keys()})
combinations = config["msmc2"]["combinations"] if "combinations" in config["msmc2"] else {}
input_dir = config["msmc2"]["input_msmc2_dir"]
prefix_dir = config["msmc2"]["prefix_dir"]
bootstrap_prefix = config["msmc2"]["bootstrap_prefix"]

PARAMS = {
    # "A": "1*2+20*1+1*10",
    "B": "1*1+20*1+1*5",
    "C": "1*1+25*1+1*10"
    # "D": "1*1+25*1+1*5",
    # "E": "1*2+15*1+1*10",
    # "F": "1*2+15*1+1*5",
    # "G": "1*2+20*1+1*2",
    # "H": "1*2+1*3+20*1+1*5",
    # "I": "1*2+1*5+15*1+1*2",
    # "J": "1*1+1*3+15*1+1*5",
    # "original": "1*2+25*1+1*2+1*3"
}

popNames = list(INDICES.keys())
paramNames = list(PARAMS.keys())

wildcard_constraints:
    popName = "|".join(popNames)

def get_indices_crosspop(combination):
    [pop1,pop2] = combination.split("_")
    n_1 = len(pops[pop1])
    n_2 = len(pops[pop2])
    index_pop1 = [ i for i in range( 0, (2*n_1) ) ]
    index_pop2 = [ i for i in range( 2*n_1, 2*(n_2+n_1) ) ]
    indices=",".join([f"{i1}-{i2}" for i1 in index_pop1 for i2 in index_pop2])
    return(indices)

INDICES_crosspop ={
    f"{comb}":get_indices_crosspop(comb) for comb in combinations
}


# run it once per pop
[f"msmc2 -t {threads} -I {INDICES[popName]} -o {popName}_4hap input_files_msmc2/chr*.txt" for popName in INDICES.keys()]


bootstrap_dir = "bootstrapped"
reps = [str(i) for i in range(1, 21)]
chromosomes = [str(i) for i in range(1,23)]

rule all:
    input:
        # onepop = expand("onepop/rep{rep}_{popName}_param{p}_4hap.final.txt", 
        #                 popName = popNames,
        #                 p = paramNames,
        #                 rep = reps
        #                  )#,
        # crosspop = expand("crosspop/rep{rep}_{combinations}_param{p}_crosspop.final.txt",
        #                     p = paramNames,
        #                     combinations = combinations,
        #                     rep = reps
        # ),
        crosspop = expand("crosspop/rep{rep}.comb_{combinations}_param{p}_crosspop.final.txt",
            #"crosspop/rep{rep}.comb.{combinations}_param{p}_crosspop.final.txt",
                            p = paramNames,
                            combinations = combinations,
                            prefix_dir = prefix_dir, rep = reps
        ),
        combined = expand(
            "combined/rep{rep}.comb_{combinations}_param{p}.final.txt",
            p = paramNames,
            combinations = combinations, prefix_dir = prefix_dir, rep = reps
        )

rule msmc2_onepop:
    input:
        lambda wildcards: [f"{bootstrap_dir}/pop.{wildcards.popName}_{wildcards.rep}/bootstrap_multihetsep.chr{chr}.txt" 
         for chr in chromosomes]
    output:
        "onepop/rep{rep}_{popName}_param{p}_4hap.final.txt"
    threads:
        int(threads)
    params:
        index = lambda wildcards: INDICES[wildcards.popName],
        paramTime = lambda wildcards: PARAMS[wildcards.p],
        basename =  lambda wildcards: f"onepop/rep{wildcards.rep}_{wildcards.popName}_param{wildcards.p}_4hap",
        basename_input = lambda wildcards: f"{bootstrap_dir}/pop.{wildcards.popName}_{wildcards.rep}/bootstrap_multihetsep.chr*.txt" 
    resources:
        mem = 1024 * 32,
        runtime = 60 * 2
    shell:
        """
        ./msmc2 -p {params.paramTime} -t {threads} -o {params.basename} {params.basename_input}
        """


rule msmc2_crosspop:
    input:
        lambda wildcards: [f"{bootstrap_dir}/comb.{wildcards.pop1}_{wildcards.pop2}_{wildcards.rep}/bootstrap_multihetsep.chr{chr}.txt" 
                            for chr in chromosomes]
    output:
        "crosspop/rep{rep}.comb_{pop1}_{pop2}_param{p}_crosspop.final.txt"
    threads:
        int(threads)
    params:
        index = lambda wildcards: INDICES_crosspop[f"{wildcards.pop1}_{wildcards.pop2}"],
        basename = "crosspop/rep{rep}.comb_{pop1}_{pop2}_param{p}_crosspop",
        basename_input = lambda wildcards: f"{bootstrap_dir}/comb.{wildcards.pop1}_{wildcards.pop2}_{wildcards.rep}/bootstrap_multihetsep.chr*.txt" ,
        paramTime = lambda wildcards: PARAMS[wildcards.p],
    resources:
        mem = 1024 * 64,
        runtime = 60 * 12
    shell:
        """
        ./msmc2 -s -p {params.paramTime} -t {threads} -I {params.index} -o {params.basename} {params.basename_input}
        """


rule combineCrossCoal:
    input:
        pop1 = "onepop/rep{rep}_{pop1}_param{p}_4hap.final.txt",
        pop2 = "onepop/rep{rep}_{pop2}_param{p}_4hap.final.txt",
        crosspop = "crosspop/rep{rep}.comb_{pop1}_{pop2}_param{p}_crosspop.final.txt"
    output:
        "combined/rep{rep}.comb_{pop1}_{pop2}_param{p}.final.txt"
    params:
    resources:
    shell:
        """
        python combineCrossCoal.py \
            {input.crosspop} \
            {input.pop1} \
            {input.pop2} \
            > {output}      
        """
