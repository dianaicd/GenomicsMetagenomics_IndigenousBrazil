configfile: "multiple_purposes.yaml"
#configfile: "test.yml"
threads = config["msmc2"]["threads"]

inds = config["msmc2"]["inds"]
INDICES = {f"{inds[i-1]}": "{x},{y}".format(x=2*i-2,y=2*i-1) for i in range(1, len(inds)+1)}
pops = config["msmc2"]["populations"] if "populations" in config["msmc2"] else {}
INDICES.update({p: ",".join(INDICES[ind] for ind in pops[p]) for p in pops.keys()})
combinations = config["msmc2"]["combinations"] if "combinations" in config["msmc2"] else {}
input_dir = config["msmc2"]["input_msmc2_dir"]
prefix_dir = config["msmc2"]["prefix_dir"]

PARAMS = {
    # "A": "1*2+20*1+1*10",
    "B": "1*1+20*1+1*5",
    "C": "1*1+25*1+1*10",
    # "D": "1*1+25*1+1*5",
    # "E": "1*2+15*1+1*10",
    # "F": "1*2+15*1+1*5",
    # "G": "1*2+20*1+1*2",
    # "H": "1*2+1*3+20*1+1*5",
    # "I": "1*2+1*5+15*1+1*2",
    # "J": "1*1+1*3+15*1+1*5",
    "original": "1*2+25*1+1*2+1*3"
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
    #index_pop1 = INDICES[pop1].split(",")
    #index_pop2 = INDICES[pop2].split(",")
    indices=",".join([f"{i1}-{i2}" for i1 in index_pop1 for i2 in index_pop2])
    return(indices)

INDICES_crosspop ={
    f"{comb}":get_indices_crosspop(comb) for comb in combinations
}
# run it once per pop
[f"msmc2 -t {threads} -I {INDICES[popName]} -o {popName}_4hap input_files_msmc2/chr*.txt" for popName in INDICES.keys()]




rule all:
    input:
        # onepop = expand("{prefix_dir}/onepop/{popName}_param{p}_4hap.final.txt", 
        #                 popName = popNames,
        #                 p = paramNames,
        #                 prefix_dir = prefix_dir
        #                 ),
        crosspop = expand("{prefix_dir}/crosspop/comb.{combinations}_param{p}_crosspop.final.txt",
                            p = paramNames,
                            combinations = combinations,
                            prefix_dir = prefix_dir
        ),
        combined = expand(
            "{prefix_dir}/combined/combined_{combinations}_param{p}.final.txt",
            p = paramNames,
            combinations = combinations, prefix_dir = prefix_dir
        )



rule msmc2_crosspop:
    input:
        lambda wildcards: [f"{input_dir}/comb.{wildcards.pop1}_{wildcards.pop2}.chr{chr}.txt" for chr in [str(i) for i in range(1,23)]]
    output:
        "{prefix_dir}/crosspop/comb.{pop1}_{pop2}_param{p}_crosspop.final.txt"
    threads:
        int(threads)
    params:
        index = lambda wildcards: INDICES_crosspop[f"{wildcards.pop1}_{wildcards.pop2}"],
        basename = "{prefix_dir}/crosspop/comb.{pop1}_{pop2}_param{p}_crosspop",
        paramTime = lambda wildcards: PARAMS[wildcards.p],
        input_data = lambda wildcards: f"{input_dir}/comb.{wildcards.pop1}_{wildcards.pop2}.chr*.txt"
    resources:
        mem = 1024 * 120,
        runtime = 60 * 6
    shell:
        """
        ./msmc2 -p {params.paramTime} -t {threads} \
            -I {params.index} -o {params.basename} \
            {params.input_data}
        """

rule msmc2_onepop:
    input:
        lambda wildcards: [f"{input_dir}/pop.{wildcards.popName}.chr{chr}.txt" for chr in [str(i) for i in range(1,23)]]
    output:
        "{prefix_dir}/onepop/pop.{popName}_param{p}_4hap.final.txt"
    threads:
        int(threads)
    params:
        #index = lambda wildcards: INDICES[wildcards.popName],
        paramTime = lambda wildcards: PARAMS[wildcards.p],
        basename =  lambda wildcards: f"{prefix_dir}/onepop/pop.{wildcards.popName}_param{wildcards.p}_4hap",
        input_data = lambda wildcards: f"{input_dir}/pop.{wildcards.popName}.chr*.txt"
    resources:
        mem = 1024 * 32,
        runtime = 60 * 2
    shell:
        """
        ./msmc2 -p {params.paramTime} -t {threads}  -o {params.basename} {params.input_data}
        """


rule combineCrossCoal:
    input:
        pop1 = "{prefix_dir}/onepop/pop.{pop1}_param{p}_4hap.final.txt",
        pop2 = "{prefix_dir}/onepop/pop.{pop2}_param{p}_4hap.final.txt",
        crosspop = "{prefix_dir}/crosspop/comb.{pop1}_{pop2}_param{p}_crosspop.final.txt"
    output:
        "{prefix_dir}/combined/combined_{pop1}_{pop2}_param{p}.final.txt"
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
# python combineCrossCoal.py SuruiB_Bot crosspop.final.txt SuruiB_4hap.final.txt Botocudo_4hap.final.txt > combined_SuruiB_Botocudo.final.txt
# msmc2 -t $threads -I $Karitiana -o Karitiana_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Botocudo -o Botocudo_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Maya -o Maya_4hap multihetsep/chr*.txt
# msmc2 -t $threads -I $Surui -o Surui_4hap multihetsep/chr*.txt