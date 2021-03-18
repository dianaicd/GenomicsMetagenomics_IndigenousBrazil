# pipeline for fastngsadmix

# The idea is to run this for a few panels and a few values of K for many independent individuals
configfile: "multiple_purposes.yaml"

panel_dict = config["fastNGSadmix"]["panels"]
panels = list(panel_dict.keys())

min_K = {panel:panel_dict[panel]["min_K"] for panel in panels}
max_K = {panel:panel_dict[panel]["max_K"] for panel in panels}
replicates = {panel:panel_dict[panel]["replicates"] for panel in panels}
individuals = { panel:list(panel_dict[ panel ][ "individuals" ].keys()) for panel in panels }
# print( individuals )  

include: "run_admixture.snakefile"

wildcard_constraints:
    panel =  "(" + "|".join([p for p in panels])+ ")"#(?!_" + "|_".join([ind for panel in panels for ind in individuals[panel]]) + ")" ,
    # ind = ")".join(["(?<!"+p for p in panels])  + "_)" + "|".join([ind for panel in panels for ind in individuals[panel]])


# Making a reference panel in NGSadmix to use with fastNGSadmix
# ./ngsadmix32 -likes $panel.beagle -K $K -printInfo 1 -outfiles $panel.beagle
rule all_fastNGSadmix:
    input:
        qopt = ["{K}/{panel}_{ind}_k{K}.qopt".format( 
                                                      panel = panel,
                                                      K = K,
                                                      ind = ind)
                    for panel in panels
                    for ind in individuals[ panel ]
                    for K in range( min_K[ panel ], max_K[ panel ] )
        ]
        
rule gunzip:
    input: 
        "{file}.gz"
    output:
        "{file}"
    wildcard_constraints:
        file = ".*.(fopt)"
    resources:
        runtime = lambda wildcards, attempt: get_runtime_alloc( "gunzip_time", 
                                                                attempt, 
                                                                1
                                                                ),
        memory = lambda wildcards, attempt: get_memory_alloc("gunzip_mem", attempt, 1)
    shell:
        "gunzip {input}"

# although we have filter files for every K, they are virtually the
# same across the different values of K as they depend
# on the number of sites that were effectively analysed (i.e., after applying a maf filter)
rule find_overlap:
    input:
        panel = lambda wildcards: panel_dict[wildcards.panel]["path"],
        filter = lambda wildcards: "{K}/{panel}_k{K}.filter".format(panel = wildcards.panel,
                                                                    K = min_K[ wildcards.panel ])
    output:
        positions = temp("positions_{panel}.txt"),
        overlaping = temp("overlaping_{panel}.txt")  
    resources:
        runtime = lambda wildcards, attempt: get_runtime_alloc( "overlap_time", attempt, 1),
        memory = lambda wildcards, attempt: get_memory_alloc("overlap_mem", attempt, 1)   
    log:
        "logs/{panel}.log"  
    shell:
        """
        ## find overlap of .beagle and .filter files and creates tmp.ref with 6 first column of refPanel_
        cut -f1,2,3 {input.panel} > {output.positions}
        Rscript example/makeRefNGSadmix.R {output.positions} {input.filter} 
        mv {output.positions}.tmp.ref {output.overlaping}
        """

# We can then create a reference panel from the .fopt.gz file, with n populations:
rule create_refPanel:
    input:
        panel = lambda wildcards: panel_dict[wildcards.panel]["path"], 
        fopt = "{K}/{panel}_k{K}.fopt",
        filter = "{K}/{panel}_k{K}.filter",
        positions = "positions_{panel}.txt",
        overlaping = "overlaping_{panel}.txt"
    output:
        refPanel = "{K}/refPanel_{panel}_{K}.txt"
    resources:
        runtime = lambda wildcards, attempt: get_runtime_alloc( "NGSAdmix_time", attempt, 1),
        memory = lambda wildcards, attempt: get_memory_alloc("refPanel_mem", attempt, 1)
    log:
        "logs/{K}_{panel}.log"
    shell:
        """
        k_seq=( $( for i in $( seq 1 {wildcards.K} ) ; do echo "K$i" ; done ) )
        echo "id chr pos name A0_freq A1 ${{k_seq[@]}}" > {output.refPanel}
        paste {input.overlaping}  {input.fopt} >> {output.refPanel}
        """

# And then we can create a file with the number of individuals in each reference population,
# summing each column of the .qopt file:
# rule create_nInd:
#     input:
#         qopt = "{K}/{panel}_k{K}.qopt"
#     output:
#         nInd = "{K}/nInd_{panel}_k{K}.txt"
#     shell:
#         """
#         k_seq=( $( for i in $( seq 1 {wildcards.K} ) ; do echo "K$i" ; done ) )

#         sum_cols=( $( for i in $( seq 1 {wildcards.K} ) 
#                     do  
#                         cut -f $i -d " " {input.qopt} | paste -sd+ | bc 
#                     done 
#                     ) 
#                 )
#         echo "${{k_seq[@]}}" > {output.nInd}
#         echo ${{sum_cols[@]}} >> {output.nInd}
#         """

# Then we can run fastNGSadmix using this reference panel
rule run_fastNGSadmix:
    input:
        geno_like = lambda wildcards: panel_dict[ wildcards.panel ][ "individuals" ][wildcards.ind],
        refPanel = "{K}/refPanel_{panel}_{K}.txt",
        nInd = "nInd_{panel}.txt"
    output:
        qopt = "{K}/{panel}_{ind}_k{K}.qopt"
    log:
        "logs/{K}_{panel}_{ind}.log"
    resources:
        runtime = lambda wildcards, attempt: get_runtime_alloc( "fastNGSadmix_time", 
                                                                attempt, 
                                                                1
                                                                ),
        memory = lambda wildcards, attempt: get_memory_alloc("fastNGSadmix_mem", attempt, 4)
    shell:
        """
        fastNGSadmix -likes {input.geno_like} \
            -fname {input.refPanel} \
            -Nname {input.nInd} \
            -out {wildcards.K}/{wildcards.panel}_{wildcards.ind}_k{wildcards.K}  \
            -boot 100 \
            -whichPops all
        """
