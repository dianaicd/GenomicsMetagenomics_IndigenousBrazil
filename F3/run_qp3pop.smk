# Snakefile to run qp3pop
configfile: "multiple_purposes.yaml"

import itertools
#=============================================================================#
def read_pops(fileName):
    with open(fileName, 'r') as myFile:
        myPops = [line.replace("\n", "") for line in myFile.readlines() if line != "\n"]
    return(myPops)
#=============================================================================#


panels = [config["f3"]["panel"][d]["path"].split(".")[0] for d in list(config["f3"]["panel"].keys())]

panel = list(config["f3"]["panel"].keys())[0]
ind_pop = config["f3"]["ind_pop"]

include: "format_eigenstrat.smk"
# In qp3pop output, the populations are
# named as "SourceN"
source1 = config["f3"]["Source1"] if "Source1" in config["f3"].keys() else False
source2 = config["f3"]["Source2"] if "Source2" in config["f3"].keys() else False
outgroup = config["f3"]["Outgroup"] if "Outgroup" in config["f3"].keys() else "Yorubas"

rule all:
    input:
        pop_list = expand("{panel}.3pop.combinations", panel = panels),
        par = expand("{panel}_qp3pop.par", panel = panels),
        results = expand("{panel}.qp3test.results", panel = panels)



rule make_3pop_list:
    input:
        # source1 = source1,
        # source2 = source2,
        pops = "{panel}.mod.ind"
    output:
        pop_list = "{panel}.3pop.combinations"
    run:
        def read_source(source):
            if source:
                pops = read_pops(source)
            else:
                with open(input.pops) as myFile:
                    pops = [line.replace("\n", "").split()[2] for line in myFile.readlines()]
                    pops = [pop for pop in set(pops)]
            return(pops)
        
        pops1 = read_source(source1)
        pops2 = read_source(source2)
        p1_p2 = [subset for subset in itertools.product(pops1,pops2) if subset[0] != subset[1] ]

        with open(output.pop_list, 'w') as outFile:
            [outFile.write("\t".join([combination[0], combination[1], outgroup])+"\n") 
            for combination in p1_p2 ]

rule make_par:
    input:
        geno = "{panel}.eigenstratgeno",
        snp = "{panel}.snp",
        ind = "{panel}.mod.ind",
        pop_list = "{panel}.3pop.combinations"
    output:
        par = "{panel}_qp3pop.par"
    shell:
        """
        parops=(genotypename snpname indivname popfilename)
        pargs=({input.geno} {input.snp} {input.ind} {input.pop_list})
        
        for i in $(seq 0 3)
        do
            echo "${{parops[$i]}}: ${{pargs[$i]}}"    
        done    >{output.par}
        """

rule run_qp3pop:
    input:
        par = "{panel}_qp3pop.par",
        geno = "{panel}.eigenstratgeno",
        snp = "{panel}.snp",
        ind = "{panel}.mod.ind",
        pop_list = "{panel}.3pop.combinations"
    output:
        results = "{panel}.qp3test.results"
    shell:
        """
        qp3Pop -p {input.par} > {output.results}
        """

