configfile: "multiple_purposes.yaml"
import itertools
include: "parse_resources.smk"
def param_default(myPar, myDict, myDef):
    if myPar in myDict.keys():
        param = myDict[myPar]
    else:
        param = myDef
    return(param)

anc_fa = param_default("ancestral", config["relatedness"], "chimpHg19.fa")
minMapQ = param_default("minMapQ", config["relatedness"], 30)
minQ = param_default("minQ", config["relatedness"], 20)
sites = [s for s in list(config["relatedness"]["Sites"].keys())]
nThreads = param_default("nThreads", config["relatedness"], 8)
ref = param_default("ref_genome", config, "/scratch/axiom/FAC/FBM/DBC/amalaspi/popgen/reference_human/hs.build37.1/hs.build37.1.fa")
inds = [ind for ind in list(config["relatedness"]["Samples"].keys())]

def combine_inds(inds):
    myCombs = list(itertools.combinations(inds, 2))
    output = ["_".join(list(c))  + ".2D.SFS" for c in myCombs]
    return(output)
#=============================================================================#

rule all:
    input:
        SFS_2d = combine_inds(inds)

# Restrict analysis to chromosomes in ancestral sample
rule get_regions:
    input:
        fai = anc_fa + ".fai"
    output:
        regions = anc_fa + ".regions"
    shell:
        """
        cat {input.fai} | awk '{{print $1":1-"$2}}' > {output.regions}
        """
# get site allele frequency likelihood
rule make_saf:
    input:
        bam = lambda wildcards: config["relatedness"]["Samples"][wildcards.ind],
        regions = anc_fa + ".regions"
       # sites = lambda wildcards: config["relatedness"]["Sites"][wildcards.sites]
    output:
        idx = "{ind}.saf.idx"
    params:
        ancestral = anc_fa
    threads:
        get_threads("make_saf_threads", nThreads)
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("make_saf_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("make_saf_time", attempt, 24)
    log:
        "logs/{ind}_make_saf.log"
    shell:
        """
        angsd -gl 1 -anc {params.ancestral} -dosaf 1 -baq 1 -C 50 \
        -minMapQ {minMapQ} -minQ {minQ} -i {input.bam} \
        -rf {input.regions} \
         -out {wildcards.ind} -P {threads}  -ref {ref} \
         -checkBamHeaders 0 &>{log}
        """

# optimization of the .saf file which will estimate the SFS
rule optim_sfs:
    input:
        idx1 = "{ind1}.saf.idx",
        idx2 = "{ind2}.saf.idx"
    output:
        SFS_2d = "{ind1}_{ind2}.2D.SFS"
    threads:
        nThreads
    log:
        "logs/{ind1}_{ind2}_optim_sfs.log"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("optim_sfs_mem", attempt, 10),
        runtime=lambda wildcards, attempt: get_runtime_alloc("optim_sfs_time", attempt, 4)
    shell:
        """ 
        realSFS {input.idx1} {input.idx2} -P {threads} > {output.SFS_2d} 2>{log}
        """