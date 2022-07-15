configfile: "multiple_purposes.yaml"
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
nThreads = param_default("nThreads", config["relatedness"], 8)
ref = param_default("ref_genome", config, "/scratch/axiom/FAC/FBM/DBC/amalaspi/popgen/reference_human/hs.build37.1/hs.build37.1.fa")

#=============================================================================#

rule all:
    input:
        # saf = expand("saf/group_{n}_pop{i}.saf.idx", i = [1,2], n = [n for n in range(1, 256)]),
        # prior = expand("prior/group_{n}.ml", n = [n for n in range(1, 256)]),
        fst = expand("fst/group_{n}.fst.txt", n = [n for n in range(1, 256)])

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
        bamlist = "partitions/group_{n}_pop{i}.txt",
        regions = anc_fa + ".regions"
       # sites = lambda wildcards: config["relatedness"]["Sites"][wildcards.sites]
    output:
        idx = temp("saf/group_{n}_pop{i}.saf.idx"),
        saf = temp("saf/group_{n}_pop{i}.saf.gz"),
        pos = temp("saf/group_{n}_pop{i}.saf.pos.gz"),
    params:
        ancestral = anc_fa,
        basename = "saf/group_{n}_pop{i}"
    threads:
        4
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("make_saf_mem", attempt, 16),
        runtime=60*24 #lambda wildcards, attempt: get_runtime_alloc("make_saf_time", attempt, 24)
    log:
        "logs/group_{n}_pop{i}.make_saf.log"
    shell:
        """
        angsd -b {input.bamlist} \
        -anc {params.ancestral} -out {params.basename} -dosaf 1 -gl 1 \
        -baq 1 -C 50 \
        -minMapQ {minMapQ} -minQ {minQ} \
        -rf {input.regions} \
        -P {threads} \
        -ref hs.build37.1.fa \
        -noTrans \
        -checkBamHeaders 0 &>{log}
        """

#this is with 2pops
#first calculate per pop saf for each populatoin
#angsd -b list1  -anc hg19ancNoChr.fa -out pop1 -dosaf 1 -gl 1
#angsd -b list2  -anc hg19ancNoChr.fa -out pop2 -dosaf 1 -gl 1

# #calculate the 2dsfs prior
rule sfs2d_prior:
    input:
        pop1 = "saf/group_{n}_pop1.saf.idx",
        pop2 = "saf/group_{n}_pop2.saf.idx",
        saf_1 = temp("saf/group_{n}_pop1.saf.gz"),
        pos_1 = temp("saf/group_{n}_pop1.saf.pos.gz"),
        saf_2 = temp("saf/group_{n}_pop2.saf.gz"),
        pos_2 = temp("saf/group_{n}_pop2.saf.pos.gz"),
    output:
        original = "prior/group_{n}.ml",
        one_line = "prior/group_{n}_new.ml"
    resources:
        runtime = 24*60
    shell:
        """
        realSFS {input.pop1} {input.pop2} -nSites 1000000 > {output.original}
        awk '{{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}}END{{print}}' {output.original} > {output.one_line}
        """


rule index_sfs:
    input:
        pop1 = "saf/group_{n}_pop1.saf.idx",
        pop2 = "saf/group_{n}_pop2.saf.idx",
        saf_1 = temp("saf/group_{n}_pop1.saf.gz"),
        pos_1 = temp("saf/group_{n}_pop1.saf.pos.gz"),
        saf_2 = temp("saf/group_{n}_pop2.saf.gz"),
        pos_2 = temp("saf/group_{n}_pop2.saf.pos.gz"),
        sfs = "prior/group_{n}_new.ml"
    output:
        index = temp("index/group_{n}.fst.idx"),
        fst = temp("index/group_{n}.fst.gz"),
    resources:
        memory = 1024*32,
        runtime = 2*60
    params:
        basename = "index/group_{n}"
    shell:
        """
        realSFS fst index {input.pop1} {input.pop2} -sfs {input.sfs} -fstout {params.basename}
        """

rule fst:
    input:
        index = "index/group_{n}.fst.idx",
        fst = temp("index/group_{n}.fst.gz")
    output:
        fst = "fst/group_{n}.fst.txt"
    resources:
        runtime= 2*60
    shell:
        """
        realSFS fst stats {input.index} > {output.fst}
        """

# ../misc/realSFS pop1.saf.idx pop2.saf.idx -nSites 10000 >pop1.pop2.ml
# #prepare the fst for easy window analysis etc
# ../misc/realSFS fst index pop1.saf.idx pop2.saf.idx -sfs pop1.pop2.ml -fstout here
# #get the global estimate
# ../misc/realSFS fst stats here.fst.idx 
# -> FST.Unweight:0.069395 Fst.Weight:0.042349
# #below is not tested that much, but seems to work
# ../misc/realSFS fst stats2 here.fst.idx -win 50000 -step 10000 >slidingwindow