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

inds = ["MN00010","MN00013","MN00016","MN00019","MN00021","MN00022",
"MN00023","MN00039","MN0003","MN00045","MN00056","MN00064","MN00066","MN00067","MN00068",
"MN00069","MN0008","MN0009","MN00118","MN00119","MN00316","MN00346","MN01701","MN1943", "9Botocudos", "22Botocudos"]
#=============================================================================#
localrules: get_regions
rule all:
    input:
        # saf = expand("saf/group_{n}_pop{i}.saf.idx", i = [1,2], n = [n for n in range(1, 256)]),
        # prior = expand("prior/group_{n}.ml", n = [n for n in range(1, 256)]),
        fst = expand("het/{ind}_new.ml", ind = inds),
        window = expand("thetas/{ind}.thetasWindow.gz.pestPG", ind = inds)

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
        #bam = "bams/{ind}.hg19.bam",
        bam = "{ind}.txt",
        regions = anc_fa + ".regions"
       # sites = lambda wildcards: config["relatedness"]["Sites"][wildcards.sites]
    output:
        idx = "saf/{ind}.saf.idx",
        saf = "saf/{ind}.saf.gz",
        pos = "saf/{ind}.saf.pos.gz",
    params:
        ancestral = anc_fa,
        basename = "saf/{ind}"
    threads:
        4
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("make_saf_mem", attempt, 32),
        runtime=60*24 #lambda wildcards, attempt: get_runtime_alloc("make_saf_time", attempt, 24)
    log:
        "logs/make_saf_{ind}.log"
    shell:
        """
        angsd -b {input.bam} \
        -anc {params.ancestral} -out {params.basename} -dosaf 1 -gl 1 \
        -baq 1 -C 50 \
        -minMapQ {minMapQ} -minQ {minQ} \
        -rf {input.regions} \
        -P {threads} \
        -ref hs.build37.1.fa \
        -noTrans 1 \
        -checkBamHeaders 0 &>{log}
        """

#this is with 2pops
#first calculate per pop saf for each populatoin
#angsd -b list1  -anc hg19ancNoChr.fa -out pop1 -dosaf 1 -gl 1
#angsd -b list2  -anc hg19ancNoChr.fa -out pop2 -dosaf 1 -gl 1

# #calculate het
rule het_estimation:
    input:
        saf = "saf/{ind}.saf.idx",
        saf_1 = "saf/{ind}.saf.gz",
        pos_1 = "saf/{ind}.saf.pos.gz",
    output:
        original = "het/{ind}.ml",
        one_line = "het/{ind}_new.ml"
    resources:
        runtime = 24*60,
        memory = 1024*8
    shell:
        """
        realSFS {input.saf}  -nSites 1000000 > {output.original}
        awk '{{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}}END{{print}}' {output.original} > {output.one_line}
        """

rule theta:
    input:
        saf = "saf/{ind}.saf.idx",
    output:
        sfs = "sfs/{ind}.sfs",
        thetasIdx = "thetas/{ind}.thetas.idx",
        window = "thetas/{ind}.thetasWindow.gz.pestPG"
    resources:
        runtime = 12*60,
        memory = 1024*8
    shell:
        """
        realSFS {input.saf} -nSites 1000000 > {output.sfs}
        mv {output.sfs} sfs/{wildcards.ind}_old.sfs
        awk '{{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}}END{{print}}' sfs/{wildcards.ind}_old.sfs > {output.sfs}

        realSFS saf2theta {input.saf} -outname thetas/{wildcards.ind} -sfs {output.sfs}

        #Estimate for every Chromosome/scaffold
        thetaStat do_stat {output.thetasIdx}

        #Do a sliding window analysis based on the output from the make_bed command.
        thetaStat do_stat {output.thetasIdx} -win 50000 -step 10000  -outnames thetas/{wildcards.ind}.thetasWindow.gz
        """