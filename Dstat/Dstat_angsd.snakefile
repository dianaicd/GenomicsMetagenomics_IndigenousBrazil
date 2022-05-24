# Script to run abba baba using angsd

# Have a magic file with paths to BAM files to use, and their population name
configfile: "multiple_purposes.yaml"

def param_or_default(key, default):
    param = config[key] if key in config.keys() else default
    return param

ancestral_fasta = param_or_default(key = "Outgroup_file", default="Yoruba.fa")
ancestral_fasta = ancestral_fasta.split(".")[0]
perfect_fasta = param_or_default(key = "Perfect_file", default="Yoruba.fa")
perfect_fasta = perfect_fasta.split(".")[0]

chromosomes = [i for i in range(1,23)]

baseQ       = param_or_default(key = "BaseQuality", default = 20)
mapQ        = param_or_default(key = "MapQuality", default = 30)
threads     = param_or_default(key = "abbababa_threads", default = 1)
blockSize   = param_or_default(key = "BlockSize", default = 5e6)

wildcard_constraints:
    bamlist = "|".join([bamlist for bamlist in config["bamlists"].keys()])

# You must have an error estimation per population
# defined in your initial sets
#include: "error_angsd.smk"

rule ALL:
    input:
        # ancError = set(["{group}/{group}_ancErr.ancError".format(bamlist = bamlist, group = group) 
        #         for bamlist in list(config["bamlists"].keys())
        #         for group in config["bamlists"][bamlist]["paths"]]),
        abbabba = expand("{bamlist}.abbababa2", 
                        bamlist = config["bamlists"].keys()),
        errorCorr = expand("{bamlist}.ErrorCorr.txt", bamlist = config["bamlists"].keys())
        
# You need a file with population size and population name
# If everything went ok while runnning error_angsd.smk,
# the bamlist has the populations in a certain order

# rule files_size_name_ancErr_per_pop:
#     input:
#         bamlist = "{bamlist}.txt",
#         bamlist_group = lambda wildcards: ["Bams_group/{group}.txt".format(group = group, bamlist = "{bamlist}") 
#                             # for bamlist in list(config["bamlists"].keys()) 
#                             for group in config["bamlists"][wildcards.bamlist]["paths"]],
#         ancError_group = lambda wildcards: ["{group}/{group}_ancErr.ancError".format(group = group) for group in config["bamlists"][wildcards.bamlist]["paths"]]

#     output:
#         pop_size = "{bamlist}.popsize",
#         pop_name = "{bamlist}.popname",
#         anc_error = "{bamlist}.ancErr"
#     run:
#         bamlist = input.bamlist.split(".")[0]
#         all_groups = list(config["bamlists"][bamlist]["paths"].keys())
#         with open(output.pop_size, 'w') as pop_size, open(output.pop_name, 'w') as pop_name, open(output.anc_error, 'w') as anc_error:
#         # This file has the actual order of the groups (populations)
#             for group in input.bamlist_group:
#                 group_name = group.split("/")[-1].split(".")[0]
#                 with open(group, 'r') as file:
#                     size = len(file.readlines())
                    
#                     pop_size.write(str(size) + "\n")
#                     pop_name.write(group_name + "\n")
#                     anc_error.write(group_name+"/"+group_name+"_ancErr.ancError\n")
#             # Don't forget to append the ancestral sample
#             pop_size.write("1")
#             pop_name.write(ancestral_fasta)
#             anc_error.write("NA")

# Get chromosome size and break it in blocks
with open(perfect_fasta + ".fa.fai", 'r') as index:
    chr_size = {}
    for line in index.readlines():
        chr,size = line.split()[0:2]
        chr_size[chr] = size
    lower_chr = {}
    upper_chr = {}
    for chr in chromosomes:
        chr = str(chr)
        lower_chr[chr] = [i for i in range(1, int(chr_size[chr]), int(float(blockSize)))]
        upper_chr[chr] = [i - 1 for i in lower_chr[chr][1:]]
        upper_chr[chr].append(chr_size[chr])

rule run_abbababa2:
    input:
        bamlist = "{bamlist}.txt",
        ancestral = ancestral_fasta + ".fa",
        pop_size = "{bamlist}.popsize",
        # region = "{chr}"
    output:
        abbababa = "{chr}/{chr}_{start}_{end}_{bamlist}.abbababa2"
    threads:
        4
    params:
        mapQ = mapQ,
        baseQ = baseQ
    resources:
        runtime = 60 * 12,
        mem = 1024 * 12
    shell:
        """
        angsd -doAbbababa2 1 -bam {input.bamlist} \
            -sizeFile {input.pop_size} -doCounts 1 \
            -out {wildcards.chr}/{wildcards.chr}_{wildcards.start}_{wildcards.end}_{wildcards.bamlist} \
            -anc {input.ancestral} \
            -r {wildcards.chr}:{wildcards.start}-{wildcards.end} \
            -useLast 0 \
            -checkBamHeaders 0 \
            -minQ {params.baseQ} -minMapQ {params.mapQ} -p {threads}
        """

def expand_abbababa2(bamlist):
    outputs = ["{chr}/{chr}_{start}_{end}_{bamlist}.abbababa2".format(chr = chr,
    start = lower_chr[str(chr)][i], 
    end = upper_chr[str(chr)][i], bamlist = bamlist) for chr in chromosomes    for i in range(0, len(lower_chr[str(chr)]))]

    return(outputs)

rule merge_abbababa2:
    input:
        lambda wildcards: expand_abbababa2(wildcards.bamlist)   #expand("{chr}/{chr}_{bamlist}.abbababa2", chr = chromosomes, bamlist = "{bamlist}")
    output:
        "{bamlist}.abbababa2"
    shell:
        """
        myFiles=({input})
        head -n1 ${{myFiles[0]}} >{output}
        for file in ${{myFiles[@]}}
        do
            tail -n+2 $file
        done >> {output}
        """

rule correct_error:
    input:
        abbababa = "{bamlist}.abbababa2",
        pop_size = "{bamlist}.popsize",
        pop_name = "{bamlist}.popname",
        ancError = "{bamlist}.ancErr"    
    output:
        errorCorr = "{bamlist}.ErrorCorr.txt",
        ErrorCorrTransRem = "{bamlist}.ErrorCorr.TransRem.txt",
        TransRem = "{bamlist}.TransRem.txt"
    shell:
        """
        Rscript DSTAT angsdFile={wildcards.bamlist} \
        out={wildcards.bamlist} \
        sizeFile={input.pop_size} \
        errFile={input.ancError} \
        nameFile={input.pop_name}
        """