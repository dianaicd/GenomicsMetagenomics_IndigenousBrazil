# Script to run abba baba using angsd

# Have a magic file with paths to BAM files to use, and their population name
configfile: "multiple_purposes.yaml"

def param_or_default(key, default):
    param = config["abbababa2"][key] if key in config["abbababa2"].keys() else default
    return param

# As the error estimation must have been run at this point,
# the following files must be available in fasta format

# Genome used to identify ancestral states for error correction
outgroup_err_name = [key for key,value in config["abbababa2"]["Outgroup_file_error"].items()][0]
outgroup_err_fasta = config["abbababa2"]["Outgroup_file_error"][outgroup_err_name] 
perfect_name = [key for key,value in config["abbababa2"]["Perfect_file"].items()][0]
perfect_fasta = config["abbababa2"]["Perfect_file"][perfect_name]

# Outgroup population
ancestral_name = [key for key,value in config["abbababa2"]["Outgroup_file_ancestral"].items()][0]
ancestral_fasta = config["abbababa2"]["Outgroup_file_ancestral"][ancestral_name] 

# a long bamlist will have the bams listed per group
long_bamlist = [key for key,value in config["abbababa2"]["bamlists"].items()][0]

groups = [key for key,value in config["abbababa2"]["bamlists"][long_bamlist]["paths"].items()]


chromosomes = [i for i in range(1,23)]

wildcard_constraints:
    chr =   "|".join([str(i) for i in chromosomes])

baseQ       = param_or_default(key = "BaseQuality", default = 30)
mapQ        = param_or_default(key = "MapQuality", default = 30)
threads     = param_or_default(key = "abbababa_threads", default = 1)
blockSize   = param_or_default(key = "BlockSize", default = 5e6)

# You must have an error estimation per population
# defined in your initial sets
# Therefore, first run: "error_angsd.smk"
localrules: all,merge_bamlists,files_size_name_ancErr_per_pop,merge_abbababa2
rule all:
    input:
        errorCorr = f"{long_bamlist}/{long_bamlist}.ErrorCorr.txt",
        D = f"{long_bamlist}/{long_bamlist}.D"
        
# You need a file with population size and population name
# If everything went ok while runnning error_angsd.smk,
# you have a bamlist per group that you can merge.
# As every time I run this pipeline the groups are obtained
# from a dictionary (config), the order of the groups may change.
# Thus, I save the order in a file.
rule merge_bamlists:
    input: 
        expand("{group}/{group}.txt", group = groups)
    output:
        long_bamlist = "{long_bamlist}.txt",
        group_order = "{long_bamlist}.order"
    wildcard_constraints:
        long_bamlist = long_bamlist
    shell:
        """
        cat {input} > {output.long_bamlist}
        for group in {groups} ; do echo $group  ; done > {output.group_order}
        """

rule files_size_name_ancErr_per_pop:
    input:
        bamlist = "{long_bamlist}.txt",
        group_order = "{long_bamlist}.order",

    output:
        pop_size = "{long_bamlist}.popsize",
        pop_name = "{long_bamlist}.popname",
        anc_error = "{long_bamlist}.ancErr"
    run:
        with open(output.pop_size, 'w') as pop_size, open(output.pop_name, 'w') as pop_name, open(output.anc_error, 'w') as anc_error:
            with open(input.group_order, "r") as group_list:
                for line in group_list.readlines():
                    group = line.replace("\n", "")
                    group_name = group
                    if group in config["abbababa2"]["correct_error"]:
                        group_error = f"{group}/{group}_perfect.{perfect_name}_outgroup.{outgroup_err_name}_ancErr.ancError"
                    else:
                        group_error = "NA"

                    with open(f"{group}/{group}.txt", 'r') as file:
                        group_size = len(file.readlines())

                    pop_size.write(str(group_size) + "\n")
                    pop_name.write(group_name + "\n")
                    anc_error.write(group_error + "\n")
            # Don't forget to append the ancestral sample
            pop_size.write("1")
            pop_name.write(ancestral_name)
            anc_error.write("NA")


# Get chromosome size and break it in blocks
with open(perfect_fasta + ".fai", 'r') as index:
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
        ancestral = ancestral_fasta,
        pop_size = "{bamlist}.popsize",
        # region = "{chr}"
    output:
        abbababa = temp("{bamlist}/{chr}/{start}_{end}.abbababa2")
    threads:
        threads
    resources:
        runtime = 60,
        mem = 2*1024
    params:
        mapQ = mapQ,
        baseQ = baseQ,
        prefix = "{bamlist}/{chr}/{start}_{end}"
    shell:
        """
        angsd -doAbbababa2 1 -bam {input.bamlist} \
        -sizeFile {input.pop_size} -doCounts 1 \
        -out {params.prefix} \
        -anc {input.ancestral} \
        -r {wildcards.chr}:{wildcards.start}-{wildcards.end} \
        -useLast 0 \
        -checkBamHeaders 0 \
        -minQ {params.baseQ} -minMapQ {params.mapQ} -p {threads}
        """

def expand_abbababa2(bamlist):
    outputs = [
        "{bamlist}/{chr}/{start}_{end}.abbababa2".format(
        chr = chr,
        start = lower_chr[str(chr)][i], 
        end = upper_chr[str(chr)][i],
        bamlist = bamlist
        ) 
        for chr in chromosomes    for i in range(0, len(lower_chr[str(chr)]))
        ]

    return(outputs)

rule merge_abbababa2:
    input:
        lambda wildcards: expand_abbababa2(wildcards.bamlist)   #expand("{chr}/{chr}_{bamlist}.abbababa2", chr = chromosomes, bamlist = "{bamlist}")
    output:
        "{bamlist}/{bamlist}.abbababa2"
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
        abbababa = "{bamlist}/{bamlist}.abbababa2",
        pop_size = "{bamlist}.popsize",
        pop_name = "{bamlist}.popname",
        ancError = "{bamlist}.ancErr"    
    output:
        errorCorr = "{bamlist}/{bamlist}.ErrorCorr.txt",
        ErrorCorrTransRem = "{bamlist}/{bamlist}.ErrorCorr.TransRem.txt",
        TransRem = "{bamlist}/{bamlist}.TransRem.txt"
    params:
        basename="{bamlist}/{bamlist}"
    shell:
        """
        Rscript /software/UHTS/Analysis/ANGSD/0.931/R/estAvgError.R \
        angsdFile={params.basename} \
        out={params.basename} \
        sizeFile={input.pop_size} \
        errFile={input.ancError} \
        nameFile={input.pop_name}
        """

rule get_D_values:
    input:
        abbababa = "{bamlist}/{bamlist}.abbababa2",
        errorCorr = "{bamlist}/{bamlist}.ErrorCorr.txt",
    output:
        D = "{bamlist}/{bamlist}.D"
    shell:
        """
        cut -f4-6 {input.abbababa} > {input.abbababa}.mini
        Rscript ~/data/Git/Botocudos-scripts/Dstat/D_values_from_angsd.R \
            {input.abbababa}.mini {input.errorCorr} {output.D}
        rm {input.abbababa}.mini
        """