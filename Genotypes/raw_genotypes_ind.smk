# configfile: "multiple_purposes.yaml"
configfile: "genotype_calling.yaml"
include: "parse_resources.smk"
import os.path
# Call genotypes on a medium/high coverage sample
# output all possible sites
#-----------------------------------------------------------------------------#
# Variables and functions to begin with
def param_is_defined(name, default_par = False):
    myParameter = config[name] if name in config.keys() else default_par
    return(myParameter)

def get_range(chr = False):
    if(chr):
        ranges = ["_".join([str(chr), 
                    str(lower_chr[str(chr)][i]), 
                    str(upper_chr[str(chr)][i])]) for i in range(0, len(lower_chr[str(chr)]))]
    else:
        ranges = ["_".join([str(chr), 
                    str(lower_chr[str(chr)][i]), 
                    str(upper_chr[str(chr)][i])]) for chr in chromosomes for i in range(0, len(lower_chr[str(chr)]))]
    return(ranges)
#-----------------------------------------------------------------------------#

# 
# bamlists = list(config["geno_calls"]["bamlists"].keys())
individuals = list(config["geno_calls"]["individuals"].keys())

ref_genome = param_is_defined(name = "ref_genome", 
                            default_par = "/scratch/axiom/FAC/FBM/DBC/amalaspi/popgen/reference_human/hs.build37.1/hs.build37.1.fa"
)

blockSize = param_is_defined(name = "blockSize", default_par = 5e6)

with open(ref_genome + ".fai", 'r') as index:
    chr_size = {}
    for line in index.readlines():
        chr,size = line.split()[0:2]
        chr_size[chr] = size
    chromosomes = list(chr_size.keys())
    chromosomes = [str(x) for x in range(1,23)] 
    lower_chr = {}
    upper_chr = {}
    for chr in chromosomes:
        chr = str(chr)
        lower_chr[chr] = [i for i in range(1, int(chr_size[chr]), int(float(blockSize)))]
        upper_chr[chr] = [i - 1 for i in lower_chr[chr][1:]]
        upper_chr[chr].append(chr_size[chr])

chromosomes = [str(i) for i in range(1, 23)]

wildcard_constraints:
    chr = "(" + "|".join([i for i in chromosomes]) + ")"
#-----------------------------------------------------------------------------#
# Get chromosome size and break it in blocks
# include: "break_blocks.smk"
# include: "make_bamlist.smk"
#-----------------------------------------------------------------------------#

rule all:
    input:
        expand(
            "Raw/{ind}_chr{chr}.bcf",
            ind = individuals,
            chr = chromosomes
            )


rule print_positions:
    output:
        positions = "Raw/{ind}_positions.txt"
    run:
        with open(output.positions, "w") as out:
            [out.write(line.replace("_", "\t")+"\n") for chr in chromosomes for line in get_range(chr)]

rule call_genos:
    input:
        ind = lambda wildcards: os.path.expanduser(config["geno_calls"]["individuals"][wildcards.ind])
    output:
        raw_genos = temp("Raw/chr{chr}/{ind}_{chr}_{start}_{end}.bcf")
    log:
        "logs/call_{ind}_{chr}_{start}_{end}.log"
    params:
        minMapQ=30,
        minBaseQ=20,
        ref =param_is_defined(name = "ref_genome"),
        threads = 4,
        moreno2019 = param_is_defined("Genos_Moreno2019", "Yes")
    resources:
        memory=2*1024,
        runtime= 60
    shell:
        """
            if [ {params.moreno2019} == "Yes" ]
            then
               samtools mpileup -q {params.minMapQ} \
                    -t DP -C50 -uf {params.ref} \
                    --region {wildcards.chr}:{wildcards.start}-{wildcards.end} \
                 {input.ind} | bcftools call -f GQ -c \
                    --threads {params.threads} \
                    -Ob -o {output.raw_genos} - \
                2>{log}
            else
                bcftools mpileup -C 50 -q {params.minMapQ} -Q {params.minBaseQ} -a FMT/DP,SP \
                        -f {params.ref} --threads {params.threads} \
                        --regions {wildcards.chr}:{wildcards.start}-{wildcards.end} \
                        -Ou  {input.ind} | \
                    bcftools annotate -c RPB | \
                    bcftools call --threads {params.threads} -c -V indels \
                        -Ob -o {output.raw_genos} - \
                2>{log}
            fi
        """

rule concat_genos_chr:
    input:
        positions = "Raw/{bamfile}_positions.txt",
        bcf = lambda wildcards: expand("Raw/chr{chr}/{bamfile}_{range}.bcf", 
                    bamfile = "{bamfile}",
                    range =  get_range(wildcards.chr),
                    chr = "{chr}")
    wildcard_constraints:
        chr = "|".join([str(x) for x in chromosomes])
    output:
        genos_list = "Raw/{bamfile}_{chr}.list",
        genos_all_depths = "Raw/{bamfile}_chr{chr}.bcf"
    params:  
        threads = 4
    log:
        "logs/concat_{bamfile}_{chr}.log"
    run:
        def write_line(line, myChr):
            chr,start,end = line.replace("\n", "").split()
            myNewLine = "Raw/chr{chr}/{bamfile}_{chr}_{start}_{end}.bcf\n".format(
                bamfile = wildcards.bamfile, 
                chr = chr, start = start, 
                end = end
                )

            if chr == myChr:
                genos_list.write(myNewLine)

        with open(input.positions, 'r') as positions, open(output.genos_list, 'w') as genos_list:
            [write_line(line, wildcards.chr) for line in positions.readlines()]

        myCommand = f"bcftools concat -f {output.genos_list} -Ob -o {output.genos_all_depths} --threads {params.threads}"
        os.system(myCommand)
