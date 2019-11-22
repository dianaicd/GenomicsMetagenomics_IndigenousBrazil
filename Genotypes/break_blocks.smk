import os,glob
# Break chromosomes in bam header into
# many small blocks of blockSize (bp)
# For simplicity, we take the chromosomes defined in the reference genome

def param_is_defined(name, default_par = False):
    myParameter = config[name] if name in config.keys() else default_par
    return(myParameter)
#-----------------------------------------------------------------------------#
# Variables and functions to begin with
#samples_names = [list(myPaths)[0] for mySamples in dict_bamlists.values() for myPaths in mySamples.values()]
#bamfile = list(dict_bamlists.keys())
# for bam in bamfile:
#     for value in dict_bamlists[bam]["paths"].values():

files = [myFile for bam in bamfile for myFile in list(dict_bamlists[bam]["paths"].values())]

blockSize = param_is_defined(name = "blockSize", default_par = 5e6)
ref_genome = param_is_defined(name = "ref_genome", 
                            default_par = "/scratch/axiom/FAC/FBM/DBC/amalaspi/popgen/reference_human/hs.build37.1/hs.build37.1.fa"
)


# Awful regular expresions to deal with files bearing bamfile and panel name
wildcard_constraints:
    bamfile =  "|".join([b for b in bamfile]) ,
    file = "|".join([f.replace(".bam", "") for f in files])

def expand_path(myBamfile):
    mySample = list(dict_bamlists[myBamfile]["paths"].keys())[0]
    path = dict_bamlists[myBamfile]["paths"][mySample]
    full_path = os.path.expanduser(path)
    return(full_path)
#-----------------------------------------------------------------------------#
# Get chromosome size and break it in blocks
chromosomes = [str(x) for x in range(1,23)] 
with open(ref_genome + ".fai", 'r') as index:
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

def expand_ranges(bamlist):
    outputs = ["{chr}\t{start}\t{end}".format(chr = chr,
    start = lower_chr[str(chr)][i], 
    end = upper_chr[str(chr)][i], bamlist = bamlist) for chr in chromosomes    for i in range(0, len(lower_chr[str(chr)]))]

    return(outputs)
#-----------------------------------------------------------------------------#
rule index_bam:
    input:
        bam = "{file}.bam"
    output:
        bai = "{file}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """
def match_file(bam):
    for f in files:
        if f.find(bam) > -1:
            return(os.path.expanduser(f))

rule get_header:
    input:
        bam = lambda wildcards: match_file(wildcards.bamfile),
        bai = lambda wildcards: match_file(wildcards.bamfile) + ".bai"

    output:
        header = "Raw/{bamfile}.header"
    shell:
        """
        samtools view -H {input.bam} > {output.header}
        """

# Break chromosomes into small blocks to
# call genotypes in parallel

rule get_positions:
    input:
        header = "Raw/{bamfile}.header"
    output:
        positions = "Raw/{bamfile}_positions.txt"
    run:
        with open(input.header, 'r') as header, open(output.positions, 'w') as positions:
            chr_size = {}
            for line in header.readlines():
                if line.split()[0] == "@SQ":
                    chr = line.split()[1].split(":")[1]
                    size = line.split()[2].split(":")[1].replace("\n", "")
                    chr_size[chr] = size

            lower_chr = []
            upper_chr = []

            chromosomes = list(chr_size.keys())
            for chr in chromosomes:
                chr = str(chr)
                lower_chr = [i for i in range(1, int(chr_size[chr]), int(float(blockSize)))]
                upper_chr = [i - 1 for i in lower_chr[1:]]
                upper_chr.append(chr_size[chr])
            
                [positions.write("\t".join([str(chr), str(lower_chr[i]), str(upper_chr[i])])+"\n")  
                for i in range(0, len(lower_chr))]
