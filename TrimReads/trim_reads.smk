configfile: "multiple_purposes.yaml"
import glob, os, pysam, getopt, copy, sys

myDict = config["Trim"]["bams"]

def expand_path(bamlist):
    paths = list(myDict[bamlist].values())
    full_paths = [os.path.expanduser(p) for p in paths]
    bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

bam_path = {path.split("/")[-1]:path for bamlist in myDict.keys() 
            for path in expand_path(bamlist)}
path_bam = {bam_path[key]:key for key in bam_path.keys()}

# trim_bp = config["Trim"]["trim_bp"]
trim_bp = [i for i in range(5,6)]

rule all:
    input:
        expand("{sample}_trim{trim}.bam",  
                sample = [path_bam[path].replace(".bam", "") 
                            for bamlist in myDict.keys() 
                            for path in expand_path(bamlist)],
                trim = trim_bp )

rule trim_reads:
    input:
        bam = lambda wildcards: bam_path[wildcards.sample+".bam"]
    output:
        bam = "{sample}_trim{trim}.bam"
    run:
        untrimmed_bam = pysam.AlignmentFile(input.bam, "rb")
        trimmed_bam = pysam.AlignmentFile(output.bam, "wb", template = untrimmed_bam)
        trim = int(wildcards.trim)
        for read in untrimmed_bam:
            read.qual = '"'*trim + ''.join(list(read.qual)[trim:-trim]) + '"'*trim
            trimmed_bam.write(read)

        untrimmed_bam.close()
        trimmed_bam.close()