configfile: "to_trim.yaml"
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

trim_bp = config["Trim"]["trim_bp"]

rule all:
    input:
        expand("{folder}/{sample}_trim{trim}.bam", 
                folder = ["me", "victor", "sn"], 
                sample = [path_bam[path].replace(".bam", "") 
                            for bamlist in myDict.keys() 
                            for path in expand_path(bamlist)],
                trim = trim_bp )

rule trim_reads_bam:
    input:
        bam = lambda wildcards: bam_path[wildcards.sample+".bam"]
    output:
        bam = "me/{sample}_trim{trim}.bam"
    benchmark:
        repeat("benchmarks/me_{sample}_{trim}.benchmark.txt", 100)
    shell:
        """
        python trim_reads_bam.py --bam_in {input.bam} --bam_out {output.bam} --trim {wildcards.trim}
        """

rule trim_snake:
    input:
        bam = lambda wildcards: bam_path[wildcards.sample+".bam"]
    output:
        bam = "sn/{sample}_trim{trim}.bam"
    benchmark:
        repeat("benchmarks/sn_{sample}_{trim}.benchmark.txt", 100)
    run:
        untrimmed_bam = pysam.AlignmentFile(input.bam, "rb")
        trimmed_bam = pysam.AlignmentFile(output.bam, "wb", template = untrimmed_bam)
        trim = int(wildcards.trim)
        for read in untrimmed_bam:
            read.qual = '"'*trim + ''.join(list(read.qual)[trim:-trim]) + '"'*trim
            trimmed_bam.write(read)

        untrimmed_bam.close()
        trimmed_bam.close()


rule EndQTrimPipe:
    input:
        bam = lambda wildcards: bam_path[wildcards.sample + ".bam"]
    output:
        bam = "victor/{sample}_trim{trim}.bam"
    benchmark:
        repeat("benchmarks/victor_{sample}_{trim}.benchmark.txt", 100)
    shell:
        """
        ./EndQTrimPipe.pl -i {input.bam} -trim {wildcards.trim} -bed myBed.bed >{output.bam}
        """
    