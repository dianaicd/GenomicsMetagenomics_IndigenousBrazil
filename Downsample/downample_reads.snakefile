# Snakefile to downsample reads
import numpy as np 
from io import StringIO

configfile: "todownsample.yaml"

samples = [sample for sample in config["samples"].split(" ")]
number_reads = [3*pow(10, n) for n in range(2, 9)]
depths = [0.01, 0.05, 0.1, 0.5, 1, 2, 9]

wildcard_constraints:
    n = "\d+"

rule all:
    input:
        downsampled = expand("Downsampled/{depth}/{sample}_{depth}x.bam", 
                            sample = samples,
                            depth = depths)

rule idxstats:
    input:
        bam = "{sample}.bam",
        bai = "{sample}.bam"
    output:
        "{sample}_idxstats.txt"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule length:
    input:
        bam="{file}.bam",
        bai="{file}.bam.bai"
    output:
        length="{file}.length"
    shell:
        "cat {input.bam} |"
        "python ~/data/Git/Botocudos-scripts/DataQuality/read_length.py -o {output.length} ;"

# rule genome_coverage:
#     input:
#         bam = "{file}.bam",
#         idxstats = "{file}_idxstats.txt",
#         length = "{file}.length"
#     output:
#         downsample = "{file}_nreads.txt"
#     run:
#         # Get genome length
#         with open({input.idxstats}, 'r') as file:
#             genome_len = np.nansum([np.genfromtxt(StringIO(line))[1] for line in file.readlines()])
#         # Number of mapped bases
#         mapped_len = np.genfromtxt({input.length}, delimiter = '\t')
#         total_reads = np.nansum(mapped_len[:,1])
#         total_bases = np.nansum(mapped_len[:,0]*total_reads)
#         coverage = total_bases/genome_len

#         reads_to_sample = [(d, int(d*total_reads/coverage)) for d in depths] 
#         np.savetxt(fname = {output.downsample}, X = reads_to_sample, fmt = "%.5f %d")

rule shuf:
    input:
        bam = "{sample}.bam",
        idxstats = "{sample}_idxstats.txt",
        length = "{sample}.length"
    output:
        bam = "{folder}/{sample}_{depth}x.bam",
        sam = temp("{folder}/{depth}/{sample}_{depth}x.sam")
    shell:
        """
        numReads=$(python ~/data/Git/Botocudos-scripts/Downsample/calc_nReads.py \
          -i {input.idxstats} -l {input.length} -d {wildcards.depth})

        samtools view -H {input.bam} | \
            sed "s/SM:{wildcards.sample}/SM:{wildcards.sample}_{wildcards.depth}/" \
            > {output.sam} ;
        samtools view {input.bam} |shuf -n $numReads >> {output.sam} 
        samtools view -b {output.sam} |samtools sort > {output.bam}
        """