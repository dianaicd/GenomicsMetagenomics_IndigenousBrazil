# Snakefile to downsample reads
import numpy as np 
from io import StringIO

configfile: "multiple_purposes.yaml"

samples = [sample for sample in list(config["downsample"]["groups"].keys())]
#number_reads = [3*pow(10, n) for n in range(2, 9)]
#depths = [0.01, 0.05, 0.1, 0.5, 1, 2, 9]
depths = [10,15]
wildcard_constraints:
    n = "\d+"

rule all:
    input:
        downsampled = expand("{depth}/{sample}_{depth}x.bam", 
                            sample = samples,
                            depth = depths),

        coverage = expand("{depth}/{sample}_{depth}x_depth.txt", 
                            sample = samples,
                            depth = depths)

rule idxstats:
    input:
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai"
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
        """
        cat {input.bam} | python ~/data/Git/Botocudos-scripts/DataQuality/read_length.py -o {output.length} 
        """

rule check_coverage:
    input:
        bam = "{file}.bam",
        idxstats = "{file}_idxstats.txt",
        length = "{file}.length"
    output:
        depth = "{file}_depth.txt"
    run:
        # Get genome length
        with open(input.idxstats, 'r') as file:
            genome_len = np.nansum([np.genfromtxt(StringIO(line))[1] for line in file.readlines()])
        # Number of mapped bases
        mapped_len = np.genfromtxt(input.length, delimiter = '\t')
        total_reads = np.nansum(mapped_len[:,1])
        total_bases = np.nansum(mapped_len[:,0]*total_reads)
        coverage = total_bases/genome_len

        # reads_to_sample = [(d, int(d*total_reads/coverage)) for d in depths] 
        np.savetxt(fname = output.depth, X = [coverage], fmt = "%.5f")

        # """

rule samtools_sample:
    input:
        bam = "{sample}.bam",
        idxstats = "{sample}_idxstats.txt",
        length = "{sample}.length"
    output:
        bam = "{depth}/{sample}_{depth}x.bam" 
    shell:
        """
        mkdir -p {wildcards.depth}
        numReads=$(python ~/data/Git/Botocudos-scripts/Downsample/calc_nReads.py \
          -i {input.idxstats} -l {input.length} -d {wildcards.depth})

        n=$( awk '{{sum += $3}} END {{print sum}}' {input.idxstats} )
        frac_reads=$(echo $numReads/$n |bc -l)

        samtools view -bs $frac_reads {input.bam} > {output.bam}

        """


# rule bam_chr:
#     input:
#         bam = "{sample}.bam"
#     output:
#         bam = temp("{sample}.{chr}.bam")
#     shell:
#         """
#         samtools view -b {input.bam} {wildcards.chr} > {output.bam}
        # """

# rule fred:
#     input:
#         bam = "{sample}.bam",
#         idxstats = "{sample}_idxstats.txt",
#         length = "{sample}.length"
#     output:
#         bam = "{depth}/{sample}_{depth}x.bam" 
#     shell:
#         """
#         mkdir -p {wildcards.depth}
#         numReads=$(python ~/data/Git/Botocudos-scripts/Downsample/calc_nReads.py \
#           -i {input.idxstats} -l {input.length} -d {wildcards.depth})

#         n=$( awk '{{sum += $3}} END {{print sum}}' {input.idxstats} )
#         python simulate_damage.py -i {input.bam} -b {output.bam} -ns $numReads -n $n
# rule shuf:
#     input:
#        # folder = "{depth}",
#         bam = "{sample}.{chr}.bam",
#         idxstats = "{sample}_idxstats.txt",
#         length = "{sample}.{chr}.length"
#     wildcard_constraints:
#         depth = "|".join([str(x) for x in depths])
#     output:
#         sam = temp("{depth}/{sample}_{depth}x_{chr}.sam")
#     shell:
#         """
#         mkdir -p {wildcards.depth}
#         numReads=$(python ~/data/Git/Botocudos-scripts/Downsample/calc_nReads.py \
#           -i {input.idxstats} -l {input.length} -d {wildcards.depth})

#         samtools view -H {input.bam} | \
#             sed "s/SM:{wildcards.sample}/SM:{wildcards.sample}_{wildcards.depth}/" \
#             > {output.sam} ;
#         samtools view {input.bam} | shuf -n $numReads >> {output.sam} 
#         """

# rule sort:
#     input:
#         sam = "{depth}/{sample}_{depth}x_{chr}.sam"
#     wildcard_constraints:
#         depth = "|".join([str(x) for x in depths])
#     output:
#         bam = "{depth}/{sample}_{depth}x_{chr}.bam"
#     params:
#         mem_sort = config["downsample"]["mem_sort"]
#     shell:
#         """
#         samtools view -b {input.sam} |samtools sort -m {params.mem_sort} > {output.bam}
#         """

# rule merge:
#     input:
#         bam = expand("{depth}/{sample}_{depth}x_{chr}.bam", chr = [str(x) for x in range(1,23)],
#                     depth = "{depth}", sample = "{sample}")
#     output:
#         bam = "{depth}/{sample}_{depth}x.bam"
#     shell:
#         """
#         samtools merge {output.bam} {input.bam}
#         """
