# Downsample a male and a female genome to
# different depths of coverage and see the ratio of correct sex assignations

# Snakefile to downsample reads
import numpy as np 
from io import StringIO
include: "parse_resources.smk"
configfile: "sex_config.yaml"


genomes = config.get("genome", "")
n_replicates = config.get("replicates", 10)
females = [female for genome in genomes for female in config.get("BAM_female")]
males = [male for genome in genomes for male in config.get("BAM_male")]
depths = [pow(10,n) for n in range(-7,0)]

wildcard_constraints:
    n = "\d+"

rule all:
    input:
        sex = expand("sex/{depth}/{sample}_{depth}x_rep{rep}.sex",
                            depth = depths,
                            sample = females + males,
                            rep = range(0,n_replicates)
                            )
        # downsampled = expand("{depth}/{sample}_{depth}x.bam", 
        #                     sample = females + males,
        #                     depth = depths),

        # coverage = expand("{depth}/{sample}_{depth}x_depth.txt", 
        #                     sample = females + males,
        #                     depth = depths),

rule idxstats:
    input:
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai"
    output:
        "{sample}_idxstats.txt"
    log:
        "logs/idxstats/{sample}_idxstats.txt"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    log:
        "logs/index/{sample}.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("downsample_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("downsample_time", attempt, 8)
    shell:
        "samtools index {input}"

rule length:
    input:
        bam="{file}.bam",
        bai="{file}.bam.bai"
    output:
        length="{file}.length"
    log:
        "logs/length/{file}.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("downsample_mem", attempt, 1),
        runtime=lambda wildcards, attempt: get_runtime_alloc("downsample_time", attempt, 4)
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
    log:
        "logs/check_coverage/{file}.txt"
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
        bam = lambda wildcards: config["BAM_female"].get(wildcards.sample) or config["BAM_male"].get(wildcards.sample),
        idxstats = "{sample}_idxstats.txt",
        length = "{sample}.length"
    output:
        bam = temp("{depth}/{sample}_{depth}x_rep{rep}.bam") 
    resources:
        memory=8*1024,
        runtime=1*60
    log:
        "logs/samtools_sample/{sample}_{depth}x_rep{rep}.txt"
    shell:
        """
        mkdir -p {wildcards.depth}
        numReads=$(python3 ~/data/Git/Botocudos-scripts/Downsample/calc_nReads.py \
          -i {input.idxstats} -l {input.length} -d {wildcards.depth})

        n=$( awk '{{sum += $3}} END {{print sum}}' {input.idxstats} )
        frac_reads=$(echo $numReads/$n |bc -l)

        samtools view -bs $frac_reads {input.bam} > {output.bam}

        """

rule genomecov:
    input:
        bam = "{depth}/{sample}_{depth}x_rep{rep}.bam" 
    output:
        genomecov = temp("genomecov/{sample}_{depth}x_rep{rep}.genomecov")
    resources:
        memory = 2*1024,
        runtime = 1*60
    log:
        "logs/genomecov/{sample}_{depth}x_rep{rep}.txt"
    shell:
        """
        bedtools genomecov -ibam {input.bam} > {output.genomecov}
        """

rule assign_sex:
    input:
        genomecov = "genomecov/{sample}_{depth}x_rep{rep}.genomecov"
    output:
        sex="sex/{depth}/{sample}_{depth}x_rep{rep}.sex",
    params:
        sex_params = ""
    log:
        "logs/sex/{sample}_{depth}x_rep{rep}.txt"
    shell:
        """
        Rscript ~/mapache/workflow/scripts/sex_assignation.r \
            --genomecov={input.genomecov} \
            --out={output.sex} \
            {params.sex_params}
        """
