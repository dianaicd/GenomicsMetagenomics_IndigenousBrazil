# We will try to benchmark the command 'shuf' 
# And Frédéric's scripts to downsample data from BAM files
configfile: "todownsample.yaml"

samples = config["samples"]
number_reads = [3*pow(10, n) for n in range(2, 9)]
pipelines = ["shuf", "fred", "me"]

wildcard_constraints:
    n = "\d+"

rule all:
    input:
        downsampled = expand("{folder}/{sample}_{n}.bam", 
                            folder = pipelines, sample = samples,
                            n = number_reads)

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

rule mine:
    input:
        bam = "{sample}.bam",
        idxstats = "{sample}_idxstats.txt"
    output:
        bam = "me/{sample}_{n}.bam"
    benchmark:
        repeat("benchmarks/me_{sample}_{n}.benchmark.txt", 10)
    shell:
        "samtools view -b {input.bam} | "
        "python downsample_reads.py -b {input.bam} -n {wildcards.n} -o {output.bam} -i {input.idxstats}"

rule fred:
    input:
        bam = "{sample}.bam",
        idxstats = "{sample}_idxstats.txt"
    output:
        bam = "fred/{sample}_{n}.bam"
    benchmark:
        repeat("benchmarks/fred_{sample}_{n}.benchmark.txt", 10)
    shell:
        "n=$( awk '{{sum += $3}} END {{print sum}}' {input.idxstats} );"
        "python simulate_damage.py -i {input.bam} -b {output.bam} -ns {wildcards.n} -n $n"

rule shuf:
    input:
        bam = "{sample}.bam"
    output:
        bam = "shuf/{sample}_{n}.bam",
        sam = temp("shuf/{sample}_{n}.sam")
    benchmark:
        repeat("benchmarks/shuf_{sample}_{n}.benchmark.txt", 10)
    shell:
        "samtools view -H {input.bam} > {output.sam} ; "
        "samtools view {input.bam} |shuf -n {wildcards.n} >> {output.sam} ;"
        "samtools view -b {output.sam} |samtools sort > {output.bam}"