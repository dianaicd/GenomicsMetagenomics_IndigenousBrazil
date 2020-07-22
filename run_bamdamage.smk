configfile: "config.yaml"

def get_samples(line):
    line = str(line)
    line = line.rstrip('\n').rsplit(sep="\t")
    return(line[5])

with open(config['samples'], 'r') as samples_file:
    next(samples_file)
    samples = [get_samples(line) for line in samples_file.readlines()]

samples = set(samples)

rule all:
    input:
        expand("{sample}/bamdamage/quality_{q}/{ref}.dam_5prime.csv", sample=samples, 
            q=config['mapping_quality'], ref=config['refs'])

# In case there is no .bai this is useful:
# Get temporal bam & bai of certain quality
rule get_bam_bai:
    input:
        "../{sample}/bam_merge/{sample}.{ref}.bam"
    output:
        bam=temp("{sample}/bamdamage/bams/{sample}.{ref}.mapq{q}.bam"),
        bai=temp("{sample}/bamdamage/bams/{sample}.{ref}.mapq{q}.bai")
    params:
        q= lambda wildcards: wildcards.q
    log:
        "logs/{sample}/bamdamage/bams/{sample}.{ref}.mapq{q}.bam_bai.log"
    shell:
        '''
        module add UHTS/Analysis/samtools/1.4;
        samtools view -b -q {params.q} {input} -o {output.bam} 2> {log};
        samtools index {output.bam} {output.bai} 2> {log}
        '''

rule bamdamage:
    """
    Run bamdamage to quantify the deamination pattern
    """
    input:
        # bam="../{sample}/bam_merge/{sample}.{ref}.bam",
        # bai="../{sample}/bam_merge/{sample}.{ref}.bai"
        bam="{sample}/bamdamage/bams/{sample}.{ref}.mapq{q}.bam",
        bai="{sample}/bamdamage/bams/{sample}.{ref}.mapq{q}.bai"

    output:
        dam="{sample}/bamdamage/quality_{q}/{ref}.dam.pdf",
        length="{sample}/bamdamage/quality_{q}/{ref}.length.pdf",
        length_table="{sample}/bamdamage/quality_{q}/{ref}.length.csv",
        dam_5prime_table="{sample}/bamdamage/quality_{q}/{ref}.dam_5prime.csv",
        dam_3prime_table="{sample}/bamdamage/quality_{q}/{ref}.dam_3prime.csv"      
    log:
        "logs/bamdamage/{sample}/bamdamage/quality_{q}/{ref}.log"
    threads: 1       
    message: "--- BAMDAMAGE {input.bam}"
    # params:
    #     q= lambda wildcards: wildcards.q
    shell:
        """
        module add UHTS/Analysis/samtools/1.4;
        module add R/3.5.1;
        ../scripts/bamdamage --output {output.dam} --output_length {output.length} {input.bam} 2> {log};
        """

        # ../scripts/bamdamage --mapquality {params.q} --output {output.dam} --output_length {output.length} {input.bam} 2> {log};