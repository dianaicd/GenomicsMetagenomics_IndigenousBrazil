# Create counts table from the output of "map_bowtie2.smk"
import pandas as pd

# Same config as "map_bowtie2.smk" 
configfile: 'config.yaml'

# Get samples using the sample.txt file from Sam's mapping pipeline.
samples_matrix = pd.read_csv(config['samples'], sep="\s+")
samples = list(set(samples_matrix["SM"]))
samples.sort()

rule all:
    input:
        expand("mapping_counts/mapped_reads_q{q}.tsv", q=config['mapping_quality'])

rule count:
    input:
        "{sample}/q_{q}/{sample}.{virus}.q_{q}.bam"
    output:
        temp("{sample}/q_{q}/{sample}.{virus}.q_{q}.txt")
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools view -c {input} > {output}
        '''

rule merge_sample_counts:
    input:
        expand("{sample}/q_{q}/{sample}.{virus}.q_{q}.txt", sample="{sample}", q="{q}", 
        virus=config['refs'])
    output:
        temp("{sample}/q_{q}/{sample}.q_{q}.counts.txt")
    shell:
        '''
        cat {input} > {output}
        '''

rule merge_by_quality:
    input:
        expand("{sample}/q_{q}/{sample}.q_{q}.counts.txt", sample=samples, q="{q}")
    output:
        "mapping_counts/mapped_reads_q{q}.tsv"
    shell:
        '''
        paste {input} > {output}
        '''
