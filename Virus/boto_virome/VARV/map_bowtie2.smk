# Mapping using Bowtie 2

import pandas as pd

# To specify reference genomes, mapping quality, and Bowtie2 parameters.
configfile: 'config.yaml'

# Save samples file into a matrix
samples_matrix = pd.read_csv(config['samples'], sep="\s+")
# Create a dictionary (samples as keys) of dictionaries (libraries as keys) with the fastqs paths
master_dic={sm:{lb:samples_matrix[(samples_matrix["SM"] == sm) & 
            (samples_matrix["LB"] == lb)]["Data"].values[0] for lb in 
            samples_matrix[(samples_matrix["SM"] == sm)]["LB"].values} for sm in 
            set(samples_matrix["SM"])}
# Get sample ID
samples = [*master_dic]
samples.sort()

rule all:
    input:
        expand("bowtie2_mapping/{sample}/q_{q}/{sample}.{ref}.q_{q}.bam", sample=samples, 
        q=config['mapping_quality'], ref=config['refs'])

rule index:
    input:
        fasta = "reference/fasta/{ref}.fasta"
    output:
        multiext("reference/{ref}/{ref}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", 
        ".rev.2.bt2")
    log:
        "logs/reference/{ref}.log"
    params:
        runtime = "2:00:00"
    resources:
        memory = 4000
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Aligner/bowtie2/2.3.4.1;
        bowtie2-build {input} reference/{wildcards.ref}/{wildcards.ref} 2> {log}
        '''

rule run_bowtie2:
    input:
        multiext("reference/{ref}/{ref}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", 
        ".rev.2.bt2"),
        fastq = lambda wildcards: master_dic[wildcards.sample][wildcards.lib]
    output:
        temp("bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.sam")
    log:
        "logs/bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.log"
    threads: 8
    params:
        runtime = '10:00:00',
        mode = config['bowtie2_params']
    resources:
        memory = 5000
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Aligner/bowtie2/2.3.4.1;
        bowtie2 -x reference/{wildcards.ref}/{wildcards.ref} -U {input.fastq} {params.mode} \
        --no-unal -S {output}
        '''

rule select_quality:
    input:
        "bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.sam"
    output:
        temp("bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.q_{q,\d+}.bam")
    log:
        "logs/bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.q_{q}.bam.log"
    threads: 8
    params:
        runtime = '10:00:00',
        q = lambda wildcards: wildcards.q
    resources:
        memory = 4000
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools view -bq{params.q} {input} > {output}
        '''

rule sort:
    input:
        "bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.q_{q}.bam"
    output:
        temp("bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.q_{q,\d+}.sort.bam")
    log:
        "logs/bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.q_{q}.sort.log"
    threads: 8
    params: 
        runtime = '10:00:00'
    resources:
        memory = 4000
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools sort -o {output} {input}
        '''

rule merge_libraries:
    input:
        lambda wildcards: 
        expand("bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.q_{q}.sort.bam", 
        sample="{sample}", lib=master_dic[wildcards.sample].keys(), ref="{ref}", q="{q}")
    output:
        "bowtie2_mapping/{sample}/q_{q}/{sample}.{ref}.q_{q}.bam"
    log:
        "logs/bowtie2_mapping/{sample}/{sample}.{ref}.q_{q}.log"
    threads: 8
    params:
        runtime = '10:00:00'
    resources:
        memory = 4000
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools merge {output} {input}
        '''




