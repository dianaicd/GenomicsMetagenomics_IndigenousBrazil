# Map using Bowtie 2

import pandas as pd

# localrules: create_dirs_index
# ruleorder: create_dirs_index > index

configfile: 'config.yaml'

# def create_dicts(line):
#     line = str(line)
#     line = line.rstrip('\n').rsplit(sep="\t")
#     if sample_libs.get(line[5]):
#             sample_libs[line[5]] = [sample_libs[line[5]]] + [line[3]]
#             fastqs[line[5]] = [fastqs[line[5]]] + [line[1]]
#     else:
#             sample_libs[line[5]] = line[3]
#             fastqs[line[5]] = line[1]
#     return(sample_libs, fastqs)

# def get_fastqFiles(sample_matrix, sample, library):
    # fastqFiles = glob.glob(config['samples'][sample][library][fileID])
    # fastqFile = sample_matrix[(sample_matrix["SM"]==sample) & (db["LB"]==library)][col].values
    # return(fastqFiles)

# def create_ref_basename(ref):
#     basename_path = "reference/" + str(ref) + "/" + str(ref)
#     return(basename_path)

# Save samples file into a matrix
samples_matrix = pd.read_csv(config['samples'], sep="\s+")
# Create a dictionary (samples as keys) of dictionaries (libraries as keys) with the fastqs paths
master_dic={sm:{lb:samples_matrix[(samples_matrix["SM"] == sm) & (samples_matrix["LB"] == lb)]["Data"].values[0] for lb 
                in samples_matrix[(samples_matrix["SM"] == sm)]["LB"].values} for sm in set(samples_matrix["SM"])}

samples = [*master_dic]
samples.sort()

# sample_libs = {}
# fastqs = {}
# with open(config['samples'], 'r') as samples_file:
#     samples_file.readline()
#     # sample_libs, fastqs are global
#     [create_dicts(line) for line in samples_file.readlines()]

# samples = [*sample_libs]
# samples.sort()

rule all:
    input:
        expand("bowtie2_mapping/{sample}/q_{q}/{sample}.{ref}.q_{q}.bam", sample=samples, q=config['mapping_quality'], 
        ref=config['refs'])

# rule create_dirs_index:
#     # input:
#     output:
#         # directory("reference/{ref}")
#         index_dirs = directory(expand("reference/{ref}", ref=config['refs'])),
#         check = "reference/index_dirs_created.done"
#     shell:
#         '''
#         mkdir {output.index_dirs}
#         touch {output.check}
#         '''

rule index:
    input:
        # directory(rules.create_dirs_index.output),
        # "reference/index_dirs_created.done",
        fasta = "reference/fasta/{ref}.fasta"
    output:
        multiext("reference/{ref}/{ref}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")#, 
        #ref=config['refs'])
        #reference/{ref}/{ref}
    log:
        "logs/reference/{ref}.log"
    params:
        runtime = "2:00:00"
    resources:
        memory = 4000
    shell:
        '''
        # module load Bioinformatics/Software/vital-it;
        # module add UHTS/Aligner/bowtie2/2.3.4.1;

        if [[ ! -d reference/{wildcards.ref} ]]; then
            mkdir -p reference/{wildcards.ref};
        fi

        bowtie2-build {input.fasta} reference/{wildcards.ref}/{wildcards.ref} 2> {log}
        '''

rule run_bowtie2:
    input:
        #index = "reference/{ref}/{ref}",
        #index = "reference/{ref}/{ref}.1.bt2",
        #index = expand("reference/{ref}/{ref}{ext}", ext=[".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]),
        multiext("reference/{ref}/{ref}", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        # index = 
        fastq = lambda wildcards: master_dic[wildcards.sample][wildcards.lib]
    output:
        temp("bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.sam")
    log:
        "logs/bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.log"
    threads: 8
    params:
        runtime = '10:00:00',
        mode = config['bowtie2_params'] #,
        # basename_ref = "{wildcards.ref}/{wildcards.ref}"
        # basename_ref = lambda wildcards : create_ref_basename(wildcards.ref)
    resources:
        memory = 2000
    shell:
        '''
        module add UHTS/Aligner/bowtie2/2.3.4.1;
        bowtie2 -x reference/{wildcards.ref}/{wildcards.ref} -U {input.fastq} {params.mode} --no-unal -S {output}
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
        memory = 2000
    shell:
        '''
        module add UHTS/Analysis/samtools/1.4;
        samtools view -bq{params.q} {input}
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
        memory = 2000
    shell:
        '''
        module add UHTS/Analysis/samtools/1.4;
        samtools sort -o {output} {input}
        '''

rule merge_libraries:
    input:
        lambda wildcards: expand("bowtie2_mapping/{sample}/{lib}/{sample}.{lib}.{ref}.q_{q}.sort.bam", sample="{sample}", 
        lib=master_dic[wildcards.sample].keys(), ref="{ref}", q="{q}")
    output:
        "bowtie2_mapping/{sample}/q_{q}/{sample}.{ref}.q_{q}.bam"
    log:
        "logs/bowtie2_mapping/{sample}/{sample}.{ref}.q_{q}.log"
    threads: 8
    params:
        runtime = '10:00:00'
    resources:
        memory = 2000
    shell:
        '''
        module add UHTS/Analysis/samtools/1.4;
        samtools merge {output} {input}
        '''




