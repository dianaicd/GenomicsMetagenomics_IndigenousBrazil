# Get stats from the mapping output of "snakemake_aDNA_mapping"

configfile: 'config.yaml'

wildcard_constraints:
    query = 'MN\d+',
    reference = '\w+(_\d+)?_\d',
    lib = 'L\d+'

# Get libraries per sample using the sample.txt file from Sam's mapping pipeline
samples_file = open(config['samples'], 'r')
# Read the header
line = samples_file.readline()
sample_libs = {}
while True:
    line = samples_file.readline()
    if line:
        line = line.rstrip('\n').rsplit(sep="\t")
        if sample_libs.get(line[5]):
            sample_libs[line[5]] = [sample_libs[line[5]]] + [line[3]]
        else:
            sample_libs[line[5]] = line[3]
    else:
        break
samples_file.close()

# Count number of reads after rmdup
rule count_final_reads:
    input:
        "../{query}/{query}.{reference}.bam"
    output:
        expand("{query}/count_final_reads/quality_{q}/{query}.{reference}_final_reads.txt",
        query="{query}", q=config['mapping_quality'], reference="{reference}")
    params:
        q=config['mapping_quality']
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools view -c -q {params.q} {input} > {output}
        '''
# Join the final reads counts in a single file
rule join_counts_final_reads:
    input:
        expand('{query}/count_final_reads/quality_{q}/{query}.{reference}_final_reads.txt',
        query="{query}", q=config['mapping_quality'], reference= config['refs'])
    output:
        expand("{query}/count_final_reads/quality_{q}/counts_final_rds_all_refs.txt",
        query="{query}", q=config['mapping_quality'])
    shell:
        "cat {input} > {output}"

# Number of reads mapped before rmdup (per library)
rule count_rds_b4_rmdup:
    input:
        "../{query}/{lib}/library_bam_merge/{lib}.{reference}.bam"
    output:
        expand("{query}/mapped_before_rmdup/quality_{q}/per_lib/{lib}.{reference}_counts.txt",
        query="{query}", q=config['mapping_quality'], lib="{lib}", reference="{reference}")
    params:
        q=config['mapping_quality']
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools view -c -q {params.q} {input} > {output}
        '''

# Sum the counts before rmdup according to sample and reference
rule merge_counts_rds_b4_rmdup_lib:
    input:
        lambda wildcards:
        expand("{query}/mapped_before_rmdup/quality_{q}/per_lib/{lib}.{reference}_counts.txt",
        query="{query}", q=config['mapping_quality'], lib=sample_libs[wildcards.query],
        reference=wildcards.reference)
    output:
        expand("{query}/mapped_before_rmdup/quality_{q}/merge_counts_lib/{query}.{reference}"
            "_counts_b4_rmdup.txt", query="{query}", q=config['mapping_quality'],
            reference="{reference}")
    shell:
        "cat {input} | awk '{{sum+=$1}} END {{print sum}}' > {output}"

# Join counts of reads before rmdup
rule join_counts_rds_b4_rmdup:
    input:
        expand("{query}/mapped_before_rmdup/quality_{q}/merge_counts_lib/{query}.{reference}"
            "_counts_b4_rmdup.txt", query="{query}", q=config['mapping_quality'],
            reference=config['refs'])
    output:
        expand("{query}/mapped_before_rmdup/quality_{q}/counts_b4_rmdup_all_refs.txt",
            query="{query}", q=config['mapping_quality'])
    shell:
        "cat {input} > {output}"

# Get the length of the reference genome
rule get_ref_length:
    input:
        "../reference/{reference}/{reference}.fasta.fai"
    output:
        "reference_length/{reference}_lgth.txt"
    shell:
        "cut -f2 {input} > {output}"

rule join_ref_lengths:
    input:
        expand("reference_length/{reference}_lgth.txt", reference=config['refs'])
    output:
        "reference_length/refs_length.txt"
    shell:
        "cat {input} > {output}"

# rule create_table:
#     input:
#         expand('{query}/{query}.NC_004295_1_final_reads.txt', query=queries)
#     output:
#         "tables/{reference}.txt"
#     shell:
#         "cat {input} > {output}"
