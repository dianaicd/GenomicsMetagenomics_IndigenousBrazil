# queries = ["MN00010", "MN1943"]
# refs = ["AJ717293_1", "AY083234_1", "DQ333427_1", "DQ357065_1", "FN669502_1", "HQ340602_1",
# "NC_000883_2", "NC_001540_1", "NC_004295_1"]

# configfile: 'config.yaml'

# Count number of reads after rmdup

rule count_nb_reads:
    input:
        "../{query}/{query}.{reference}.bam"
    output:
        "{query}/{query}.{reference}_nb_reads.txt"
    wildcard_constraints:
        query = 'MN\d+',
        reference = '\w+(_\d+)?_\d'
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools view -c {input} > {output}
        '''

# Number of reads mapped before rmdup
# rule count_unique_reads:
#     input:
#         "../MN00010/MN00010/library_bam_merge/{library}/{library}.{reference}.bam"
#     output:
#         "MN00010/mapped_before_rmdup/{library}.{reference}.bam"
#     params:
#         q=0
#     shell:
#         "samtools view -c -q {params.q} {input} > {output}"




# rule create_table:
#     input:
#         expand('{query}/{query}.NC_004295_1_nb_reads.txt', query=queries)
#     output:
#         "tables/{reference}.txt"
#     shell:
#         "cat {input} > {output}"
