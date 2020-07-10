# Get stats from the mapping output of "snakemake_aDNA_mapping"

configfile: 'config.yaml'

wildcard_constraints:
    query = 'MN\d+',
    reference = '\w+(_\d+)?_\d',
    lib = 'L\d+'

rule all:
    input:
        expand("{query}/stats/quality_{q}/{query}_stats.txt", query="{query}", 
            q=config['mapping_quality']),
        "all_samples/input_reads_samples.txt"

# Get libraries per sample using the sample.txt file from Sam's mapping pipeline.
# Get the paths to the original fastq files.
samples_file = open(config['samples'], 'r')
# Read the header
line = samples_file.readline()
sample_libs = {}
fastqs = {}
while True:
    line = samples_file.readline()
    if line:
        line = line.rstrip('\n').rsplit(sep="\t")
        if sample_libs.get(line[5]):
            sample_libs[line[5]] = [sample_libs[line[5]]] + [line[3]]
            fastqs[line[5]] = [fastqs[line[5]]] + [line[1]]
        else:
            sample_libs[line[5]] = line[3]
            fastqs[line[5]] = line[1]
    else:
        break
samples_file.close()

# Count the number of input reads
rule count_reads_input:
    input:
        lambda wildcards: expand("{fastq}", fastq=fastqs[wildcards.query])
    output:
        expand("{query}/input_reads/counts_input_reads.txt", query="{query}")
    shell:
        '''
        zgrep -c ^ {input} > output_zgrep;
        nb_libs=$(wc -l output_zgrep | cut -f1 -d' ');

        if [[ "$nb_libs" -gt 1 ]]; then
            cut -f2 -d':' output_zgrep | awk '{{sum+=$1}} END {{print sum}}' > {output}
        else
            cat output_zgrep > {output}
        fi
        rm output_zgrep;
        '''

# Join the counts of input reads of all samples
rule join_reads_input:
    input:
        expand("{query}/input_reads/counts_input_reads.txt", query="{query}")
    output:
        "all_samples/input_reads_samples.txt"
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
rule merge_lib_counts_rds_b4_rmdup:
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

# Get the length of the reference genome
rule get_ref_length:
    input:
        "../reference/{reference}/{reference}.fasta.fai"
    output:
        "reference_length/{reference}_lgth.txt"
    shell:
        "cut -f2 {input} > {output}"

# Join the lengths of the reference genomes in a single file
rule join_ref_lengths:
    input:
        expand("reference_length/{reference}_lgth.txt", reference=config['refs'])
    output:
        "reference_length/refs_length.txt"
    shell:
        "cat {input} > {output}"

# Average read length of mapped reads according to a certain mapping quality
rule get_avg_rd_lgth:
    input:
        "../{query}/{query}.{reference}.bam"
    output:
        expand("{query}/average_read_length/quality_{q}/{reference}_avg_rd_lgth.txt", 
            query="{query}", q=config['mapping_quality'], reference="{reference}")
    params:
        q=config['mapping_quality'],
        d_script=config['path2_read_length']
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        samtools view -bq{params.q} {input} | python {params.d_script} \
        --out read_length_output.txt;
        
        if [[ ! -s read_length_output.txt ]]; then
            echo "0" > {output}
        else
            awk '{{n_reads+=$2; lgth_read+=($1*$2)}} END {{print lgth_read/n_reads}}' \
            read_length_output.txt > {output};
        fi

        rm read_length_output.txt    
        '''

# Join the average read lengths in a single file
rule join_avg_rd_lgth:
    input:
        expand("{query}/average_read_length/quality_{q}/{reference}_avg_rd_lgth.txt", 
            query="{query}", q=config['mapping_quality'], reference=config['refs'])
    output:
        expand("{query}/average_read_length/quality_{q}/avg_rd_lgth_all_refs.txt", query="{query}",
            q=config['mapping_quality'])
    shell:
        '''
        cat {input} > {output}
        '''
# Get genome coverage
rule get_genome_coverage:
    input:
        "../{query}/{query}.{reference}.bam"
    output:
        expand("{query}/genome_coverage/quality_{q}/{reference}_genome_cov.txt", 
            query="{query}", q=config['mapping_quality'], reference="{reference}")
    params:
        q=config['mapping_quality']
    shell:
        '''
        module load Bioinformatics/Software/vital-it;
        module add UHTS/Analysis/samtools/1.4;
        module add UHTS/Analysis/BEDTools/2.29.2;
        samtools view -bq{params.q} {input} | bedtools genomecov -ibam - | grep -P 'genome\\t0' | \
        awk '{{print 1-$5}}' > {output}
        '''

# Join the genome coverage
rule join_genome_coverage:
    input:
        expand("{query}/genome_coverage/quality_{q}/{reference}_genome_cov.txt", 
            query="{query}", q=config['mapping_quality'], reference=config['refs'])
    output:
        expand("{query}/genome_coverage/quality_{q}/genome_coverage_all_refs.txt", query="{query}",
            q=config['mapping_quality'])
    shell:
        '''
        cat {input} > {output}
        '''

# rule create_table:
#     input:
#         expand('{query}/{query}.NC_004295_1_final_reads.txt', query=queries)
#     output:
#         "tables/{reference}.txt"
#     shell:
#         "cat {input} > {output}"
