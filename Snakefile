# Get stats from the mapping output of "snakemake_aDNA_mapping"

configfile: 'config.yaml'

wildcard_constraints:
    query = 'MN\d+',
    reference = '\w+(_\d+)?_\d',
    lib = 'L\d+'

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
samples = [*sample_libs]
samples_file.close()

rule all:
    input:
        expand("{query}/stats/quality_{q}/{query}_stats.txt", query=samples, 
            q=config['mapping_quality']),
        "all_samples/input_reads_samples.txt"

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
        expand("{query}/input_reads/counts_input_reads.txt", query=samples)
    output:
        "all_samples/input_reads_samples.txt"
    shell:
        "cat {input} > {output}"

# Number of reads mapped before rmdup (per library)
rule count_rds_b4_rmdup:
    input:
        "../{query}/{lib}/library_bam_merge/{lib}.{reference}.bam"
    output:
        "{query}/mapped_before_rmdup/quality_{q}/per_lib/{lib}.{reference}_counts.txt"
    params:
        q= lambda wildcards: wildcards.q
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
            query="{query}", q="{q}", lib=sample_libs[wildcards.query], 
            reference=wildcards.reference)
    output:
        "{query}/mapped_before_rmdup/quality_{q}/merge_counts_lib/{query}.{reference}"
        "_counts_b4_rmdup.txt"
    shell:
        "cat {input} | awk '{{sum+=$1}} END {{print sum}}' > {output}"

# Join counts of reads before rmdup
rule join_counts_rds_b4_rmdup:
    input:
        expand("{query}/mapped_before_rmdup/quality_{q}/merge_counts_lib/{query}.{reference}"
            "_counts_b4_rmdup.txt", query="{query}", q="{q}", reference=config['refs'])
    output:
        "{query}/mapped_before_rmdup/quality_{q}/counts_b4_rmdup_all_refs.txt"
    shell:
        "cat {input} > {output}"

# Count number of reads after rmdup
rule count_final_reads:
    input:
        "../{query}/{query}.{reference}.bam"
    output:
        "{query}/count_final_reads/quality_{q}/{query}.{reference}_final_reads.txt"
    params:
        q= lambda wildcards: wildcards.q
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
            query="{query}", q="{q}", reference= config['refs'])
    output:
        "{query}/count_final_reads/quality_{q}/counts_final_rds_all_refs.txt"
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
        "{query}/average_read_length/quality_{q}/{reference}_avg_rd_lgth.txt"
    params:
        q=lambda wildcards: wildcards.q,
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
            query="{query}", q="{q}", reference=config['refs'])
    output:
        "{query}/average_read_length/quality_{q}/avg_rd_lgth_all_refs.txt"
    shell:
        '''
        cat {input} > {output}
        '''
# Get genome coverage
rule get_genome_coverage:
    input:
        "../{query}/{query}.{reference}.bam"
    output:
        "{query}/genome_coverage/quality_{q}/{reference}_genome_cov.txt"
    params:
        q= lambda wildcards: wildcards.q
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
        expand("{query}/genome_coverage/quality_{q}/{reference}_genome_cov.txt", query="{query}", 
            q="{q}", reference=config['refs'])
    output:
        "{query}/genome_coverage/quality_{q}/genome_coverage_all_refs.txt"
    shell:
        '''
        cat {input} > {output}
        '''

# Get the deamination damage
rule get_deam_damage:
    input:
        three_prime="{query}/bamdamage/quality_{q}/{reference}.dam_3prime.csv",
        five_prime="{query}/bamdamage/quality_{q}/{reference}.dam_5prime.csv"
    output:
        damage="{query}/bamdamage/quality_{q}/{reference}_damage.txt",
        value_3prime=temp("{query}/bamdamage/quality_{q}/{reference}_3prime_value"),
        value_5prime=temp("{query}/bamdamage/quality_{q}/{reference}_5prime_value")
    shell:
        '''
        cut -f3 -d',' {input.three_prime} | sed -n 2p > {output.value_3prime}
        cut -f12 -d',' {input.five_prime} | sed -n 2p > {output.value_5prime}
        paste {output.value_3prime} {output.value_5prime} > {output.damage}
        '''

rule join_damage:
    input:
        expand("{query}/bamdamage/quality_{q}/{reference}_damage.txt", query="{query}", 
            q="{q}", reference=config['refs'])
    output:
        "{query}/bamdamage/quality_{q}/damage_all_refs.txt"
    shell:
        '''
        cat {input} > {output}
        '''

# Create table with statistics (per sample)
rule create_table:
    input:
        b4_rmdup="{query}/mapped_before_rmdup/quality_{q}/counts_b4_rmdup_all_refs.txt",
            # query="{query}", q="{q}",
        after_rmdup="{query}/count_final_reads/quality_{q}/counts_final_rds_all_refs.txt",
            # query="{query}", q="{q}",
        ref_lgth="reference_length/refs_length.txt",
        avg_rd_lgth="{query}/average_read_length/quality_{q}/avg_rd_lgth_all_refs.txt", 
            # query="{query}", q="{q}",
        g_cov="{query}/genome_coverage/quality_{q}/genome_coverage_all_refs.txt", 
            # query="{query}", q="{q}",
        damage="{query}/bamdamage/quality_{q}/damage_all_refs.txt"
    output:
        "{query}/stats/quality_{q}/{query}_stats.txt"
    shell:
        '''
        paste {input.b4_rmdup} {input.after_rmdup} {input.avg_rd_lgth} {input.g_cov} \
        {input.ref_lgth} {input.damage} > {output}
        '''

