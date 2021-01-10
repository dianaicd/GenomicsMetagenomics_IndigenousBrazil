# Snakemake to split the genome into chunks
# that will be analysed in GLIMPSE

configfile: "multiple_purposes.yaml"

# Split the genome into chunks
rule split_genome:
    input:
        vcf_ref_panel = "reference_panel/1000GP.chr{chr}.noNA12878.sites.vcf.gz"
    output:
        chunks = "chunks.chr{chr}.txt"
    params:
        window_size = config['GLIMPSE']['window_size'] if 'window_size' in config['GLIMPSE'].keys() else 2000000,
        buffer_size = config['GLIMPSE']['buffer_size'] if 'buffer_size' in config['GLIMPSE'].keys() else 200000
    shell:
        """

        GLIMPSE_chunk --input {input.vcf_ref_panel} \
        --region {chr} \
        --window-size {params.window_size} \
        --buffer-size {params.buffer_size} \
        --output {output.chunks}

        """
