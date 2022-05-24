configfile: "bedtools_config.yaml"

input_files = config["input_files"]

rule all:
    input:
        genomecov = expand(
            "genomecov/{sample}.genomecov", 
            sample = list(input_files.keys()))

rule genomecov:
    input:
        bam = lambda wildcards: input_files[wildcards.sample]
    output:
        genomecov = "genomecov/{sample}.genomecov"
    resources:
        memory = 2*1024,
        runtime = 4*60
    shell:
        """
        bedtools genomecov -ibam {input.bam} > {output.genomecov}
        """