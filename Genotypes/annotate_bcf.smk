# annotate genotypes 
configfile: "genotype_calling.yaml"
include: "parse_resources.smk"
samples = list(config["geno_calls"]["bamlists"].keys())

chromosomes = [str(i) for i in range(1, 23)]

rule all:
    input: 
        bcf = expand("Masked/{sample}_masked_chr{chr}.bcf",
        chr = chromosomes,
        sample = samples
        )


rule annotate:
    input:
        annotation = "masks/chr{chr}.bed.gz",
        header = "masks/header.txt",
        bcf = "Raw/{sample}_chr{chr}.bcf"
    output:
        bcf = "Masked/{sample}_masked_chr{chr}.bcf"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("annotate_genos_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("annotate_genos_time", attempt, 4)
    log:
        "logs/annotate_{sample}_chr{chr}.log"
    shell:
        """
        bcftools annotate -a {input.annotation} \
            -m MASK=strict \
            -h {input.header} \
            -c CHROM,POS,FROM,TO \
             {input.bcf} | \
            bcftools view -i 'MASK=1' \
                -Ob -o {output.bcf} 2>{log}
        """
