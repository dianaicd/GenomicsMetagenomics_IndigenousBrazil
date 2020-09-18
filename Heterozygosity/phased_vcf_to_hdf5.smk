# Snakefile to prepare 1000G haplotypes panel in hdf5 python objects
configfile: "multiple_purposes.yaml"
include: "parse_resources.smk"

n_haplotypes = 2504


rule all:
    input:
        # haplotypes = expand("Haplotypes/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_haplo{haplo}.txt.gz",
        #                     chr = [chr for chr in range(1,23)],
        #                     haplo = [1,2]
        #                     ),
        # metadata = expand("Haplotypes/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_metadata.txt.gz",
        #                     chr = [chr for chr in range(1,23)]
        #                     )
        hdf5 = expand( "hdf5/1000G.{chr}.hdf5", chr = [chr for chr in range(1,23)] )

# Prepare annotation with info relative to genetic distance in centimorgans
# the annotation file should be bgzipped and tabix-indexed
# include a header file with info for the columns
rule prepare_annotation:
    input:
        genetic_map = "HapmapII_GRCh37_RecombinationHotspots/genetic_map_GRCh37_chr{chr}.txt"
    output:
        annotation = "annotation/{chr}.txt.gz",
    shell:
        """
        # Remove first line and chr prefix
        tail -n+2 {input.genetic_map} | sed 's/^chr//' | bgzip > {output.annotation}
        tabix -s1 -b2 -e2 {output.annotation}
        """

# Extract di-allelic SNPs positions and annotate with genetic distances
rule extract_annotate_diallelic_VCFs:
    input:
        vcf = "VCF/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        annotation = "annotation/{chr}.txt.gz",
        header = "annotation/header.txt"
    output:
        vcf = "diallelic/{chr}.diallelic.vcf.gz"
    threads:
        8
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("extract_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("extract_time", attempt, 8),
    shell:
        """
         bcftools view -m 2 -M 2 -v snps -Ob  {input.vcf} | \
         bcftools annotate \
            --threads {threads} \
            -a {input.annotation} \
            -h {input.header} \
            -c CHROM,POS,MAP,MAP2 \
            -Oz -o {output.vcf} 

        """

rule convert_to_hdf5:
    input:
        vcf = "diallelic/{chr}.diallelic.vcf.gz"
    output:
        hdf5 = "hdf5/1000G.{chr}.hdf5"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("process_hdf5_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("process_hdf5_time", attempt, 8),
    run:
        import allel
        allel.vcf_to_hdf5(input.vcf, output.hdf5, fields = "*")

# rule get_haplotypes:
#     input:
#         vcf = "VCF/ALL.chr{chr}.phase3{prefix}.vcf.gz",
#         positions = "Haplotypes/ALL.chr{chr}.phase3{prefix}.diallelic.snps.positions"
#     output:
#         haplotype = "Haplotypes/ALL.chr{chr}.phase3{prefix}_haplo{haplo}.txt.gz"
#     params:
#         n_cols = 2 * n_haplotypes
#     resources:
#         memory=lambda wildcards, attempt: get_memory_alloc("get_haplotypes_mem", attempt, 8),
#         runtime=lambda wildcards, attempt: get_runtime_alloc("get_haplotypes_time", attempt, 12),
#     shell:
#         """
#         # Get mom's haplotype

#         if [ {wildcards.haplo} = 1 ]
#         then
#             bcftools view -r {wildcards.chr} -T {input.positions} {input.vcf} | \
#                 cut -f10- | \
#                 sed 's/|[0-9]//g' | \
#                 grep -v "#" | \
#                 gzip - > {output.haplotype}
#         else
#             bcftools view -r {wildcards.chr} -T {input.positions} {input.vcf} | \
#                 cut -f10- | \
#                 sed 's/[0-9]|//g'| \
#                 grep -v "#" | \
#                 gzip - > {output.haplotype}
#         fi
#         """

# rule get_metadata:
#     input:
#         vcf = "VCF/ALL.chr{chr}.phase3{prefix}.vcf.gz",
#         positions = "Haplotypes/ALL.chr{chr}.phase3{prefix}.diallelic.snps.positions"
#     output:
#         metadata = "Haplotypes/ALL.chr{chr}.phase3{prefix}_metadata.txt.gz"
#     shell:
#         """
#         bcftools view -r {wildcards.chr} -T {input.positions} {input.vcf} | \
#             bcftools query -f "%ALT\\t%CHROM\\t%FILTER\\t%ID\\t%POS\\t%QUAL\\t%REF\\n" | \
#             gzip - > \
#             {output.metadata}
#         """