# pipeline to merge a bcf split into chromosomes to a single PLINK file
configfile: "merge_bcf_plink.yaml"

panel = list(config["panel"].keys())[0]
ind = list(config["ind"].keys())[0]

rule all:
    input:
        merged = f"{panel}/{panel}_{ind}.bed"

rule get_positions_chr:
    input:
        panel = lambda wildcards: config["panel"][wildcards.panel]["path"]
    output:
        panel_chr = temp("{panel}.chr{chr}.bed"),
        positions_chr = temp("{panel}.positions.chr{chr}.bed")
    params:
        basename_input = lambda wildcards: config["panel"][wildcards.panel]["path"].replace(".bed", ""),
        basename_output = "{panel}.chr{chr}"
    shell:
        """
        plink --bfile {params.basename_input} \
            --chr {wildcards.chr} --out {params.basename_output} \
            --snps-only 'just-acgt' --make-bed

        awk 'BEGIN {{OFS="\t"}} {{print $1,$4-1,$4}}' \
        {params.basename_output}.bim  > {output.positions_chr}

        """


rule index_bcf:
    input:
        bcf = "{file}.bcf"
    output:
        tbi = "{file}.bcf.csi"
    shell:
        """
        tabix {input.bcf}
        """

rule extract_positions_vcf:
    input:
        bcf = lambda wildcards: config["ind"][wildcards.ind]["prefix"] + f"{wildcards.chr}.bcf",
        tbi = lambda wildcards: config["ind"][wildcards.ind]["prefix"] + f"{wildcards.chr}.bcf.csi",
        positions_chr = "{panel}.positions.chr{chr}.bed"
    output:
        vcf = temp("{panel}/{ind}.chr{chr}.vcf")
    shell:
        """
        bcftools view -r {wildcards.chr} -R {input.positions_chr} \
            -Ov -e 'INDEL=1' -o {output.vcf} \
            {input.bcf}
        """

rule set_id_vcf:
    input:
        vcf = "{panel}/{file}.chr{chr}.vcf"
    output:
        vcf = temp("{panel}/{file}.chr{chr}.new_id.vcf")
    shell:
        """
        bcftools annotate -Ov -o {output.vcf}  --set-id +'%CHROM\_%POS' {input.vcf}
        """

rule concat_vcf:
    input:
        vcfs = expand(
            "{panel}/{file}.chr{chr}.new_id.vcf", 
            panel = "{panel}",
            file = "{file}",
            chr = [str(chr) for chr in range(1, 23)]
            )
    output:
        vcf = "{panel}/{file}.new_id.vcf"
    shell:
        """
        bcftools concat -Ov -o {output.vcf} {input.vcfs}
        """

rule bcf_to_bed:
    input:
        vcf = "{panel}/{file}.new_id.vcf"
    output:
        bed = "{panel}/{file}.bed"
    params:
        basename = "{panel}/{file}"
    shell:
        """
        plink --make-bed --vcf {input.vcf} --out {params.basename}
        """

rule merge_beds:
    input:
        panel_bed = lambda wildcards: config["panel"][wildcards.panel]["path"],
        ind_bed = "{panel}/{ind}.bed"
    output:
        merged_bed = "{panel}/{panel}_{ind}.bed"
    params:
        panel_basename = lambda wildcards: config["panel"][wildcards.panel]["path"].replace(".bed", ""),
        ind_basename = "{panel}/{ind}",
        out_basename = "{panel}/{panel}_{ind}"
    shell:
        """
        plink --bfile {params.panel_basename} \
            --bmerge {params.ind_basename} --make-bed \
            --out {params.out_basename} --merge-equal-pos
        """