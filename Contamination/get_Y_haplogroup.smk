configfile: "multiple_purposes.yaml"


rule angsd_haploid_call:
    input:
        bamlist = "{bamlist}.txt"
    output:
    params:
        chr_y = 
    shell:
    """
    angsd -bam {input.bamlist} -dohaplocall 1 -doCounts 1 -r {params.chr_y}: -minMinor 0 -out {wildcards.bamlist}
    """

rule format_genos:
    input:
    output:
    shell:

rule query_haplogroup:
    input:
    output:
    shell:
