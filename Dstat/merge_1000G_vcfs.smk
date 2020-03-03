# Snakemake to merge VCFs from 1000 Genomes
configfile = "multiple_purposes.yaml"

chromosomes = [str(i) for i in range(1,23)]
chromosomes.append("X")
chromosomes.append("Y")

rule merge_vcfs:
    input:
        expand("ALL.chr{chr}.phase3_integrated_v2a.20130502.genotypes.vcf.gz",
        chr = chromosomes)
    output:
        bcf = "1000G_merged.bcf",
        myList = "1000G_list.txt"
    run:
        with open(myList, "w") as bcf_list:
            [bcf_list.write(line + "\n") for line in input]
        
        myCommand = "bcftools merge --force-samples -Ob -o " + output.bcf + " -l " + output.myList

rule find_non_diallelic:
