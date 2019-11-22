# Snakefile to subset variants for enhanced D-statistics
configfile = "multiple_purposes.yaml"

def param_default(name, myDict, default = False):
    if name in myDict.keys():
        parameter = myDict[name]
    else:
        parameter = default
    return(parameter)

max_missing = param_default(name = "max_missing", myDict = config["enhanced"], default = 0.5)

#=============================================================================#
rule all:
    input:
        enhanced_panels = expand("panel_enhanced_h3_{pop_h3}_h4_{pop_h4}_anc_{ancestral}.bed")
        
rule get_anc_variants:
    input:
        fasta = "{ancestral}.fa",
        bim = "{panel}.bim"
    output:
        variants = "{ancestral}.variants"
    shell:
        """
        python3.6 subset_chimp.py --bim {input.bim} \
        --ancestral {input.fasta} --out {output.variants}
        """


rule subset_H3_tped:
    input:
        bed = "{panel}.bed",
        h3 = "{pop_h3}.txt"
    output:
        tped = "{pop_h3}.txt"
    shell:
        """
        plink --recode trasnpose --bfile {wildcards.panel} \
            --out {wildcards.pop_h3} \
            --keep {input.h3}
        """

rule subset_H3_free_tped:
    input:
        bed = "{panel}.bed",
        h4 = "{pop_h4}.txt"
    output:
        tped = "{pop_h4}.tped"
    shell:
        """
        plink --recode transpose --keep {input.h4} \
        --bfile {wildcards.panel} \
        --geno {max_missing} \
        --out {wildcards.pop_h4}
        """

rule get_enhanced_vars:
    input:
        anc_variants = "{ancestral}.variants",
        h3_tped = "{pop_h3}.tped",
        h4_tped = "{pop_h4}.tped"
    output:
        enhanced_variants = "enhanced_h3_{pop_h3}_h4_{pop_h4}_anc_{ancestral}.variants"
    shell:
        """
        python3.6 subset_enhance.py \
        --ancestral {input.anc_variants} \
         --H3 {input.h3_tped} \
         --H3_free {input.h4_tped} \
         --out {output.enhanced_variants}
        """

rule subset_panel_enhanced:
    input:
        bed = "{panel}.bed",
        enhanced_variants = "enhanced_h3_{pop_h3}_h4_{pop_h4}_anc_{ancestral}.variants"
    output:
        panel = "panel_enhanced_h3_{pop_h3}_h4_{pop_h4}_anc_{ancestral}.bed"
    shell:
        """
        plink --extract {input.enhanced_variants} \
            --bfile {wildcards.panel} \
            --out panel_enhanced_h3_{wildcards.pop_h3}_h4_{wildcards.pop_h4}_anc_{wildcards.ancestral} \
            --make-bed
        """