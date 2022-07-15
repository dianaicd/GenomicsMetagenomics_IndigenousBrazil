prefix="unmapped_fnas"

samples = [
"MN00010",
"MN00013",
"MN00016",
"MN00019",
"MN00021",
"MN00022",
"MN00023",
"MN00039",
"MN0003",
"MN00045",
"MN00056",
"MN00064",
"MN00066",
"MN00067",
"MN00068",
"MN00069",
"MN0008",
"MN0009",
"MN00118",
"MN00119",
"MN00316",
"MN00346",
"MN01701",
"MN1943"
]

rule all:
    input:
        classification = expand("DIAMOND_output/{sample}/{sample}_tax.txt", sample = samples),
        fast = expand("DIAMOND_output/{sample}/{sample}_tax_fast.txt", sample = samples),
        reads_taxon = expand("DIAMOND_output/{sample}/{sample}_ReadspTaxon{mode}.txt", sample = samples, mode = ["", "_fast"])

rule run_DIAMOND:
    input:
        lambda wildcards: f"{prefix}/{wildcards.sample}_unmapped.fna"
    output:
        "DIAMOND_output/{sample}/{sample}_tax.txt"
    params:
        mode = 'blastx',
        database = "/users/dcruzdav/popgen/Panels/phibase/2021_12_07/phi.dmnd"
    resources:
        memory = 1024 * 8,
        runtime = 6*60 
    threads: 16
    log:
        "logs/DIAMOND/{sample}/{sample}_run_DIAMOND.log"
    shell:
        '''
        module load gcc/9.3.0 diamond/2.0.9
        mkdir -p DIAMOND_output/{wildcards.sample}
        diamond {params.mode} -p {threads} -d {params.database} -q {input} \
        -o {output} -f 102
        '''

rule run_DIAMOND_fast:
    input:
        lambda wildcards: f"{prefix}/{wildcards.sample}_unmapped.fna"
    output:
        "DIAMOND_output/{sample}/{sample}_tax_fast.txt"
    params:
        mode = 'blastx',
        database = "/users/dcruzdav/popgen/Panels/phibase/2021_12_07/phi.dmnd"
    resources:
        memory = 1024 * 8,
        runtime = 6*60
    threads: 16
    log:
        "logs/DIAMOND/{sample}/{sample}_run_DIAMOND.log"
    shell:
        '''
        module load gcc/9.3.0 diamond/2.0.9
        mkdir -p DIAMOND_output/{wildcards.sample}
        diamond {params.mode} -p {threads} -d {params.database} -q {input} \
        -o {output} --id 90
        '''


rule reads_per_taxon_DIAMOND:
    input:
        "DIAMOND_output/{sample}/{sample}_tax{type}.txt"
    output:
        "DIAMOND_output/{sample}/{sample}_ReadspTaxon{type}.txt"
    params:
        runtime = "02:00:00",
        empty_ReadsTaxon_NamesRanks = "DIAMOND_output/{sample}/{sample}_ReadsTaxon_NamesRanks.tsv",
        empty_Correct_Incorrect = "DIAMOND_output/{sample}/{sample}_Correct_Incorrect.tsv"
    resources:
        memory = 3000
    shell:
        '''
        set +e
        
        cut -f2 {input} | grep -wv '0' | sort | uniq -c | awk '{{print $1"\t"$2}}' > {output}
        
        exitcode=${{PIPESTATUS[1]}}
        if [[ exitcode -eq 1 ]]; then
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0}}' > {params.empty_ReadsTaxon_NamesRanks};
            awk 'BEGIN{{print 0"\\t"0"\\t"0"\\t"0"\\t"0"\\tUnclassified\\t"0}}' > \
            {params.empty_Correct_Incorrect}
            exit 0
        fi
        '''

