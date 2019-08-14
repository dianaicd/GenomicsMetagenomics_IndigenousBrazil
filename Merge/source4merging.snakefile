# Snakefile version to merge datasets
configfile: "samples_panel.yaml"

bamlists = {key:value for key,value in config["samples"]} #config["samples"].values()
chromosomes = config["chromosomes"]
panel = config["panel"]["name"]

rule all:
    input:
        bamlist = expand("{bamlist}.txt", bamlist = bamlists.keys()),
        mpileup = expand("{bamlist}_{chr}.mpileup", bamlist = bamlists.keys(), chr = chromosomes)


rule make_bamlist:
    output:
        "{group}.txt"
    shell:
        "ls {bamlists[wildcards.group]}/*.bam > {output}"

rule mpileup:
    input:
        bamlist = "{bamlist}",
        panel = "{panel}"
    output:
        "{bamlist}_{chr}.mpileup"
    params:
        baseQ = config["BaseQuality"]    
    log:
        "out_mpileup_{chr}.txt"
    shell:
        "samtools mpileup -r {wildcards.chr} -Bl {wildcards.panel}_{wildcards.chr}_sites.bed "
        "   -b {wildcards.bamlist} -a -o {wildcards.bamlist}_{wildcards.chr}.mpileup -Q {params.baseQ} "
        "   2> {log}  "


panelExtension={"vcf":"vcf", "bed":"bim", "bim":"bim"}
columns = {"vcf":"1,2,4,5", "bim":"1,4,5,6"}

rule refalt:
    input:
        panel = "{panel}.{extension}"
    params:
        extension = panelExtension[config["PanelType"]],
        col = columns[extension]
    output:
        refalt = "{panel}.refalt"
    shell:
        "cut -f {params.col} {input.panel} "
        "   | sed 's/\t/_/ ;s/A/0/g ; s/C/1/g; s/G/2/g ; s/T/3/g' > {output.refalt}"

rule refalt_chr:
    input:
        refalt = "{panel}.refalt"
    params:
        chromosomes = config["chromsomes"]
    output:
        refalt = "{panel}_{chr}.refalt",
        sites = "{panel}_{chr}_sites.bed"
    shell:
        'cut -f 1,2 {wildcards.panel}.refalt |grep -P "^{wildcards.chr}\_"| sed "s/_/\t/"| '
        'awk "{print($1"\t"$2-1"\t"$2)} > {wildcards.panel}_{wildcards.chr}_sites.bed ; '
        'grep -P "^{wildcards.chr}\_" {wildcards.panel}.refalt > {wildcards.panel}_{wildcards.chr}.refalt '

rule count_and_sample:
    input:
        mpileup = "{bamlist}_{chr}.mpileup",
        refalt = "{panel}_{chr}.refalt"
    output:
        counts = "{bamlist}_{chr}.counts.gz",
        sampled = "{bamlist}_{chr}.sampled.gz"
    params:
        ped = config[],
        allmutations = config["allmutations"]
    log:
    shell:
        "python count_and_sample.py --mpileup {input.mpileup} "
        " --counts {output.counts} "
        " --sampled {output.sampled} "
        " --refalt {input.refalt} {params.ped} {params.allmutations}"


rule merge_sampled:
    input:
        sampled = expand("{bamlist}_{chr}.sampled.gz", 
                         bamlist = wildcards.bamlist, chr = config["chromosomes"])
    output:
        merged = "{bamlist}.counts.sampled.txt"
    log:
    shell:
        "zcat {input.sampled} > {output.merged}"

#=============================================================================#
# If you want bases to be trimmed at the end of the reads,
# you might want to use ANGSD to count bases
#=============================================================================#

rule sample_from_vcf:
    input:
        panel = "{panel}.vcf"
    output:
        counts = "{panel}.counts",
        sampled = "{panel}.counts.sampled.txt"
    log:
    shell:
        "python vcf_sample_random.py {input.panel} {output.counts} {output.sampled}"

rule genotype_likelihoods:
    input:
        panel = "{panel}.{type}",
        bamlist = "{bamlist}"
    output:
        panel = "{panel}.beagle",

    params:
    log:
    shell:
        "./workflow_genolike.sh --panel {input.panel} --homozygous $homo "
        "   --bamlist {input.bamlist} --rmdamage $rmdamage"

rule bed_to_tped:
    input:
        panel = "{panel}.bed"
    output:
        panel = "{panel}.tped",
        tfam = "{panel}.tfam"
    log:
    shell:
        "plink --recode transpose --bfile {wildcards.panel} --out {wildcards.panel} "

rule bed_to_vcf:
    input:
        panel = "{panel}.bed"
    output:
        panel = "{panel}.vcf"
    log:
    shell:
        "plink --recode vcf --bfile {wildcards.panel} --out {wildcards.panel}"

rule merge_bam_to_bed:
    input:
        panel_tped = "{panel}.tped",
        panel_tfam = "{panel}.tfam",
        counts = "{bamlist}.counts.sampled.txt",
        bamlist = "{bamlist}.txt"
    output:
        merged_tped = "{panel}.{bamlist}.tped",
        merged_tfam = "{panel}.{bamlist}.tfam",
        temp_fam = temp("{panel}.{bamlist}.temp.fam"),
        merged_fam = "{panel}.{bamlist}.fam"
    log:
    shell:
        "paste {input.panel_tped} {input.counts} -d ' ' > {output.merged_tped} "
        "cp {input.panel_tfam} {output.merged_tfam} "    

        "while read line ;"
        "do "
        "    name=$(basename $line .bam) ;"
        "    echo '$name $name 0 0 0 -9' >> {output.merged_tfam} ;"
        "done < {input.bamlist} ;"

        "plink --recode --make-bed --tfile {wildcards.panel}.{wildcards.bamlist} --out {wildcards.panel}.{wildcards.bamlist} ;"
        "awk -F'\t' '{print $1,$2,$3,$4,$5,1}' "
        "   {output.merged_fam} > {output.temp_fam} ; mv {output.temp_fam} {output.merged_fam} "

rule tped_to_bed:

rule bed_to_ped:
    input:
        bed = "{panel}.{bamlist}.bed"
    output:
        ped = "{panel}.{bamlist}.ped"
    log:
    shell:
        "plink --recode --bfile {wildcards.panel}.{wildcards.bamlist} --out {wildcards.panel}.{wildcards.bamlist}"

rule ped_to_haploid:
    input:
        ped = "{panel}.{bamlist}.ped",
        map = "{panel}.{bamlist}.map"
    output:
        haploid = "{panel}.{bamlist}.haploid.ped",
        map = "{panel}.{bamlist}.haploid.map"
    log:
    shell:
        "perl sample_ped.pl -ped {input.ped} -out {output.haploid} ;"
        "cp {input.map} {output.map}"

rule calc_dist:
    input:
        ped = "{panel}.{bamlist}.haploid.ped"
    output:
        dist = "{panel}.{bamlist}.mind{mind}.dist.gz"
    wildcard_constraints:
    log:
    params:

    shell:
        " plink --distance square gz 'flat-missing' 
        "    --file {wildcards.panel}.{wildcards.bamlist}.haploid "
        "    --mind {params.mind} "
        "    --out {wildcards.panel}.{wildcards.bamlist}.mind{params.mind}.haploid"

rule ped_to_fred:
    input:
        ped = "{panel}.{bamlist}.haploid.ped"
    output:
        ped = "{panel}.{bamlist}.haploid.fred.ped",
        tped = "{panel}.{bamlist}.haploid.fred.tped",
        mind_tped = "{panel}.{bamlist}.haploid.fred.mind0.95.tped"
    log:
    shell:
        "plink --recode 12 --file {wildcards.panel}.{wildcards.bamlist}.haploid "
        "   --out {wildcards.panel}.{wildcards.bamlist}.haploid.fred ;"
    
        "plink --recode transpose --file {wildcards.panel}.{wildcards.bamlist}.haploid.fred "
        "   --out {wildcards.panel}.{wildcards.bamlist}.haploid.fred"

        "mind=0.95 ;"
        "plink --tfile {wildcards.panel}.{wildcards.bamlist}.haploid.fred --make-bed "
        "   --mind $mind --out {wildcards.panel}.{wildcards.bamlist}.haploid.fred.mind${{mind}}"

        "plink --recode transpose --bfile {wildcards.panel}.{wildcards.bamlist}.haploid.fred.mind${{mind}}"
        " --out {wildcards.panel}.{wildcards.bamlist}.haploid.fred.mind${{mind}}"

        "lastField=$(head -n1 {output.mind_tped}|wc -w) ;"
        "columns=$(echo $(seq 5 2 $lastField) | sed 's/ /,/g') ;"
        "cut -f $columns -d ' ' {output.mind_tped} "
        " > {wildcards.panel}.{wildcards.bamlist}.haploid.fred.mind${mind}"

# rule par_ped_to_eigenstrat:
#     input:
#         ped = "{panel}.{bamlist}.ped",
#         bim = "{panel}.{bamlist}.bim"
#     output:
#         par = "{panel}.{bamlist}_ped2eigenstrat.par"
#     log:
#     shell:
#         'argument=(genotype snpname indivname outputformat genotypeoutname '
#         '            snpoutname indivoutname familynames) ;'

#         'options=({input.ped} {input.bim} {input.ped} EIGENSTRAT {wildcards.panel}.{wildcards.bamlist}.eigenstratgeno '
#         '        {wildcards.panel}.snp {wildcards.panel}.ind YES) ;'
        
#         'for i in $(seq 0 7) ;'
#         'do '
#         '  echo "${argument[$i]}: ${options[$i]}" >> {output.par} ;'
#         'done'

# rule par_eigenstrat_to_ped:
