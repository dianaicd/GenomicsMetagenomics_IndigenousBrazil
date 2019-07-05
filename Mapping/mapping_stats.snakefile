configfile: "MN1943.yaml"
mito=config["Mitochondrial"]
import re
# Verify that all necessary files are present:
# BAM:
# 1 final per sample
# 1 rmdup per library
# 1 including dup per library
# Stats from AdapterRemoval:
# 1 per ID
def get_values_column(input, column, parent = False, parentColumn = False):
    # Function to return unique values per column
    myValues = {}
    with open(input, 'r') as file:
        header = file.readline()
        indexColumn = [index for index,element in enumerate(header.split()) 
                if(column == element)][0]
        for line in file.readlines():
            if(parent):
                parentIndex = [index for index,element in enumerate(header.split()) 
                                if(parentColumn == element)][0]
                if(parent == line.split()[parentIndex]):
                    myValues[line.split()[indexColumn]] = 1
            else:
                myValues[line.split()[indexColumn]] = 1
    return(list(myValues.keys()))
def expand_input(input, depth = "ID"):
    # returns prefix for an input
    myPrefix = {}
    #print(type(input))
    #print(input)
    #print('.txt'.join(input))
    input = input+".txt" #; print(input)
    with open(input, 'r') as file:
        header = file.readline()
        for line in file.readlines():
            sm = line.split()[5]
            lb = line.split()[3]
            id = line.split()[0]
    
            if depth == "SM":
                myPrefix[sm+"/"+sm] = 1     
            elif depth == "LB":
                myPrefix[sm+"/"+lb+"/"+lb] = 1  
            elif depth == "ID":
                myPrefix[sm+"/"+lb+"/"+id+"/"+id] = 1 
    return(list(myPrefix.keys()))
def identMe(text, nIdent, colon = True):
    ident = "  " #ident with two spaces
    myString = ""
    for i in range(0, nIdent):
        myString = myString + ident
    myString = myString + text 
    if colon: 
        myString = myString + ":\n"
    else: 
        myString = myString + "\n"
    return(myString)
def input2yaml(input):

    yaml = input.replace(".txt", ".yaml")
    SMs = get_values_column(input, "SM")
    with open(yaml, 'w') as YAML:
        YAML.write(identMe("Samples", 0))
        for sm in SMs:
            YAML.write(identMe(sm, 1))
            LBs = get_values_column(input = input, column = "LB",
                                   parent = sm, parentColumn = "SM")
            for lb in LBs:
                YAML.write(identMe(lb, 2))
                IDs = get_values_column(input = input, column = "ID",
                                   parent = lb, parentColumn = "LB")
                for id in IDs:
                    YAML.write(identMe(id, 3, False))
        YAML.close()

IDs = [expand_input("{file}".format(file=myInput), depth = "ID")
                    for myInput in config["Samples"].keys()][0]
LBs = [expand_input("{file}".format(file=myInput), depth = "LB") 
                    for myInput in config["Samples"].keys()][0]
SMs = [expand_input("{file}".format(file=myInput), depth = "SM") 
                    for myInput in config["Samples"].keys()][0]


rule all:
    input:
        # bamSample = expand("{prefix}.bam", prefix = SMs),
        # baiSample = expand("{prefix}.bam.bai", prefix = SMs),
        # bamLibrary = expand("{prefix}.bam", prefix = LBs),
        # baiLibrary = expand("{prefix}.bam.bai", prefix = LBs),
        # bamIds = expand("{prefix}.bam", prefix = IDs),
        # baiIds = expand("{prefix}.bam.bai", prefix = IDs),
        # smSettings = expand("{prefix}.settings", prefix = SMs),
        # lbSettings = expand("{prefix}.settings", prefix = LBs),
        # idSettings = expand("{prefix}.settings", prefix = IDs)
        # settings = expand("{prefix}.settings", prefix = SMs + LBs + IDs),
        # flagstat = expand("{prefix}_flagstat.txt", prefix = SMs + LBs + IDs),
        # idxstats = expand("{prefix}_idxstats.txt", prefix = SMs + LBs + IDs),
        # bedtools = expand("{prefix}.genomecov", prefix = SMs + LBs + IDs),
        # length = expand("{prefix}.length", prefix = SMs + LBs + IDs),
        # rmdup = expand("{prefix}_rmdup.stats", prefix = SMs)
        summary = expand("{prefix}.summary", prefix = SMs)

rule index_bam:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    shell:
        "samtools index {input}"

def theirChildren(wildcards):
    myStr = wildcards.sample+"/"+wildcards.file
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
    elif myStr in LBs:
        children = expand_input(wildcards.sample, depth = "ID")
    else:
        children = "notWanted"
    #print(myStr)
    #print(children)
    #children = [re.sub(".*/", "", child) for child in children]
    return(children)

def theirChildrenSettings(wildcards):
    myStr = wildcards.sample+"/"+wildcards.file
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
        children = expand("{c}.settings", c = children)
        return(children)
    elif myStr in LBs:
        children = expand_input(wildcards.sample, depth = "ID")
        children = expand("{c}.settings", c = children)
        return(children)
    else:
        return("")

    #return(children)

rule adapterRemoval_settings:
    input:
        bam="{sample}/{file}.bam",
        children=theirChildrenSettings
    output:
        protected("{sample}/{file}.settings")
    wildcard_constraints:
        sample = "\w+"
    params:
        descendants = theirChildren
    shell:
        '~/Projects/Botocudos/Scripts/Mapping/summary_settings.sh {wildcards.sample}/{wildcards.file} {params.descendants}'

rule flagstat:
    input:
        "{file}.bam"
    output:
        "{file}_flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"

rule idxstats:
    input:
        bam = "{file}.bam",
        bai = "{file}.bam.bai"
    output:
        "{file}_idxstats.txt"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule get_chrs:
    input:
        "{file}.bam"
    run:
        chromosomes = shell('samtools view -H {input} |grep "^@SQ" |sed "s/SN://" |cut -f 2')
        print(chromosomes)
        print(type(chromosomes))

rule bedtools:
    input:
        "{file}.bam"
    output:
        "{file}.genomecov"
    shell:
        "bedtools genomecov -ibam {input} > {output}"

rule length:
    input:
        nuc="{file}.bam"
    output:
        nuc="{file}.length",
        mito="{file}.mito.length"
    shell:
        "cat {input} |"
        "python ~/Projects/Botocudos/Scripts/DataQuality/read_length.py -o {output.nuc} ;"
        "samtools view -b {input} {mito} |"
        "python ~/Projects/Botocudos/Scripts/DataQuality/read_length.py -o {output.mito} "

def theirChildren(wildcards):
    myStr = wildcards.sample+"/"+wildcards.sample
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
    else:
        children = "notWanted"
    return(children)

def theirChildrenRmDup(wildcards):
    myStr = wildcards.sample+"/"+wildcards.sample
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
        children = expand("{c}_rmdup.stats", c = children)
        return(children)
    else:
        return("")

rule rmdup:
    input:
        bam="{sample}/{sample}.bam",
        children=theirChildrenRmDup
    output:
        sample=protected("{sample}/{sample}_rmdup.stats")
    wildcard_constraints:
        sample = "\w+"
    params:
        descendants = theirChildren
    shell:
        'for file in {params.descendants} ;'
        'do '
        'rmdup=$(grep -v "#" MN1943/L1/L1_rmdup.stats '
        '|grep -v "^$"|cut -f6 |tail -n1) ;'
        'echo "$file\t$rmdup" >> {output.sample} ; done'


def theirChildrenSettings(wildcards):
    myStr = wildcards.sample+"/"+wildcards.sample
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
        children = expand("{c}.settings", c = children)
        return(children)

def theirChildrenGenomecov(wildcards):
    myStr = wildcards.sample+"/"+wildcards.sample
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
        children = expand("{c}.genomecov", c = children)
        return(children)

def theirChildrenIdxstats(wildcards):
    myStr = wildcards.sample+"/"+wildcards.sample
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
        children = expand("{c}_idxstats.txt", c = children)
        return(children)

def theirChildrenNucLen(wildcards):
    myStr = wildcards.sample+"/"+wildcards.sample
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
        children = expand("{c}.length", c = children)
        return(children)

def theirChildrenMitoLen(wildcards):
    myStr = wildcards.sample+"/"+wildcards.sample
    if myStr in SMs:
        children = expand_input(wildcards.sample, depth = "LB")
        children = expand("{c}.mito.length", c = children)
        return(children)    

rule summary:
    input:
        bam="{sample}/{sample}.bam",
        adapter=theirChildrenSettings,
        rmdup = theirChildrenRmDup,
        genomecov=theirChildrenGenomecov,
        idxstats=theirChildrenIdxstats,
        nucLength=theirChildrenNucLen,
        mitoLength=theirChildrenMitoLen
    output:
        "{sample}/{sample}.summary"
    wildcard_constraints:
        sample = "\w+"
    shell:
        "python ~/Projects/Botocudos/Scripts/DataQuality/summary_alignment.py "
        "--sample {wildcards.sample} --output {output} --mitochondrial {mito}"
