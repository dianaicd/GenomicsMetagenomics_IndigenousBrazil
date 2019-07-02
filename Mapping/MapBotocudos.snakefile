configfile: "2019_06_13.yaml"
        #lambda wildcards: config["samples"][wildcards.sample]
rule all:
    input:
        bam=expand("{value}/{value}.bam", value = config["samples"].values()),

#print(config["samples"])
def get_info(x):
    for key, value in config["samples"].items():
        if value == x:
            return(key)

rule listSamples:
    output:
        "{value}.txt"
    params:
        key=lambda wildcards: get_info(wildcards.value)
    shell:
        'Rscript ~/data/Git/Botocudos-scripts/Mapping/makeList.r {params.key} {wildcards.value}'
        # 'echo "ID\tData\tMAPQ\tLB\tPL\tSM" > {output};'
        # 'libs=($(ls FASTQ/{params.key}*/*  |rev |cut -f1 -d M | rev |cut -f1 -d. |cut -f2- -d_ |cut -f1 -d_))'
        # 'for l in ${libs[@]} ; do'
        # '   '
        # 'libs=($(ls FASTQ/{params.key}/ | sed "s/.*_{params.key}_L/{params.key}_L/ "|cut -f2 -d_ |sort |uniq)) ;'
        # 'for l in ${{libs[@]}} ; do'
        # '   fastq=($(ls -R | '
        # '   awk -f fullpath.awk |grep {params.key} | grep -P ".*${{l}}_.*" |sed "s/\./\.\./"|sort |uniq)) ; '
        # '   for f in ${{fastq[@]}} ; do'
        # '       id=$(echo $f |sed "s/.fastq.*// ; s|.*FASTQ/||; s|/|_|g") ;'
        # '       echo "$id\t$f\t30\t$l\tILLUMINA\t{wildcards.value}" >> {output} ; '
        # '   done ; '
        # 'done; '
            
rule map:
    input:
        "{value}.txt"
    threads: 12
    output: 
        "{value}/{value}.bam"
    params:
        ref=config["ref"]
    log:
        "{value}/logs/mapping.txt"
    shell:
        "mapping_aDNA.sh -i {input} --ref {params.ref} -p {threads} --cpuPerJob {threads}"
