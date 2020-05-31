import glob
configfile: "pe_reads.yaml"
include: "parse_resources.smk"
ruleorder: adapter_removal > bwa_aln > select_mapQ > sort_bam > merge_LB > merge_SM

def get_fastqFiles(sample, library, fileID):
    fastqFiles = glob.glob(config['pe_reads'][sample][library][fileID])
    if "_R1_" in fastqFiles[1]:
        fastqFiles = [fastqFiles[1], fastqFiles[0]]
    return(fastqFiles)

ref = config['ref']
baseQ = config['baseQ'] if 'baseQ' in config.keys() else 20
mapQ = config['mapQ'] if 'mapQ' in config.keys() else 30
minLength = config['minLength'] if 'minLength' in config.keys() else 30

samples = config['pe_reads'].keys()

libraries = {sample:list(config['pe_reads'][sample].keys()) 
                for sample in samples}

fastqIDs = {sample+"_"+library:list(config['pe_reads'][sample][library].keys()) 
            for sample in samples
                for library in libraries[sample]
                }

rule all:
    input:
        expand("{sample}/{sample}.mapQ{mapQ}.realn.calmd.bam", sample = samples, mapQ = mapQ)

rule samtools_index:
    input:
        bam = "{file}.bam"
    output:
        bam = "{file}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

rule adapter_removal:
    input:
        fastq = lambda wildcards: get_fastqFiles(wildcards.sample, wildcards.library, wildcards.id)
    output:
        fastq = "{sample}/{library}/{id}/{id}.collapsed.gz"
    params:
        basename = "{sample}/{library}/{id}/{id}",
        minLength = minLength
    shell:
        """
        AdapterRemoval --file1 {input.fastq[0]} \
            --file2 {input.fastq[1]} \
            --basename {params.basename} \
            --trimns --trimqualities \
            --qualitybase 33 \
            --minlength {params.minLength} \
            --collapse --gzip
        """

rule bwa_aln:
    input:
        ref = ref,
        fastq = "{file}.collapsed.gz"
    output:
        sai = "{file}.sai"
    params:
    threads: 8
    shell:
        """
        bwa aln {input.ref} -l 500 -t {threads} {input.fastq} > {output.sai}
        """
rule bwa_samse:
    input:
        ref = ref,
        fastq = "{file}.collapsed.gz",
        sai = "{file}.sai"
    output:
        bam = "{file}.all.bam"
    params:
        sm = "{file}".split('/')[0],
        id = lambda wildcards: wildcards.file.split("/")[-1]
    threads: 8
    shell:
        """
        bwa samse {input.ref} \
            -r '@RG\\tID:{params.id}\\tSM:{params.sm}' \
            {input.sai} {input.fastq} | \
            samtools view -bSh - > {output.bam}
        """

rule select_mapQ:
    input:
        bam = "{file}.all.bam"
    output:
        bam = "{file}.mapQ{mapQ}.bam".format(mapQ = mapQ, file = "{file}")
    shell:
        """
        samtools view -bh -F4 -q {mapQ} {input.bam} > {output.bam}
        """

rule sort_bam:
    input:
        bam = "{file}.mapQ{mapQ}.bam".format(mapQ = mapQ, file = "{file}")
    output:
        bam = "{file}.mapQ{mapQ}.sort.bam".format(mapQ = mapQ, file = "{file}")
    shell:
        """
        samtools sort -o {output.bam} {input.bam}
        """

rule merge_LB:
    input:
        bam = lambda wildcards: expand("{sample}/{library}/{id}/{id}.mapQ{mapQ}.sort.bam", 
                                        sample = wildcards.sample, 
                                        library = wildcards.library, 
                                        id = list(config['pe_reads'][wildcards.sample][wildcards.library].keys()),
                                        mapQ = mapQ)
    output:
        bam = "{sample}/{library}/{library}.mapQ{mapQ}.merge.bam"
    shell:
        """
        samtools merge {output.bam} {input.bam}
        """

rule samtools_rmdup:
    input:
        bam = "{file}.merge.bam"
    output:
        bam = "{file}.rmdup.bam"
    shell:
        """
        samtools rmdup -s {input.bam} {output.bam}
        """

rule select_unique:
    input:
        bam = "{file}.sort.rmdup.bam"
    output:
        bam = "{file}.sort.rmdup.unique.bam"
    shell:
        """
        samtools view -h {input.bam} | \
            grep -v 'XT:A:R|XA:Z' | \
            awk '{{if($0~/X1:i:0/||$0~/^@/)print $0}}' | \
            samtools view -bS - > {output.bam}
        """

rule merge_SM:
    input:
        bam = lambda wildcards: expand("{sample}/{library}/{library}.mapQ{mapQ}.sort.rmdup.unique.bam", 
                                        sample = wildcards.sample,
                                        library =  libraries[wildcards.sample],
                                        mapQ = mapQ)
    output:
        bam = "{sample}/{sample}.mapQ{mapQ}.merge.bam".format(mapQ = mapQ, sample = "{sample}")
    shell:
        """
        samtools merge {output.bam} {input.bam}
        """

rule realign_gatk:
    """
    Realign sequence around indels.
    """
    input:
        ref = ref,
        # fai = ref + ".fai",
        # dict = ref + ".dict",
        bam = "{file}.merge.bam",
        bai = "{file}.merge.bam.bai"
    output:
        bam = "{file}.realn.bam",
        intervals = "{file}.realn.intervals"
    threads: 8
    shell:
        """
        GenomeAnalysisTK -I {input.bam} -R {input.ref} \
            -T RealignerTargetCreator \
            -o {output.intervals} \
            --num_threads {threads}

        GenomeAnalysisTK -I {input.bam} -T IndelRealigner \
            -R {input.ref} -targetIntervals {output.intervals} \
            -o {output.bam} 
        """
            
rule samtools_calmd:
    input:
        ref = ref,
        bam = "{file}.realn.bam"
    output:
        bam = "{file}.realn.calmd.bam"
    threads: 8
    shell:
        """
        samtools calmd --threads {threads} {input.bam} {input.ref} |\
            samtools view -bS - > {output.bam}
        """
