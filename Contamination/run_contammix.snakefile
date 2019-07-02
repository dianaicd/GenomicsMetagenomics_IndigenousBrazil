configfile: "2019_07_01.yaml"

mito=config["Mitochondrial"]
samples=config["samples"].split()

rule all:
    input:
        # fasta=expand('{mito}/{sample}.fa', mito = mito, sample = samples),
        # dict=expand('{mito}/{sample}.dict', mito = mito, sample = samples),
        # bam=expand('{mito}/{sample}.bam', mito = mito, sample = samples),
        # realign=expand("{mito}/{sample}_{mito}_consensus.bam", mito = mito, sample = samples),
        # maln=expand("{mito}/311_{sample}_aligned.fasta", mito = mito, sample = samples),
        fig=expand("{mito}/{sample}_{rmtrans}.pdf", mito = mito, sample = samples, rmtrans = ['all', 'rmTrans'])

def inputBam(wildcards):
    myInput=wildcards.sample+"/"+wildcards.sample+".bam"
    return(myInput)
def inputBai(wildcards):
    myInput=wildcards.sample+"/"+wildcards.sample+".bam.bai"
    return(myInput)

rule index:
    input:
        "{x}.bam"
    output:
        "{x}.bam.bai"
    shell:
        "samtools index {input}"

rule consensus:
    input:
        bam=inputBam,
        bai=inputBai
    params:
        mito=mito
    output:
        fasta='{mito}/{sample}.fa',
        dict='{mito}/{sample}.dict',
        bam="{mito}/{sample}.bam"
    shell:
        '  samtools view -b {input.bam} {mito} > {output.bam} ;'
        '  samtools index {output.bam} ;'
        '  angsd -doFasta 2 -doCounts 1 -i {output.bam} -out  {wildcards.sample}/{mito}/{wildcards.sample} ;'
        '  gunzip {output.fasta}.gz ;'
        '  bwa index  {output.fasta} ;'
        '  samtools faidx {output.fasta} ;'
        '  picard-tools CreateSequenceDictionary REFERENCE={output.fasta} OUTPUT={output.dict} ;'

def inputBam(wildcards):
    return(mito+"/"+wildcards.sample+".bam")


rule realign:
    input:
        bam="{mito}/{sample}.bam",
        bai="{mito}/{sample}.bam.bai"
    output:
        bam="{mito}/{sample}_{mito}_consensus.bam",
        fastq="{mito}/{sample}_{mito}_consensus.truncated",
        fasta="{mito}/{sample}_{mito}.fa"
    shell:
        "bedtools bamtofastq -i {input.bam} -fq {output.fastq} ;"
        "gzip {output.fastq} ;"
        "aDNA_mapping.sh --ref {output.fasta} --fastq1 {output.fastq} --skip1 1 ;"
        "samtools index {input.bam} ;"

rule maln:
    input:
        "{mito}/{sample}_{mito}.fa"
    output:
        "{mito}/311_{sample}_aligned.fasta"
    params:
        mt311=config["mt311"]
    shell:
        "cat {params.mt311} {input} > 311_{wildcards.sample}.fasta ;"
        "mafft --auto 311_{wildcards.sample}.fasta > {output} ;"


rmTrans={"all":"", "rmTrans":"--transverOnly"}

rule contammix:
    input:
        bam="{mito}/{sample}_{mito}_consensus.bam",
        bai="{mito}/{sample}_{mito}_consensus.bam.bai",
        maln="{mito}/311_{sample}_aligned.fasta"
    output:
        fig="{mito}/{sample}_{rmTrans}.pdf",
        data="{mito}/{sample}_{rmTrans}.Rda"
    params:
        transitions=lambda wildcards: rmTrans.get(wildcards.rmTrans),
        nIter=config["nIter"],
        path=config["contammix"],
        prefix="{mito}/{sample}_{rmTrans}"
    shell:
        "{params.path} --samFn {input.bam} --nIter {params.nIter} "
        "--malnFn {input.maln} --alpha 0.1 --figure {output.fig} "
        "--saveData {output.data} {params.transitions}"
