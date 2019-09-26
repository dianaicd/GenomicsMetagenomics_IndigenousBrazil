configfile: "samples.yaml"

mito=config["Mitochondrial"]
wildcard_constraints:
   # file = "(?P<sample>\w+)\/((?P=sample)$|((?P<library>\w+)\/(?P=library)$))",
    mito = mito

def output_prefix():
    sample_list = list(config["Samples"].keys())
    samples = ["{}/{}".format(sample, sample) for sample in config["Samples"].keys()]
    libs = []
    for sample in sample_list:
        libs.append([#"{sample}/{lib}/{lib}".format(sample = sample, lib = lib) 
        "{sample}/{lib}/library_rmdup/{lib}".format(sample = sample, lib = lib) 
                     for lib in list(config["Samples"][sample]["Libraries"].keys())])

    output = [item for sublist in libs for item in sublist] + samples

    return(output)

all_prefix = output_prefix()

rule all:
    input:
        # fasta=expand('{mito}/{file}_{mito}.fa',
        #                  mito = mito, file = all_prefix),
        # dict = expand('{mito}/{file}.dict', mito = mito, file = all_prefix),
        # bam = expand('{mito}/{file}.bam', mito = mito, file = all_prefix),
        # bai = expand('{mito}/{file}.bam.bai', mito = mito, file = all_prefix),
        # realn_bam = expand('{mito}/{file}_{mito}_consensus.bam', mito = mito, file = all_prefix),
        # realn_bai = expand('{mito}/{file}_{mito}_consensus.bam.bai', mito = mito, file = all_prefix),
        # fastq = expand('{mito}/{file}_{mito}_consensus.truncated.gz', mito = mito, file = all_prefix),
        # maln = expand('{mito}/{file}_311_aligned.fasta', mito = mito, file = all_prefix),
        result_data = expand('{mito}/{file}_{rmtrans}.{extension}',
                         mito = mito, file = all_prefix, 
                         rmtrans = ['all', 'rmTrans'], extension = ['pdf', 'Rdata']),
        result_txt = expand("{mito}/{file}_{rmTrans}.txt", 
                            mito = mito, file = all_prefix, rmTrans = ['all', 'rmTrans'])


def inputBam(wildcards):
    myInput=wildcards.sample+"/"+wildcards.sample+".bam"
    return(myInput)
def inputBai(wildcards):
    myInput=wildcards.sample+"/"+wildcards.sample+".bam.bai"
    return(myInput)

fasta = ["{mito}/{file}_{mito}.fa".format(mito = mito, file = file) for file in all_prefix]
bam_out = ["{mito}/{file}.bam".format(mito = mito, file = file) for file in all_prefix]
dict = ["{mito}/{file}.dict".format(mito = mito, file = file) for file in all_prefix]

def all_bam(wildcards):
    return("{file}.bam".format(file = wildcards.file))
    
def all_bai(wildcards):
    return("{file}.bam.bai".format(file = wildcards.file))

rule consensus:
    input:
        bam = protected("{file}.bam"), #all_bam, 
        bai = protected("{file}.bam.bai") #all_bai
    output:
        fasta = "{mito}/{file}_{mito}.fa", 
        dict = "{mito}/{file}.dict",
        bam = "{mito}/{file}.bam" 
    shell:
        '  mkdir -p {mito}/{wildcards.file} ;'
        '  samtools view -b {input.bam} {mito} > {output.bam} ;'
        '  samtools index {output.bam} ;'
        '  angsd -doFasta 2 -doCounts 1 -i {output.bam} -out {mito}/{wildcards.file}_{mito} ;'
        '  gunzip {output.fasta}.gz ;'
        '  bwa index  {output.fasta} ;'
        '  samtools faidx {output.fasta} ;'
        '  picard-tools CreateSequenceDictionary REFERENCE={output.fasta} OUTPUT={output.dict} ;'

rule index:
    input:
        "{x}.bam"
    output:
        "{x}.bam.bai"
    shell:
        "samtools index {input}"

rule realign:
    input:
        bam = "{mito}/{file}.bam", 
        bai = "{mito}/{file}.bam.bai", 
        fasta = "{mito}/{file}_{mito}.fa" 
    output:
        bam="{mito}/{file}_{mito}_consensus.bam",
        fastq="{mito}/{file}_{mito}_consensus.truncated.gz"
    params:
        basename="{file}_{mito}_consensus"
    log:
        "logs/{mito}/{file}_{mito}_consensus.log"
    shell:
        "bedtools bamtofastq -i {input.bam} -fq {mito}/{params.basename}.truncated ;"
        "gzip {mito}/{params.basename}.truncated ;"
        "aDNA_mapping.sh --ref {input.fasta} --fastq1 {params.basename} "
        "   --skip1 1 --override --base {params.basename} "
        "   --wd MT -o {mito}/{params.basename}.bam 2>{log};"

rule maln:
    input:
        "{mito}/{file}_{mito}.fa"
    output:
        "{mito}/{file}_311_aligned.fasta"
    params:
        mt311=config["mt311"]
    shell:
        "cat {params.mt311} {input} > {mito}/{wildcards.file}_311.fasta ;"
        "mafft --auto {mito}/{wildcards.file}_311.fasta > {output} ;"


rmTrans={"all":"", "rmTrans":"--transverOnly"}


rule contammix:
    input:
        bam="{mito}/{file}_{mito}_consensus.bam",
        bai="{mito}/{file}_{mito}_consensus.bam.bai", 
        maln="{mito}/{file}_311_aligned.fasta" 
    output:
        fig=protected("{mito}/{file}_{rmTrans}.pdf"),
        data=protected("{mito}/{file}_{rmTrans}.Rdata")
    params:
        transitions=lambda wildcards: rmTrans.get(wildcards.rmTrans),
        nIter=config["nIter"],
        path=config["contammix"],
        prefix="{mito}/{file}_{rmTrans}",
        basename="{mito}/{file}_{rmTrans}"
    # log:
    #     "{mito}/{file}_{rmTrans}.log"
    shell:
        "Rscript {params.path} --samFn {input.bam} --nIter {params.nIter} "
        "--malnFn {input.maln} --alpha 0.1 --figure {output.fig} "
        "--saveData {params.basename} {params.transitions}"

rule print_output:
    input:
        data = "{mito}/{file}_{rmTrans}.Rdata"
    output:
        txt = "{mito}/{file}_{rmTrans}.txt"
    params:
        path = config["print_contamMix"]
    shell:
        "Rscript {params.path} {input.data} > {output.txt}"