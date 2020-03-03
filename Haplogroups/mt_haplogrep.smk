configfile: "multiple_purposes.yaml"

import os,glob

def expand_path():
    paths = list(config["Haplogroups"]["mt"]["Samples"].values())
    full_paths = [os.path.expanduser(p) for p in paths]
    bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

mito = config["Haplogroups"]["mt"]["chromosome"] if "chromosome" in config["Haplogroups"]["mt"].keys() else "MT"
baseQ = config["Haplogroups"]["mt"]["baseQ"] if "baseQ" in config["Haplogroups"]["mt"].keys() else 20
mapQ = config["Haplogroups"]["mt"]["mapQ"] if "mapQ" in config["Haplogroups"]["mt"].keys() else 30
ref = config["ref_genome"]

bams = expand_path()
sample_file = {}

def add_key(sample, file):
    sample_file[sample] = file

samples = [file.split("/")[-1].replace(".hg19.bam", "") for file in bams]
[add_key(samples[i], bams[i]) for i in range(0, len(bams))]

rule all:
    input:
        fasta = expand("{mito}/{sample}/{sample}_{mito}_{type}.fa", 
                        mito = mito, sample = samples, type = ["all", "rmTrans"]),
        haplo = expand("{mito}/{sample}/{sample}_{mito}_{type}.haplo",
                         mito = mito, sample = samples, type = ["all", "rmTrans"])

rule index:
    input:
        "{x}.bam"
    output:
        "{x}.bam.bai"
    shell:
        "samtools index {input}"

rule mk_dir:
    output:
        "{mito}/{sample}/"
    shell:
        """
        mkdir {mito}/{wildcards.sample}
        """

rmTrans = {"rmTrans": 1, "all": 0}

rule filter_bam:
    input:
        bam = lambda wildcards: sample_file[wildcards.sample],
        bai = lambda wildcards: sample_file[wildcards.sample] + ".bai"
    output:
        bam = "{mito}/{sample}/{sample}_{mito}.bam"
    shell:
        """
        samtools view -b {input.bam} {mito} > {output.bam} ;
        """

rule consensus_fasta:
    input:
        bam = "{mito}/{sample}/{sample}_{mito}.bam",
        bai = "{mito}/{sample}/{sample}_{mito}.bam.bai"
    output:
        fasta = "{mito}/{sample}/{sample}_{mito}_{type}.fa", 
    params:
        baseQ = baseQ,
        mapQ = mapQ,
        ref = ref,
        rmTrans = lambda wildcards: rmTrans[wildcards.type]
    shell:
        """
          angsd -doFasta 2 -doCounts 1 -i {input.bam} \
            -out {mito}/{wildcards.sample}/{wildcards.sample}_{mito}_{wildcards.type} \
            -minMapQ {params.mapQ} -minQ {params.baseQ} -ref {params.ref} \
            -rmTrans {params.rmTrans} -r {mito} ;
          gunzip {output.fasta}.gz ;
        """

rule haplogrep:
    input:
        fasta = "{mito}/{sample}/{sample}_{mito}_{type}.fa", 
    output:
        haplo = "{mito}/{sample}/{sample}_{mito}_{type}.haplo", 
    shell:
        """
         java -jar ~/Project_popgen/americas/install/haplogrep-cmd/haplogrep-2.1.25.jar \
            --in {input.fasta} --format fasta --out {output.haplo} --extend-report
        """