configfile: "multiple_purposes.yaml"

import os,glob,subprocess

def expand_path():
    paths = list(config["Contamination"]["schmutzi"]["Samples"].values())
    full_paths = [os.path.expanduser(p) for p in paths]
    bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

bams = expand_path()
samples = [file.split("/")[-1].replace(".hg19.bam", "") for file in bams]
mito = config["Contamination"]["schmutzi"]["Mitochondrial"] if "Mitochondrial" in config["Contamination"]["schmutzi"] else "MT"
ref_genome = config["ref_genome"]

rule all:
    input:
        endo_fasta = expand("{mito}/{file}/{file}_final_endo.fa",
                            mito = mito, file = samples)

rule query_mt_fasta:
    input:
        reference = ref_genome
    params:
        mito = mito
    output:
        mt_fa = "{mito}.fa"
    wildcard_constraints:
        mito = mito
    shell:
        """
        samtools faidx {input.reference} {params.mito} > {output.mt_fa}
        """

rule index_fasta:
    input:
        fasta = "{file}.fa"
    output:
        index = "{file}.fa.fai"
    shell:
        """
        samtools faidx {input.fasta}
        """

rule bam_index:
    input:
        bam = "{file}.bam"
    output:
        bai = "{file}.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

sample_file = {}

def add_key(sample, file):
    sample_file[sample] = file

[add_key(samples[i], bams[i]) for i in range(0, len(bams))]

rule bam_mito:
    input:
        bam = lambda wildcards: sample_file[wildcards.sample],
        bai =  lambda wildcards: sample_file[wildcards.sample] + ".bai",
    output:
        bam_mito = "{mito}/{sample}/{sample}.bam"
    params:
        mito = mito
    shell:
        """
        samtools view -b {input.bam} {params.mito} >{output.bam_mito}
        """

rule get_header:
    input:
        bam = "{mito}/{sample}/{sample}.bam",
        bai = "{mito}/{sample}/{sample}.bam.bai"
    output:
        header = temp("{mito}/{sample}/{sample}_header.txt")
    run:
        myCommand = "samtools view -H " + input.bam 
        proc = subprocess.Popen(myCommand, stdout = subprocess.PIPE, shell = True)
        (out, err) = proc.communicate()
        lines = [l for l in out.decode().split("\n")]

        def write_line_h(line, myOut):
            if len(line.split()) > 2:
                first_col = line.split()[0]
                second_col = line.split()[1]
                if first_col == "@SQ" and second_col == "SN:"+mito:
                    myOut.write(line + "\n")
                elif first_col != "@SQ":
                    myOut.write(line + "\n")

        with open(output.header, "w") as header:
            [write_line_h(line, header) for line in lines]


rule reheader:
    input:
        bam = "{mito}/{sample}/{sample}.bam",
        bai = "{mito}/{sample}/{sample}.bam.bai",
        header = "{mito}/{sample}/{sample}_header.txt",
    output:
        bam = "{mito}/{sample}/{sample}_reheaded.bam"
    shell:
        """
        (cat {input.header} <(samtools view {input.bam}) | samtools view -bo {output.bam} -)
        """

rule contDeam:
    input:
        mito_bam = "{mito}/{file}_reheaded.bam",
        mito_bai = "{mito}/{file}_reheaded.bam.bai",
        mito_ref = "{mito}.fa",
        mito_fai = "{mito}.fa.fai"
    wildcard_constraints:
        mito = mito
    output:
        cont_deam = "{mito}/{file}.cont.deam",
        endo_3p_prof = "{mito}/{file}.endo.3p.prof",
        endo_5p_prof = "{mito}/{file}.endo.5p.prof",
        config = "{mito}/{file}.deam.config",
        cont_est = "{mito}/{file}.cont.est",
    shell:
        """
        contDeam.pl --library single --out {wildcards.mito}/{wildcards.file} \
            --uselength --ref {input.mito_ref} {input.mito_bam}
        """

rule schmutzi:
    input:
        mito_ref = "{mito}.fa",
        mito_fai = "{mito}.fa.fai",
        mito_bam = "{mito}/{file}_reheaded.bam",
        mito_bai = "{mito}/{file}_reheaded.bam.bai",
        cont_deam = "{mito}/{file}.cont.deam",
        cont_est = "{mito}/{file}.cont.est",
        endo_5p_prof = "{mito}/{file}.endo.5p.prof",
        endo_3p_prof = "{mito}/{file}.endo.3p.prof"
    wildcard_constraints:
        mito = "("+mito+")",
        file = "(" + "|".join([i for i in samples]) + ")/(" +  "|".join([i for i in samples]) + ")"
    output:
        endo_fasta = "{mito}/{file}_final_endo.fa",
        cont_fasta = "{mito}/{file}_final_cont.fa",
        cont_est = "{mito}/{file}_final.cont.est"

    params:
        schmutzi_path = "~/Project_popgen/americas/install/schmutzi/share/schmutzi"
    threads:
        8
    shell:
        """
        schmutzi.pl --ref {input.mito_ref} --t {threads} \
            {wildcards.mito}/{wildcards.file} \
            {params.schmutzi_path}/alleleFreqMT/197/freqs/ \
            {input.mito_bam}
        """