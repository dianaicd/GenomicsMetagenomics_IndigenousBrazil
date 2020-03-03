configfile: "multiple_purposes.yaml"

import subprocess
mito=config["Contamination"]["contamMix"]["Mitochondrial"]

wildcard_constraints:
    mito = mito

#=============================================================================#
# Parse samples and libraries per sample
samples = list(config["Contamination"]["contamMix"]["Samples"].keys())
geno_name = config["Contamination"]["Whole_genome"]

def parse_libraries(sample, geno_name):
    myBam = sample + "/" + sample + "." + geno_name + ".bam"
    myCommand = "samtools view -H " + myBam
    proc = subprocess.Popen(myCommand, stdout = subprocess.PIPE, shell = True)
    (out, err) = proc.communicate()
    
    libs = list(set([line.split()[2].split(":")[1]
     for line in out.decode("utf-8").split("\n") 
     if len(line.split()) and line.split()[0] == "@RG"]))
    return(libs)

def add_libs(myDict, sample, geno_name):
    myDict[sample] = parse_libraries(sample, geno_name) + ["All"]

sample_libs = {}
[add_libs(sample_libs, sample, geno_name) for sample in samples]
#-----------------------------------------------------------------------------#
# Function to get paths to BAM files for samples
def output_prefix():
    sample_list = list(config["Contamination"]["contamMix"]["Samples"].keys())
    samples = ["{}/{}.{}_{}".format(sample, sample, geno_name, lib) for sample in config["Contamination"]["contamMix"]["Samples"].keys() for lib in sample_libs[sample]]
    output = samples
    return(output)

all_prefix = output_prefix()
# print(all_prefix)
rule all:
    input:
        result_data = expand('{mito}/{file}_{rmtrans}.{extension}',
                         mito = mito, file = all_prefix, 
                         rmtrans = ['all', 'rmTrans'], extension = ['pdf', 'Rdata']),
        result_txt = expand("{mito}/{file}_{rmTrans}.txt", 
                            mito = mito, file = all_prefix, rmTrans = ['all', 'rmTrans']),
        merged_est = mito + "/estimates.txt"


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

ruleorder: mito_bam_library > index > contammix

rule mito_bam_library:
    input:
        bam = "{file}.bam"
    output:
        bam = "{mito}/{file}_{library}.bam"
    wildcard_constraints:
        library = "("+"|".join([lib for lib in set( [lib for myList in sample_libs.values() for lib in myList])])+")"
    shell:
        """
        if [ {wildcards.library} == "All" ]
        then
            samtools view -b {input.bam} {mito} > {output.bam} 
        else
            samtools view -l {wildcards.library} -b {input.bam} {mito} > {output.bam} 
        fi
        """

rule index:
    input:
        "{x}.bam"
    output:
        "{x}.bam.bai"
    shell:
        "samtools index {input}"

rmTrans={"all":"", "rmTrans":"--transverOnly"}

rule print_output:
    input:
        data = "{file}_{rmTrans}.Rdata"
    output:
        txt = "{file}_{rmTrans}.txt"
    params:
        path = config["Contamination"]["contamMix"]["print_contamMix"]
    shell:
        "/software/R/3.5.1/bin/Rscript {params.path} {input.data} > {output.txt}"

rule contammix:
    input:
        bam="{file}_consensus.bam",
        bai="{file}_consensus.bam.bai", 
        maln="{file}_311_aligned.fasta" 
    output:
        fig=protected("{file}_{rmTrans}.pdf"),
        data=protected("{file}_{rmTrans}.Rdata")
    params:
        transitions=lambda wildcards: rmTrans.get(wildcards.rmTrans),
        nIter=config["Contamination"]["contamMix"]["nIter"],
        path=config["Contamination"]["contamMix"]["path"],
        prefix="{file}_{rmTrans}",
        basename="{file}_{rmTrans}"
    log:
         "logs/{file}_{rmTrans}_contammix.log"
    shell:
        """
        /software/R/3.5.1/bin/Rscript {params.path} --samFn {input.bam} --nIter {params.nIter} \
        --malnFn {input.maln} --alpha 0.1 --figure {output.fig} \
        --saveData {params.basename} {params.transitions}
        """

rule realign:
    input:
        bam = "{file}.bam", 
        bai = "{file}.bam.bai", 
        fasta = "{file}.fa" 
    output:
        low_qual = "{file}_low_qual.bam",
        bam="{file}_consensus.bam",
        fastq="{file}_consensus.truncated.gz"
    params:
        basename="{file}_consensus",
        q = 30
    log:
        "logs/{file}_consensus.log"
    shell:
        """
        bedtools bamtofastq -i {input.bam} -fq {params.basename}.truncated ;
        gzip {params.basename}.truncated ;
        
        bwa aln -l 1024 -t 1 {input.fasta} -f \
            {params.basename}.sai {params.basename}.truncated.gz;
        
        bwa samse -n 3 \
        -r "@RG\\tID:{params.basename}\\tLB:{params.basename}\\tSM:{params.basename}\\tPL:Illumina"   \
        {input.fasta} {params.basename}.sai {params.basename}.truncated.gz \
        | samtools view -Sb > {params.basename}_full.bam;
        
        samtools sort --threads 1 {params.basename}_full.bam > {params.basename}_sorted.bam

        samtools view -b --threads 1 -F 4 -q {params.q} -U {output.low_qual} \
        {params.basename}_sorted.bam > {params.basename}_filtered.bam 2>> {log}
        samtools index {params.basename}_filtered.bam
        # realign
        GenomeAnalysisTK -T RealignerTargetCreator -I {params.basename}_filtered.bam \
         -R {input.fasta}  -o {params.basename}.intervals 2>>{log}
        
        GenomeAnalysisTK -T IndelRealigner -I {params.basename}_filtered.bam \
         -R {input.fasta}  -targetIntervals {params.basename}.intervals \
         -o {params.basename}_realign.bam 2>>{log}

        samtools index {params.basename}_realign.bam 
        
        # calmd
        samtools calmd --threads 1 {params.basename}_realign.bam \
            {input.fasta} 2>> {log} | samtools view -bS - > {output.bam}

        """

rule maln:
    input:
        "{file}.fa"
    output:
        "{file}_311_aligned.fasta"
    params:
        mt311=config["Contamination"]["contamMix"]["contaminants"]
    shell:
        "cat {params.mt311} {input} > {wildcards.file}_311.fasta ;"
        "mafft --auto {wildcards.file}_311.fasta > {output} ;"


rule consensus_sample:
    input:
        bam = "{file}.bam", 
        bai = "{file}.bam.bai" 
    output:
        fasta = "{file}.fa", 
        dict = "{file}.dict",
        # bam = "{file}.bam" 
    params:
        basename = "{file}"
    shell:
        """
          mkdir -p {params.basename} ;
          angsd -doFasta 2 -doCounts 1 -i {input.bam} -out {params.basename} ;
          gunzip {output.fasta}.gz ;
          bwa index  {output.fasta} ;
          samtools faidx {output.fasta} ;
          picard-tools CreateSequenceDictionary REFERENCE={output.fasta} OUTPUT={output.dict} ;
        """

rule parse_estimates:
    input:
        est = "{file}_{rmTrans}.txt"
    output:
        est = temp("{file}_{rmTrans}.estimate")
    shell:
        """
        first_str=$(echo {input.est}| sed 's|.*/|| ; s/.hg19_/,/ ; s/_/,/ ; s/.txt//')
        up_low=$(head -n9 {input.est} |tail -n1 | sed 's/ /,/')
        est=$(head -n10 {input.est} |tail -n1 | cut -f3 -d ' ')
        echo "${{first_str}},${{est}},${{up_low}}" > {output.est}
        """

rule merge_estimates:
    input:
        est= expand("{mito}/{file}_{rmTrans}.estimate", 
                            mito = mito, file = all_prefix, rmTrans = ['all', 'rmTrans'])
    output:
        est = mito+"/estimates.txt"
    run:
        def write_est(in_file, out_file):
            with open(in_file, 'r') as f:
                line = f.readline()
                out_file.write(line)

        with open(output.est, 'w') as out:
            [write_est(in_file = file, out_file = out) for file in input.est]

