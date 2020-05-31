configfile: "multiple_purposes.yaml"
import os,glob,pysam, subprocess
include: "parse_resources.smk"

myPath = os.path.expanduser(config["Microbiome"]["bam2fastq"]["Samples"])
myFiles = glob.glob(myPath) 
samples = [f.split("/")[-2] for f in myFiles]

def view_header(sample):
    myCommand = "samtools view -H " + sample + "/" + sample + ".hg19_low_qual.bam > " + sample + ".header"
    os.system(myCommand)

[view_header(sample) for sample in samples]

paired_libs = {}

def add_paired(sample, lib):
    lib = lib.decode()
    if sample in paired_libs.keys():
        paired_libs[sample].append(lib)
    else:
        paired_libs[sample] = [lib]

def get_libraries(sample):
    with open(sample + ".header", "r") as h:
        lines = [line for line in h.readlines()]
    libs = list(set([line.split()[2].split(":")[1] for line in lines if line.split()[0] == "@RG"]))
    [add_paired(sample, l) for l in is_paired_end(sample)] 

    return(libs)

def is_paired_end(sample):
    myCommand =  'cat ' + sample + '.header' + "|grep '^@PG' |grep sampe "
    proc = subprocess.Popen(myCommand, stdout = subprocess.PIPE, shell = True)
    (out, err) = proc.communicate()
    
    myLibs = set([line.split(b"\\tLB:")[1].split(b"\\t")[0] for line in out.splitlines()])

    return(myLibs)
    

lib_samples = {}

def add_key(sample):
    libs = get_libraries(sample)
    lib_samples[sample] = libs

[add_key(sample) for sample in samples ]


rule all:
    input:
        stats = expand("Microbiome/stats/{sample}_{type}.txt", 
                    sample = samples, 
                    type = ["orig_bam","unmapped_bam", "unmapped_fastq"]),
        md5 = expand("Microbiome/stats/{sample}_md5.txt", sample = samples)

# rule merge_stats:
#     input:
#     output:
#         "Microbiome/stats/unmapped_reads.txt"

rule filter_bam:
    input:
        bam = "{sample}/{sample}.hg19_low_qual.bam"
    output:
        bam = "Microbiome/BAM_unmapped/{sample}/{library}/{sample}_{library}.bam"
    log:
        "Microbiome/logs/filter_bam_{sample}_{library}.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("filter_bam_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("filter_bam_time", attempt, 48),      
    shell:
        """
        samtools view -bf4 -l {wildcards.library} {input.bam} > {output.bam} 2>{log}
        """

def paired(sample, lib):
    if sample in paired_libs.keys() and lib in paired_libs[sample]:
        return("paired")
    else:
        return("single")
    
def expand_fastq(lib, sample):
    if sample in paired_libs.keys() and lib in paired_libs[sample]:
        return(["_R1.fastq.gz", "_R2.fastq.gz"])
    else:
        return([".fastq.gz"])

fastq_sample_lib = {}
def add_fastqs(sample, lib):
    fastqs =  expand_fastq(lib, sample)
    myPrefix = "Microbiome/FASTQ_unmapped/{sample}/{library}/{sample}_{library}".format(sample = sample, library = lib)

    myInput = "{sample}/{sample}.hg19_low_qual.bam".format(sample = sample)
    fastq_sample_lib[sample] = myInput
    # myDict = {"sufix": fastqs, "input": myInput}
    # fastq_sample_lib[myPrefix] = myDict

    # if sample in fastq_sample_lib.keys():
    #     myD = fastq_sample_lib[sample]
    #     myD[lib] = ["Microbiome/FASTQ_unmapped/{sample}/{library}/{sample}_{library}{end}".format(sample = sample, library = lib, end = end) for end in fastqs]
    #     fastq_sample_lib[sample] = myD
    # else:
    #     fastq_sample_lib[sample] = {lib : ["Microbiome/FASTQ_unmapped/{sample}/{library}/{sample}_{library}{end}".format(sample = sample, library = lib, end = end) for end in fastqs]}

[add_fastqs(sample, lib) for sample in samples for lib in lib_samples[sample] ]

rule bam2fastq:
    input:
        bam =lambda wildcards: fastq_sample_lib[wildcards.sample]  # "{sample}/{sample}.hg19_low_qual.bam"
    output:
        fastq_gz = expand("Microbiome/FASTQ_unmapped/{sample}/{library}/{sample}_{library}{end}",
                             sample = "{sample}", 
                             library = "{library}",
                             end = ["{end}"]) #expand_fastq(lib="{library}", sample="{sample}"))
    log:
        "Microbiome/logs/bam2fastq_{sample}_{library}.{end}.txt"
    params:
        lib_type = lambda wildcards: paired(wildcards.sample, wildcards.library)
    wildcard_constraints:
        end = "(_R1.fastq.gz|_R2.fastq.gz|.fastq.gz)",
    threads: 4
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("bam2fastq_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("bam2fastq_time", attempt, 48), 
    shell:
        """
        if [ {params.lib_type} == 'paired' ]
        then
            # samtools view -u -f 12 -F 304 -b {input.bam} | \
            # samtools sort -n - | \
            # samtools fastq  -1 {output.fastq_gz}[0] -2 {output.fastq_gz}[1] - 
            touch {output}     2>{log}
        else
            samtools view -bl  {wildcards.library} -f4 {input.bam} | \
            samtools fastq -c9 - |gzip >{output.fastq_gz} 2>{log}
        fi
        """

rule headers:
    input:
        bam = "{sample}/{sample}.hg19_low_qual.bam"
    output:
        header = temp("{sample}.header")
    shell:
        """
        samtools view -H {input.bam} > {output.header} 2>{log}
        """

def count_unmapped(file):
    myCommand = "samtools view -cf4 " + file
    print("command count_unmapped sayssss " + myCommand)
    proc = subprocess.Popen(myCommand, stdout = subprocess.PIPE, shell = True)
    (out, err) = proc.communicate()
    unmapped = str(out)
    return(unmapped)

def write_stats(sample, lib, type, file_out):

    if type == "BAM_unmapped_original":
        if lib == "all":
            sample_lib = sample + "/" + sample + ".hg19_low_qual.bam" 
        else:
            sample_lib = "-l " + lib + " " + sample  + "/" + sample + ".hg19_low_qual.bam" 
        unmapped = str(count_unmapped(sample_lib))
        
    elif lib == "all":
        sample =  sample + "/" + sample + ".bam" 
        unmapped = count_unmapped(sample)
    else:
        sample_lib = "-l " + lib + " " + "Microbiome/BAM_unmapped/" + sample + "/" + lib + "/" + sample + "_" + lib + ".bam" 
        unmapped = str(count_unmapped(sample_lib))
    file_out.write("\t".join([sample, lib, unmapped, type]) + "\n")

rule get_stats_bam_orig:
    input:
        bam = "{sample}/{sample}.hg19_low_qual.bam"
    output:
        stats = "Microbiome/stats/{sample}_orig_bam.txt"
    log:
        "Microbiome/logs/get_stats_bam_orig_{sample}.txt"
    run:
        libs = lib_samples[wildcards.sample]
        with open(output.stats, 'w') as file_out:
            write_stats(wildcards.sample, "all", "BAM_unmapped_original", file_out)
            [write_stats(wildcards.sample, l, "BAM_unmapped_original", file_out) for l in libs]

def return_libs(sample):
    myLibs = lib_samples[sample]
    return(myLibs)

rule get_stats_bam_unmapped:
    input:
        bam = lambda wildcards: expand("Microbiome/BAM_unmapped/{sample}/{library}/{sample}_{library}.bam",
                     sample = "{sample}", library = return_libs(wildcards.sample))
    output:
        stats = "Microbiome/stats/{sample}_unmapped_bam.txt"
    log:
        "Microbiome/logs/get_stats_bam_unmapped_{sample}.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("stats_unmapped_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("stats_unmapped_time", attempt, 12), 
    run:
        libs = lib_samples[wildcards.sample]
        print("libs are" + str(len(libs)))
        with open(output.stats, 'w') as file_out:
            print("output issss " + output.stats )
            [write_stats(wildcards.sample, l, "BAM_unmapped_diana",file_out) for l in libs]

def expand_path_fastq(sample):
    libs = lib_samples[sample]
    path = [ l + "/" + sample + "_" + l + end for l in libs for end in expand_fastq(l, sample)]
    return(path)

rule get_stats_fastq_unmapped:
    input:
        fastq_gz = lambda wildcards: expand("Microbiome/FASTQ_unmapped/{sample}/{fastq}",
                    sample = "{sample}",
                    fastq = expand_path_fastq(sample=wildcards.sample))
    output:
        stats = "Microbiome/stats/{sample}_unmapped_fastq.txt"
    log:
        "Microbiome/logs/get_stats_fastq_unmapped_{sample}.txt"
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("stats_unmapped_mem", attempt, 2),
        runtime = lambda wildcards, attempt: get_runtime_alloc("stats_unmapped_time", attempt, 12), 
    run:   

        libs = lib_samples[wildcards.sample]

        def write_unmapped(f, file_out):
            myCommand = "gunzip -c " + f + " |wc -l"
            # sample = f.split("/")[2]
            lib = f.split("/")[3]
            print("COUNTING FASTQ " + myCommand)
            proc = subprocess.Popen(myCommand, stdout = subprocess.PIPE, shell = True)
            (out, err) = proc.communicate()
            unmapped = str(int(out)/4)
            file_out.write("\t".join([wildcards.sample, lib, unmapped, "FASTQ_unmapped_diana"]) + "\n")
        
        with open(output.stats, "w") as file_out:
            [write_unmapped(f, file_out) for f in input.fastq_gz]

rule md5:
    input:
        unmapped = "Microbiome/{dir}/{sample}/{library}/{sample}_{library}{extension}"
    output:
        md5 = temp("Microbiome/{dir}/{sample}/{library}/{sample}_{library}{extension}.md5sum")
    wildcard_constraints:
        extension = "(.bam|.fastq.gz|_R1.fastq.gz|_R2.fastq.gz)"
    log:
        "Microbiome/logs/md5_{dir}_{sample}_{library}_{extension}.txt"
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("md5_mem", attempt, 2),
        runtime = lambda wildcards, attempt: get_runtime_alloc("md5_time", attempt, 24), 
    shell:
        """
        md5sum {input.unmapped} > {output.md5} 2>{log}
        """
    
def expand_md5_fastq(sample):
    myMD5 = ["Microbiome/FASTQ_unmapped/{sample}/{library}/{sample}_{library}{extension}.md5sum".format(sample = sample,
    library = lib, extension = ext) for lib in lib_samples[sample] for ext in expand_fastq(sample=sample, lib = lib)]
    return(myMD5)

rule wrap_md5:
    input:
        md5_bam = lambda wildcards: expand("Microbiome/BAM_unmapped/{sample}/{library}/{sample}_{library}.bam.md5sum",
                        sample = "{sample}", library = lib_samples[wildcards.sample]),
        md5_fastq_gz = lambda wildcards: expand_md5_fastq(wildcards.sample)
        # md5_fastq_gz = lambda wildcards: expand("Microbiome/FASTQ_unmapped/{sample}/{library}/{sample}_{library}{extension}.md5sum",
        #                 zip, sample = "{sample}",library = lib_samples[wildcards.sample],
        #                 extension = [expand_fastq(sample=wildcards.sample, lib = l) for l in lib_samples[wildcards.sample]])
    output:
        md5 = "Microbiome/stats/{sample}_md5.txt"
    log:
        "Microbiome/logs/wrap_md5_{sample}.txt"
    wildcard_constraints:
        extension = "(.fastq.gz|_R1.fastq.gz|_R2.fastq.gz)"

    run:
        libs = lib_samples[wildcards.sample]

        def write_md5(md5_to_open, #sample, lib, type, 
        file_out):
            sample = md5_to_open.split("/")[2]
            lib = md5_to_open.split("/")[3]
            type = md5_to_open.split("/")[1] + "_diana"
            # f = str_to_format.format(sample = sample, library = lib)

            with open(md5_to_open, "r") as file:
                line = file.readline()
                md5 = line.split()[0]
                file_out.write("\t".join([sample, lib, md5, type]) + "\n")
                    
        with open(output.md5, "w") as file_out:
            [write_md5(md5_to_open = bam, file_out = file_out) for bam in input.md5_bam]
            [write_md5(md5_to_open = fastq, file_out = file_out) for fastq in input.md5_fastq_gz]
            # myString = "Microbiome/BAM_unmapped/{sample}/{library}/{sample}_{library}.bam.md5sum"
            # [write_md5(myString, wildcards.sample, l, "BAM_unmapped_diana", file_out) for l in libs]
            # myString = "Microbiome/FASTQ_unmapped/{sample}/{library}/{sample}_{library}.fastq.gz.md5sum"
            # [write_md5(myString, wildcards.sample, l, "FASTQ_unmapped_diana", file_out) for l in libs]

