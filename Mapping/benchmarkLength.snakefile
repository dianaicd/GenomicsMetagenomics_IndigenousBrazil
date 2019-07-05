configfile: "2019_07_01.yaml"
samples=config["samples"].split()
nRep = 10
rule all:
    input:
        myScript=expand("{sample}/{sample}.myscript.length", sample = samples),
        myPython=expand("{sample}/{sample}.python.length", sample = samples),
        theirDepth=expand("{sample}/{sample}.their.depth", sample = samples),
        paleomixCov=expand("{sample}/{sample}.coverage", sample = samples),
        paleomixDepth=expand("{sample}/{sample}.depths", sample = samples)


rule perlLength:
    input:
        "{file}.bam"
    output:
        "{file}.myscript.length"
    benchmark:
        repeat("{file}/benchmarks_myPerl.txt", nRep)
    shell:
        "samtools view {input} |"
        "perl ~/Projects/Botocudos/Scripts/DataQuality/length.pl -o {output} "

rule pythonLength:
    input:
        "{file}.bam"
    output:
        "{file}.python.length"
    benchmark:
        repeat("{file}/benchmarks_python.txt", nRep)
    shell:
        "cat {input} |"
        "python ~/Projects/Botocudos/Scripts/DataQuality/read_length.py -o {output} "

rule theirDepth:
    input:
        "{file}.bam"
    output:
        "{file}.their.depth"
    benchmark:
        repeat("{file}/benchmarks_theirdepth.txt", nRep)
    conda:
        "myEnv_python2.7.yaml"
    shell:
        #"conda init bash ; conda activate python27 ; "
        "python ~/install/sapfo/depth.py {input} > {output} ;"


rule paleomixCov:
    input:
        "{file}.bam"
    output:
        "{file}.coverage"
    benchmark:
        repeat("{file}/benchmarks_paleomixCov.txt", nRep)
    shell:
        "paleomix coverage {input} --overwrite"

rule paleomixDp:
    input:
        "{file}.bam"
    output:
        "{file}.depths"
    benchmark:
        repeat("{file}/benchmarks_paleomixDp.txt", nRep)
    shell:
        "paleomix depths {input} --overwrite "