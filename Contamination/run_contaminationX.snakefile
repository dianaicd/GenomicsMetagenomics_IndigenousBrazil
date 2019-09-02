configfile: "24samples.yaml"

minDepth = [1,2,3] 

def output_prefix():
    sample_list = list(config["Samples"].keys())
    samples = ["{}/{}".format(sample, sample) for sample in config["Samples"].keys()]
    libs = []
    for sample in sample_list:
        libs.append(["{sample}/{lib}/{lib}".format(sample = sample, lib = lib) 
                     for lib in list(config["Samples"][sample]["Libraries"].keys())])

    output = [item for sublist in libs for item in sublist] + samples

    return(output)

all_prefix = output_prefix()

rule all:
    input:
        result = expand('X/{file}_md{depth}.result',
                        file = all_prefix, depth = minDepth)


rule doCounts:
    input:
        bam = protected("{file}.bam"),
        bai = protected("{file}.bam.bai")
    output:
        icounts = "X/{file}.icnts.gz"
    shell:
        " ~/install/angsd/angsd -i {wildcards.file}.bam -r X: -doCounts 1 "
        " -iCounts 1 -minMapQ 30 -minQ 20 -out X/{wildcards.file}"


rule angsd_contamination:
    input:
        icounts = "X/{file}.icnts.gz",
        freqs = "/home/dcruzdva/install/contaminationX/HapMapFreqs/HapMapCEU.gz"
    output:
        counts = "X/{file}_md{depth}_counts"
    shell:
        " set +e ;"
        " ~/install/angsd/misc/contamination -b 5000000 -c 154900000 "
        "   -k 1 -m 0.05 -d {wildcards.depth} -e 20"
        "   -h {input.freqs}"
        "   -a {input.icounts} > {output.counts} ;"
        """
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        fi
        """


rule estimate:
    input:
        counts = "X/{file}_md{depth}_counts",
        freqs = "/home/dcruzdva/install/contaminationX/HapMapFreqs/HapMapCEU.gz"
    output:
        result = "X/{file}_md{depth}.result"
    threads:
        16
    shell:
        """
        nLines=$(wc -l {input.counts} |cut -f1 -d' ')
        if [ $nLines -eq 0 ] ; then
            touch {output.result}
        else
            Rscript ~/install/contaminationX/bin/ContaEstBoth.R \
            counts={input.counts} freqs={input.freqs} \
            maxsites=1000 nthr={threads} outfile={output.result} oneCns=1 
        fi
        """


rule index:
    input:
        "{x}.bam"
    output:
        "{x}.bam.bai"
    shell:
        "samtools index {input}"
