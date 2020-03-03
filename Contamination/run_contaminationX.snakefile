configfile: "multiple_purposes.yaml"

minDepth = [1,2,3] 

#=============================================================================#
# Parse samples and libraries per sample
samples = list(config["Contamination"]["contaminationX"]["Samples"].keys())
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
    sample_list = list(config["Contamination"]["contaminationX"]["Samples"].keys())
    samples = ["{}/{}.{}_{}".format(sample, sample, geno_name, lib) for sample in config["Contamination"]["contaminationX"]["Samples"].keys() for lib in sample_libs[sample]]
    output = samples
    return(output)

all_prefix = output_prefix()

rule all:
    input:
        result = expand('X/{file}_md{depth}.result',
                        file = all_prefix, depth = minDepth)

rule X_bam_library:
    input:
        bam = "{file}.bam"
    output:
        bam = "X/{file}_{library}.bam"
    wildcard_constraints:
        library = "("+"|".join([lib for lib in set( [lib for myList in sample_libs.values() for lib in myList])])+")"
    shell:
        """
        if [ {wildcards.library} == "All" ]
        then
            samtools view -b {input.bam} X > {output.bam} 
        else
            samtools view -l {wildcards.library} -b {input.bam} X > {output.bam} 
        fi
        """

rule doCounts:
    input:
        bam = protected("{file}.bam"),
        bai = protected("{file}.bam.bai")
    output:
        icounts = "{file}.icnts.gz"
    shell:
        """
         angsd -i {wildcards.file}.bam -r X: -doCounts 1  \
         -iCounts 1 -minMapQ 30 -minQ 20 -out {wildcards.file}
        """

rule angsd_contamination:
    input:
        icounts = "{file}.icnts.gz",
        freqs = "/users/dcruzdav/Project_popgen/americas/install/contaminationX/HapMapFreqs/HapMapCEU.gz"
    output:
        counts = "{file}_md{depth}_counts"
    shell:
        """
        set +e ;
        contamination -b 5000000 -c 154900000 \
           -k 1 -m 0.05 -d {wildcards.depth} -e 20 \
           -h {input.freqs} \
           -a {input.icounts} > {output.counts} ;
        
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 0
        fi
        """


rule estimate:
    input:
        counts = "{file}_md{depth}_counts",
        freqs = "/users/dcruzdav/Project_popgen/americas/install/contaminationX/HapMapFreqs/HapMapCEU.gz"
    output:
        result = "{file}_md{depth}.result"
    threads:
        16
    shell:
        """
        nLines=$(wc -l {input.counts} |cut -f1 -d' ')
        if [ $nLines -eq 0 ] ; then
            touch {output.result}
        else
            Rscript ~/Project_popgen/americas/install/contaminationX/bin/ContaEstBoth.R \
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
