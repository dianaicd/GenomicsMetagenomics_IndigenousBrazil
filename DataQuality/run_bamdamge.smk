configfile: "multiple_purposes.yaml"
include: "parse_resources.smk"
import os,glob

def expand_path():
    paths = list(config["bamdamage"]["Samples"].values())
    full_paths = [os.path.expanduser(p) for p in paths]
    bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

all_paths = expand_path()
samples = [sample.split("/")[1].split(".")[0] for sample in all_paths]
ref = config["ref_genome"]

rule all:
    input:
        expand("bamdamage/{file}.dam_5prime.csv",
        file = [p.replace(".bam", "") for p in all_paths])

rule bamdamage:
    """
    Run bamdamage to quantify the deamination pattern
    """
    input:
        bam="{file}.bam",
        bai="{file}.bam.bai"
    output:
        dam="bamdamage/{file}.dam.pdf",
        length="bamdamage/{file}.length.pdf",
        length_table="bamdamage/{file}.length.csv",
        dam_5prime_table="bamdamage/{file}.dam_5prime.csv",
        dam_3prime_table="bamdamage/{file}.dam_3prime.csv"      
    log:
        "logs/bamdamage/{file}.log"
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("bamdamage_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("bamdamage_time", attempt, 24),        
    message: "--- BAMDAMAGE {input.bam}"
    params:
        bamdamage_params = config["bamdamage"]["bamdamage_params"] if "bamdamage_params" in config["bamdamage"].keys() else '',
        fraction = config["bamdamage"]["bamdamage_fraction"] if "bamdamage_fraction" in config["bamdamage"].keys() else 0,
    shell:    
    	"""
    	nb=$(samtools idxstats {input.bam} | awk '{{sum += $3}} END {{print sum}}'); 
    	nth_line=1; 
    	if [ {params.fraction} -eq 0 ]; then
    		nth_line=1; 
    	elif [ {params.fraction} -lt 1 ]; then
    		nth_line=$(( {params.fraction} * $nb )); 
    	elif [ {params.fraction} -lt "$nb" ]; then
    	   nth_line=$(( $nb / {params.fraction} )); 
    	fi;
        scripts/bamdamage {params.bamdamage_params} --nth_read $nth_line --output {output.dam} --output_length {output.length} {input.bam} 2> {log};

        """

