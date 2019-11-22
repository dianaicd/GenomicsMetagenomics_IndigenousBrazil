configfile: "multiple_purposes.yaml"
#include: "genotype_calling_all.smk"
#=============================================================================#
# Some defaults
default_roh = {"homozyg-window-kb ": 5000,
        "homozyg-window-snp": 50,
        "homozyg-window-het": 1,
        "homozyg-window-missing": 5,
        "homozyg-window-threshold": 0.05,
        "homozyg-snp": 100,
        "homozyg-kb": 1000,
        "homozyg-density": 50, 
        "homozyg-gap": 1000
        }

map_filter_name = {"window_kb":"homozyg-window-kb ",
        "window_snp": "homozyg-window-snp",
        "hets_window": "homozyg-window-het",
        "window_missing": "homozyg-window-missing",
        "window_threshold": "homozyg-window-threshold",
        "segment_snp": "homozyg-snp",
        "segment_kb":"homozyg-kb",
        "segment_density": "homozyg-density", 
        "homozyg_gap": "homozyg-gap"}

def build_filter(par_name):
    myFilter = ""
    if "filters" in config.keys():
        filters = list(config["roh_plink"]["filters"].keys())
        if par_name in filters:
            myFilter = "--" + map_filter_name[par_name] + " " + str(config["roh_plink"]["filters"][par_name])
        # else:
            #myFilter = "--" + map_filter_name[par_name] + " " + str(default_roh[map_filter_name[par_name]])
            
    return(myFilter)

ref_genome = config["ref_genome"]

with open(ref_genome + ".fai", 'r') as index:
    chromosomes = [line.split()[0] for line in index.readlines()]
    chromosomes = [str(x) for x in range(1,23)]

bamlists = [l for l in list(config["roh_plink"]["bamlists"].keys())]

#=============================================================================#

rule all:
    input:
        roh = expand("ROH/{bamfile}_{chr}_{rmTrans}_roh.hom",
                     bamfile = bamlists, chr = chromosomes, rmTrans = ["all", "rmTrans"])


rule roh_plink:
    input:
        bed = "Filtered/{bamfile}_{chr}_depth_filter_{rmTrans}.bed"
    output:
        roh = "ROH/{bamfile}_{chr}_{rmTrans}_roh.hom"
    params:
        window_snp = build_filter("window_snp"),
        hets_window = build_filter("hets_window"),
        window_missing = build_filter("window_missing"),
        window_threshold = build_filter("window_threshold"),
        segment_density = build_filter("segment_density"),
        segment_snp = build_filter("segment_snp"),
        segment_kb = build_filter("segment_kb"),
        homozyg_gap = build_filter("homozyg_gap")
    shell:
        """
        plink --homozyg \
        --bfile Filtered/{wildcards.bamfile}_{wildcards.chr}_depth_filter_{wildcards.rmTrans} \
        --out ROH/{wildcards.bamfile}_{wildcards.chr}_{wildcards.rmTrans}_roh \
        {params.window_snp} \
        {params.hets_window} \
        {params.window_missing} \
        {params.window_threshold} \
        {params.segment_density} \
        {params.segment_snp} \
        {params.segment_kb} \
        {params.homozyg_gap}
        """