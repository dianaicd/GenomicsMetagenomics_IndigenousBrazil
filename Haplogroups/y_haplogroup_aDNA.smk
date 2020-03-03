configfile: "multiple_purposes.yaml"

#-----------------------------------------------------------------------------#
# Parse samples
groups = [group for group in config["Haplogroups"]["y"].keys()]
groups.remove("path_db")
groups.remove("path_goQuery")

samples = {}
for group in groups:
    samples[group] = {sample:config["Haplogroups"]["y"][group][sample]  
    for sample in config["Haplogroups"]["y"][group] }

#-----------------------------------------------------------------------------#

rule all:
    input:
        # ancestral = ["{group}/{file}.1000Y.ancestral.txt".format(group = group, file = file) for group in groups for file in samples[group].keys()],
        # derived = ["{group}/{file}.1000Y.derived.txt".format(group = group, file = file) for group in groups for file in samples[group].keys()]
        haplogroups = [ "Results/output/{group}/_{group}.1000Y.txt".format(group = group) for group in groups]

rule haploid_call:
    input:
        bam = lambda wildcards: samples[wildcards.group][wildcards.file]
    output:
        haplo = "data/{group}/{file}.haplo.gz"
    params:
        basename = "data/{group}/{file}",
        baseQ = config["BaseQuality"] if "BaseQuality" in config.keys() else 20
    shell:
        """
        angsd -i {input.bam} -dohaplocall 1 -doCounts 1 -r Y: -out {params.basename} -maxMis 0 -minQ {params.baseQ}
        """

rule format_genos:
    input:
        haplo = "data/{group}/{file}.haplo.gz"
    output:
        genos = "data/{group}/{file}.1000Y.genos.txt.bz2"
    shell:
        """
        gunzip -c {input.haplo} | cut -f 2,3  |tail -n+2| bzip2 > {output.genos}
        """

rule gather_databases:
    input:
        index = config["Haplogroups"]["y"]["path_db"] + "1000Y.index2label.pkl",
        db = config["Haplogroups"]["y"]["path_db"] + "1000Y.snp.db.pkl"
    output:
        index = "data/1000Y.index2label.pkl",
        db = "data/1000Y.snp.db.pkl"
    shell:
        """
        ln -s {input.index} {output.index}
        ln -s {input.db} {output.db}
        """


rule query_haplogroup:
    input:
        index = "data/1000Y.index2label.pkl",
        db = "data/1000Y.snp.db.pkl",
        genos = lambda wildcards: expand("data/{group}/{file}.1000Y.genos.txt.bz2", group = wildcards.group,
                        file = samples[wildcards.group].keys() )
    output:
        "Results/output/{group}/_{group}.1000Y.txt"
    params:
        path_goQuery = config["Haplogroups"]["y"]["path_goQuery"]
    shell:
        """
        mkdir -p Results ; cd Results

        {params.path_goQuery} {wildcards.group}
        """