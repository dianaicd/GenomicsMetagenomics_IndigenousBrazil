configfile: "bcftools_stats.yaml"
import yaml
import pandas as pd

config_filters = yaml.load(open("genotype_calling.yaml", "r"))["geno_calls"]["filters"]

chromosomes = [str(i) for i in range(1, 23)]
individuals = config["individuals"]
tables = ["SN", "TSTV", "QUAL", "ST", "DP"]
#-----------------------------------------------------------------------------#

sufixes = {
    "mask": "mask{mask_label}".format(mask_label=config_filters["mask"]["label"]),
    "depth": "dp",
    "gq": "GQ{gq}".format(gq=config_filters["gq"]["min"]), 
    "repeats": "repeats{repeats_label}".format(repeats_label=config_filters["repeats"]["label"])
    }

filters_order = config_filters["order"]

def build_sufix(filters, filters_order):
    sufix_file = "_".join([sufixes[f] for f in filters_order])
    return sufix_file

filters_combinations = [
    build_sufix(list(config_filters.keys()), filters_order[0:i]) 
    for i in range(0, len(filters_order)+1)
    ]

filters_combinations[0] = "raw"
filters_combinations = ["raw", "maskStrict_repeatsUCSC_dp_GQ30"]
print(filters_combinations)
print(filters_order)

# specific to tables rearrangements
headers_tables = {
    "SN": ["SN", "id", "key", "value", "sample", "filter"],
    "TSTV": ["TSTV", "id", "ts", "tv", "ts_tv", "ts_1st_Alt", "tv_1st_ALT", "ts_tv_1st_ALT", "sample", "filter"],
    "QUAL": ["QUAL", "id", "Quality", "SNPs", "ts", "tv", "indels", "sample", "filter"],
    "ST": ["ST", "id", "type", "count", "sample", "filter"],
    "DP": ["DP", "id", "bin", "genotypes", "frac_genos", "sites", "frac_sites", "sample", "filter"]
}

columns_group = {
    "SN": ["sample", "filter", "key"],
    "TSTV": ["sample", "filter"],
    "QUAL": ["Quality", "sample", "filter"],
    "ST": ["type", "sample", "filter"],
    "DP": ["bin", "sample", "filter"]
}

columns_fraction = {
    "SN": {},
    "TSTV": {
        "ts_tv":{"num":"ts", "div":"tv"}, 
        "ts_tv_1st_ALT":{"num":"ts_1st_Alt", "div":"tv_1st_ALT"}
        },
    "QUAL": {},
    "ST": {},
    "DP": {}
}

def fix_fraction(table, column, short_table):
    num = short_table[columns_fraction[table][column]["num"]]
    div = short_table[columns_fraction[table][column]["div"]]
    short_table[column] = num / div 

#-----------------------------------------------------------------------------#
rule all:
    input:
        stats = expand(
            "stats/{filter}/chr{chr}/{ind}_chr{chr}.txt", 
            ind = individuals, 
            chr = chromosomes,
            filter = filters_combinations
        ),
        tables_ind_filter = expand( 
            "concat_stats/{table}/{filter}/{ind}.txt", 
            table = tables,
            filter = filters_combinations,
            ind = individuals
            ),
        tables_merged = expand(
            "ind_stats/{table}.txt",
            table = tables
        )

rule bcftools_stats_raw:
    input:
        bcf = "Raw/{ind}_chr{chr}.bcf"
    output:
        stats = "stats/raw/chr{chr}/{ind}_chr{chr}.txt"
    log:
        "logs/stats_{ind}_raw_{chr}.txt"
    resources:
        memory = 4*1024,
        runtime = 10
    shell:
        """
        bcftools stats {input.bcf} > {output.stats} 
        """

rule bcftools_stats_filter:
    input:
        bcf = "Filtered/chr{chr}/{ind}_{filter}_chr{chr}.bcf"
    output:
        stats = "stats/{filter}/chr{chr}/{ind}_chr{chr}.txt"
    resources:
        memory = 4*1024,
        runtime = 10
    log:
        "logs/stats_{ind}_{filter}_{chr}.txt"
    shell:
        """
        bcftools stats {input.bcf} > {output.stats} 
        """

# Tables I would like to get:
# "SN", "TSTV", "QUAL", "ST", "DP"
# columns to add:
# ind, chr, filter

rule concat_tables:
    input:
        stats = expand(
            "stats/{filter}/chr{chr}/{ind}_chr{chr}.txt", 
            chr = chromosomes,
            filter = "{filter}",
            ind = "{ind}"
            )
    output:
        stats = "concat_stats/{table}/{filter}/{ind}.long.txt"
    resources:
        memory = 1*1024,
        runtime = 10
    log:
        "logs/stats_{ind}_{filter}_{table}.txt"
    shell:
        """
        for file in {input.stats}
        do
            grep -P "^{wildcards.table}" $file 
        done | \
            awk 'BEGIN {{OFS="\\t"}} {{print $0,"{wildcards.ind}","{wildcards.filter}"}}' > {output.stats} &2>{log}
        """


rule fix_tables:
    input:
        stats = "concat_stats/{table}/{filter}/{ind}.long.txt"
    output:
        stats = "concat_stats/{table}/{filter}/{ind}.txt"
    resources:
    log:
    run:
        filename = input.stats
        table = wildcards.table

        long_table = pd.read_csv(
                filename, 
                delimiter="\t", 
                header=None, 
                names=headers_tables[table]
                )

        short_table = long_table.groupby(columns_group[table]).sum()
        [fix_fraction(table, column, short_table) for column in columns_fraction[table].keys()]
        if table == "DP":
            short_table["frac_sites"] = short_table["sites"]/short_table["sites"].sum()

        short_table.to_csv(output.stats)

rule concat_inds_tables:
    input:
        stats = expand(
            "concat_stats/{table}/{filter}/{ind}.txt",
            ind = individuals,
            table = "{table}",
            filter = filters_combinations
        )
    output:
        stats = "ind_stats/{table}.txt"
    resources:
    log:
    run:
        stats = pd.concat([pd.read_csv(filename) for filename in input.stats ], axis = 0)
        stats.to_csv(output.stats)