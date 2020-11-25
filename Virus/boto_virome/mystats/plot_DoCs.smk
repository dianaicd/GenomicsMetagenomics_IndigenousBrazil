# Get depth of coverage plots

configfile: 'config_DoC.yaml'

def strip_jump(line):
    line = str(line)
    line = line.rstrip('\n')
    return(line)

with open(config['plot_names'], 'r') as plot_names_file:
    plot_names = [strip_jump(line) for line in plot_names_file.readlines()]

rule all:
    input:
        expand("plots_DoC/{plot}", plot=plot_names)

rule plot_DoC:
    input:
        "{sample}/coverage_tables/quality_{quality}/{sample}.{virus}.mapq{quality}.coverage.tsv"
    output:
        "plots_DoC/{sample}.{virus}.mapq{quality}.coverage.png"
    log:
        "plots_DoC/logs/{sample}.{virus}.mapq{quality}.coverage.log"
    shell:
        "Rscript --vanilla plot_DoC.R --cov_table {input} --plot_name {output} 2> {log}"