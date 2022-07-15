configfile: "config_merge_vcfs.yaml"
groups=config["populations_to_merge"]

rule all:
    input:
        pop_sfs = [f"SFS/pop/{population}.sfs" for population in config["individuals_to_merge"].keys()],
        group_sfs = [f"SFS/groups/{group}.sfs" for group in groups]

rule make_popfile:
    input:
        lambda wildcards: [f"poplists/{population}.txt" for population in groups[wildcards.group]]
    output:
        "popfile/{group}.txt"
    shell:
        """
        cat {input} | sed 's/ /\t/' > {output}
        """

rule sfs:
    input:
        popfile = "popfile/{group}.txt",
        vcf = "groups/{group}.vcf.gz"
    output:
        "SFS/{group}_DSF.obs"
    params:
        basename="SFS/{group}"
    resources:
        runtime = 60 * 12,
        mem = 1024 * 4
    shell:
        """
        
        # /work/FAC/FBM/DBC/amalaspi/popgen/dcruzdav/conda/envs/python2.7/bin/python2.7 vcf2sfs.py -i {input.vcf} -p {input.popfile} -o {params.basename}
        """

rule SFS_1D:
    input:
        vcf = "populations/{population}.vcf.gz",
        popfile = "poplists/{population}.txt"
    output:
        #stats = "stats/{population}_all_sites.txt",
        sfs = "SFS/pop/{population}.sfs"
    params:
        stats = "stats/pop/{population}"
    resources:
        runtime = 60 * 2
    shell:
        """
        python make_sfs_from_vcf.py \
            -i {input.vcf} \
            -p {input.popfile} \
            --prefix_stats {params.stats} \
            --in_popfile {input.popfile} \
            --out_file {output.sfs}
        """

rule SFS_2D:
    input:
        vcf = "groups/{group}.vcf.gz",
        popfile = "popfile/{group}.txt"
    output:
        #stats = "stats/{population}_all_sites.txt",
        bcf = temp("tmp_vcf/groups/{group}.bcf"),
        sfs = "SFS/groups/{group}.sfs"
    params:
        stats = "stats/groups/{group}",
        bcf = "tmp_vcf/groups/{group}"
    resources:
        runtime = 60 * 2
    shell:
        """
        python make_sfs_from_vcf.py \
            -i {input.vcf} \
            -p {input.popfile} \
            --prefix_stats {params.stats} \
            --tmp_vcf {params.bcf} \
            --in_popfile {input.popfile} \
            --out_file {output.sfs}
        """