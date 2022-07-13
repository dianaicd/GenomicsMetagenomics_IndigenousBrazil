
configfile: "config_merge_vcfs.yaml"
localrules: mv_results
replicates = config["fastsimcoal2"]["replicates"]
models = config["fastsimcoal2"]["models"]

wildcard_constraints:
    model = "|".join([key for key in models.keys()])
    
rule all:
    input:
        [
            f"sim_{group}/{group}_{model}_rep{rep}/{group}_{model}_rep{rep}.bestlhoods"
            for model in models.keys()
            for group in models[model]
            for rep in range(1, replicates + 1)
        ]

rule fastsimcoal2_2Pop:
    input:
        template = "{model}.tpl",
        est = "{model}.est",
        sfs = "SFS/groups/{group}.sfs",
    output:
        template = temp("{group}_{model}_rep{rep}.tpl"),
        est = temp("{group}_{model}_rep{rep}.est"),
        sfs = temp("{group}_{model}_rep{rep}_jointDAFpop1_0.obs"),
        likelihood = "{group}_{model}_rep{rep}/{group}_{model}_rep{rep}.bestlhoods"
    params:
        n_sim = 500000,
        loops = 50
    resources:
        runtime = 60 * 2
    shell:
        """
        cp {input.template} {output.template}
        cp {input.est} {output.est}

        echo "1 observation" > {output.sfs}
        sed 's/,/\t/g' {input.sfs} >> {output.sfs}

        ./fsc2709 \
            -t {output.template} -e {output.est} \
            -n {params.n_sim} -d -M \
            -L {params.loops} -q
        """

rule fastsimcoal2_multipop:
    input:
    output:
    params:
    resources:
    shell:
        """
        # colnames are d0_i and rownames are d1_i
        # *pop1_0 has colnames as d0_i and rownames as d1_i ; it should be pop0_pop1.sfs
        # *pop2_0 has colnames as d0_i and rownames as d2_i ; it should be pop0_pop2.sfs
        # *pop2_1 has colnames as d1_i and rownames as d2_i ; it should be pop1_pop2.sfs

        """

rule mv_results:
    input:
        likelihood = "{group}_{model}_rep{rep}/{group}_{model}_rep{rep}.bestlhoods"
    output:
        "sim_{group}/{group}_{model}_rep{rep}/{group}_{model}_rep{rep}.bestlhoods"
    params:
        input = "{group}_{model}_rep{rep}",
        output = "sim_{group}/"
    shell:
        """
        mkdir -p sim_{wildcards.group}
        mv {params.input} {params.output}
        """