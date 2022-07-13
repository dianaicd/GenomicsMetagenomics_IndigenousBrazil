
configfile: "config_merge_vcfs.yaml"
localrules: mv_results


    
rule all:
    input:
        [
            f"sim_{model}/{model}_rep{rep}/{model}_rep{rep}.bestlhoods"
            for model in ["5PopMayaAym_mig"]
            for rep in range(1, 100 + 1)
        ]


rule fastsimcoal2_multipop:
    input:
        template = "{model}.tpl",
        est = "{model}.est",
        sfs1 = "{model}_jointDAFpop1_0.obs",
        sfs2 = "{model}_jointDAFpop2_0.obs",
        sfs3 = "{model}_jointDAFpop2_1.obs",
        sfs4 = "{model}_jointDAFpop3_0.obs",
        sfs5 = "{model}_jointDAFpop3_1.obs",
        sfs6 = "{model}_jointDAFpop3_2.obs",
        sfs7 = "{model}_jointDAFpop4_0.obs",
        sfs8 = "{model}_jointDAFpop4_1.obs",
        sfs9 = "{model}_jointDAFpop4_2.obs",
        sfs10 = "{model}_jointDAFpop4_3.obs"
    output:
        template = temp("{model}_rep{rep}.tpl"),
        est = temp("{model}_rep{rep}.est"),
        sfs1 = temp("{model}_rep{rep}_jointDAFpop1_0.obs"),
        sfs2 = temp("{model}_rep{rep}_jointDAFpop2_0.obs"),
        sfs3 = temp("{model}_rep{rep}_jointDAFpop2_1.obs"),
        sfs4 = temp("{model}_rep{rep}_jointDAFpop3_0.obs"),
        sfs5 = temp("{model}_rep{rep}_jointDAFpop3_1.obs"),
        sfs6 = temp("{model}_rep{rep}_jointDAFpop3_2.obs"),
        sfs7 = temp("{model}_rep{rep}_jointDAFpop4_0.obs"),
        sfs8 = temp("{model}_rep{rep}_jointDAFpop4_1.obs"),
        sfs9 = temp("{model}_rep{rep}_jointDAFpop4_2.obs"),
        sfs10 = temp("{model}_rep{rep}_jointDAFpop4_3.obs"),
        likelihood = "{model}_rep{rep}/{model}_rep{rep}.bestlhoods"
    params:
        n_sim = 500000,
        loops = 50
    resources:
        runtime = 2 * 60
    threads: 8
    shell:
        """
        # colnames are d0_i and rownames are d1_i
        # *pop1_0 has colnames as d0_i and rownames as d1_i ; it should be pop0_pop1.sfs
        # *pop2_0 has colnames as d0_i and rownames as d2_i ; it should be pop0_pop2.sfs
        # *pop2_1 has colnames as d1_i and rownames as d2_i ; it should be pop1_pop2.sfs
        cp {input.template} {output.template}
        cp {input.est} {output.est}
        cp {input.sfs1} {output.sfs1}
        cp {input.sfs2} {output.sfs2}
        cp {input.sfs3} {output.sfs3}
        cp {input.sfs4} {output.sfs4}
        cp {input.sfs5} {output.sfs5}
        cp {input.sfs6} {output.sfs6}
        cp {input.sfs7} {output.sfs7}
        cp {input.sfs8} {output.sfs8}
        cp {input.sfs9} {output.sfs9}
        cp {input.sfs10} {output.sfs10}

        ./fsc2709 \
            -t {output.template} -e {output.est} \
            -n {params.n_sim} -d -M \
            -L {params.loops} -q -c {threads}
        """

rule mv_results:
    input:
        likelihood = "{model}_rep{rep}/{model}_rep{rep}.bestlhoods"
    output:
        "sim_{model}/{model}_rep{rep}/{model}_rep{rep}.bestlhoods"
    params:
        input = "{model}_rep{rep}",
        output = "sim_{model}/"
    shell:
        """
        mkdir -p sim_{wildcards.model}
        mv {params.input} {params.output}
        """