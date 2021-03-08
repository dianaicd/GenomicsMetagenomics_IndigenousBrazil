configfile: "multiple_purposes.yaml"

rule all:
    input:
        expand(
            "ibs/{bamlist}.chr{chr}.txt",
            bamlist = [bamlist for bamlist in config["IBS"]["bamlists"].keys()],
            chr = [str(chr) for chr in range(1, 23)]
        )

rule genotype_likelihoods:
    input:
        bamlist = "{bamlist}.txt"
    output:
        genotypes = "{bamlist}.chr{chr}.glf.gz"
    params:
        basename = "{bamlist}.chr{chr}"
    resources:
        runtime = 6*60,
        memory = 8*1024
    shell:
        """
        angsd -b {input.bamlist} \
            -minMapQ 30 -minQ 20 -GL 2 \
            -doGlf 1 \
            -r {wildcards.chr}: \
            -out {params.basename}
        """

rule ibs:
    input:
        bamlist = "{bamlist}.txt",
        genotypes = "{bamlist}.chr{chr}.glf.gz"
    output:
        results = "results_IBS/ibs.model0.{bamlist}.chr{chr}.results.ibspair"
    resources:
        runtime = 23*60,
        memory = 24*1024
    params:
        basename = "results_IBS/ibs.model0.{bamlist}.chr{chr}.results"
    shell:
        """
        nInd=$(wc -l {input.bamlist} | cut -f1 -d " ")

        /software/UHTS/Analysis/ANGSD/0.931/misc/ibs -glf {input.genotypes} \
            -model 0 \
            -nInd $nInd -allpairs 1 \
            -outFileName {params.basename}
        """

rule examine_results:
    input:
        results = "results_IBS/ibs.model0.{bamlist}.chr{chr}.results.ibspair"
    output:
        examination = "ibs/{bamlist}.chr{chr}.txt"
    shell:
        """
        Rscript \
            -e "source('./read_IBS.R')" \
            -e "res = do_derived_stats(read_ibspair_model0('{input.results}'))" \
            -e "print(res[6,c('ind1', 'ind2', 'nSites', 'Kin', 'R0', 'R1') ])" \
                > {output.examination}
        """
