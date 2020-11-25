Ne = [100, 1000, 10000, 20000]
replicates = [i for i in range(1, 101)]
chromosomes = [i for i in range(10, 23)]
rule all:
    input:
        expand('Ne_{Ne}/chr{chr}_rep{rep}.txt',
                Ne = Ne,
                chr = chromosomes,
                rep = replicates)
rule msprime:
    input:
        recomb_map = '/Users/dcruz/Projects/Botocudos/Files/HapMap_genetic_maps/HapmapII_GRCh37_RecombinationHotspots/genetic_map_GRCh37_chr{chr}.txt',
        positions = 'chr{chr}.positions'
    output:
        'Ne_{Ne}/chr{chr}_rep{rep}.txt'
    params:
        sample_size = 2
    # conda:
        # 'msprime.yaml'
    shell:
        ''' 
        set +e
        # conda activate msprime
        python simulate_chromosomes.py \
            --recombination-map {input.recomb_map} \
            --positions {input.positions} \
            --Ne {wildcards.Ne} \
            --sample-size {params.sample_size} \
            --chromosome {wildcards.chr} \
            --num-replicates 1 \
            | tail -n1 > {output}
        
        exitcode=$? 
        if [ $exitcode -eq 1 ] ; then exit 1; else exit 0 ; fi
        '''