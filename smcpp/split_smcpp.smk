configfile: 'multiple_purposes.yaml'

chromosomes = [str(i) for i in range(1, 23)]


pop_inds = config['smcpp']['populations']
use_docker = config['smcpp_docker'] if 'smcpp_docker' in config.keys() else True
mount_dir = config['smcpp_mount_dir'] if 'smcpp_mount_dir' in config.keys() else "$PWD"
populations = list(pop_inds.keys())
combinations = [f"{populations[i]}.{populations[j]}" for i in range(0, len(populations)) for j in range(i+1, len(populations))]
wildcard_constraints:
    pop1 = "|".join([pop for pop in populations]),
    population = "|".join([pop for pop in populations])

def smc_command(use_docker, mount_dir):
    if use_docker:
        command = f'docker run --rm -v {mount_dir}:/mnt terhorst/smcpp:latest '
    else:
        command = 'smc++ '
    return command 

rule all:
    input:
        Ne = ['plot/{population}.csv'.format(population=population) 
                for population in pop_inds.keys()],
        splits = [f"split/{comb}/model.final.json" for comb in combinations]

rule check_programs:
    input: 
        docker_img = 'Dockerfile',
        src = 'src'

rule index_vcf:
    input:
        vcf = '{file}.vcf.gz'
    output:
        tbi = '{file}.vcf.gz.tbi'
    shell:
        '''
        tabix -f {input.vcf}
        '''

rule vcf2smc:
    input:
        #docker_img = 'Dockerfile',
        vcf = config['smcpp']['prefix'] ,
        vcf_index = config['smcpp']['prefix'] + '.tbi',
        #mask = config["smcpp"]["mask"]
        # mask = "centromeres/chr{chr}.bed.gz"
        mask = "input_smcpp/merged_masks/{population}.chr{chr}.bed.gz"
    output:
        smc = 'smc/{population}.chr{chr}.smc.gz'
    threads: 2
    params:
        command = smc_command(use_docker, mount_dir = mount_dir),
        individuals = lambda wildcards: pop_inds[wildcards.population],
    shell:
        '''
        #inds=$( echo {params.individuals} | sed -e 's/[0-9a-zA-Z_-]*/&_&/g ; s/ /,/g ; s/,$//' )
        inds=$( echo {params.individuals} | sed -e 's/ /,/g ; s/,$//' )
        {params.command} \
            vcf2smc --cores {threads} --mask {input.mask} {input.vcf} {output.smc} {wildcards.chr} \
            {wildcards.population}:$inds 
        '''

rule estimate:
    input:
        smc = expand('smc/{population}.chr{chr}.smc.gz', chr = chromosomes,
                    population = '{population}')
    output:
        model = 'analysis/{population}/model.final.json'
    params:
        command = smc_command(use_docker, mount_dir = mount_dir),
        out_prefix = 'analysis/{population}',
        name = lambda wildcards: f"smc/{wildcards.population}.chr*.smc.gz",
        mutation_rate = 1.25e-8
    threads: 6
    shell:
        '''
        {params.command} \
            estimate -o {params.out_prefix} \
            --timepoints 33 100000 -c 50000 -rp .1 --knots 60 \
            --cores {threads} \
            --spline cubic \
            {params.mutation_rate} \
            {params.name} 
        '''

rule plot:
    input:
        model = 'analysis/{population}/model.final.json'
    output:
        csv = 'plot/{population}.csv',
        pdf = 'plot/{population}.pdf'
    params:
        generation = 25,
        command = smc_command(use_docker, mount_dir = mount_dir),
        x = "0 1e6",
        y = "0 1e5"
    shell:
        '''
        {params.command} \
            plot {output.pdf} {input.model} --csv -g {params.generation} -x {params.x} -y {params.y}
        '''





# Next, create datasets containing the joint frequency spectrum for both populations:

# $ smc++ vcf2smc my.vcf.gz data/pop12.smc.gz <contig> pop1:ind1_1,ind1_2 pop2:ind2_1,ind2_2
# $ smc++ vcf2smc my.vcf.gz data/pop21.smc.gz <contig> pop2:ind2_1,ind2_2 pop1:ind1_1,ind1_2

rule joint_vcf2smc:
    input:
        #docker_img = 'Dockerfile',
        vcf = config['smcpp']['prefix'] ,
        vcf_index = config['smcpp']['prefix'] + '.tbi',
        # bed = "centromeres/chr{chr}.bed.gz",
        bed = "input_smcpp/merged_masks/{pop1}.{pop2}.chr{chr}.bed.gz"
    output:
        smc = 'smc/{pop1}_{pop2}.chr{chr}.smc.gz'
    threads: 2
    params:
        command = smc_command(use_docker, mount_dir = mount_dir),
        inds_pop1 = lambda wildcards: pop_inds[wildcards.pop1],
        inds_pop2 = lambda wildcards: pop_inds[wildcards.pop2],
    shell:
        '''
        pop1_inds=$( echo {params.inds_pop1} | sed -e 's/ /,/g ; s/,$//' )
        pop2_inds=$( echo {params.inds_pop2} | sed -e 's/ /,/g ; s/,$//' )
        {params.command} \
            vcf2smc {input.vcf} \
            --cores {threads} \
            --mask {input.bed} \
            {output.smc} {wildcards.chr} \
             {wildcards.pop1}:${{pop1_inds}} {wildcards.pop2}:${{pop2_inds}}
        '''

rule split:
    input:
        model_pop1 = 'analysis/{pop1}/model.final.json',
        model_pop2 = 'analysis/{pop2}/model.final.json',
        smc_files = lambda wildcards: 
            expand(["smc/{pop1}_{pop2}.chr{chr}.smc.gz",
            "smc/{pop2}_{pop1}.chr{chr}.smc.gz",
            'smc/{pop1}.chr{chr}.smc.gz',
            'smc/{pop2}.chr{chr}.smc.gz'],
             chr = chromosomes, pop1 = wildcards.pop1, pop2 = wildcards.pop2
            )
    output:
        out = "split/{pop1}.{pop2}/model.final.json"
    params:
        command = smc_command(use_docker, mount_dir = mount_dir),
        out = "split/{pop1}.{pop2}/"
    shell:
        """
        files=$(ls smc/{{{wildcards.pop1}_{wildcards.pop2}.chr,{wildcards.pop1}.chr,{wildcards.pop2}.chr}}*.smc.gz)
        {params.command} \
        split -o {params.out} {input.model_pop1} {input.model_pop2} $files
        """

# Finally, run split to refine the marginal estimates into an estimate of the joint demography:

# $ smc++ split -o split/ pop1/model.final.json pop2/model.final.json data/*.smc.gz
# $ smc++ plot joint.pdf split/model.final.json

rule merge_masks_inds:
    input:
        lambda wildcards: [f"input_smcpp/called_masks/{ind}/{ind}_chr{wildcards.chr}.bed" for ind in pop_inds[wildcards.population]]
    output:
        tmp_bed = "input_smcpp/merged_masks/{population}.chr{chr}.tmp.bed",
        bed = "input_smcpp/merged_masks/{population}.chr{chr}.bed.gz"
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output.tmp_bed}
        bedtools merge -i {output.tmp_bed} | bgzip > {output.bed}
        tabix -f {output.bed}
        """

rule merge_masks_pops:
    input:
        pop1="input_smcpp/merged_masks/{pop1}.chr{chr}.bed.gz",
        pop2="input_smcpp/merged_masks/{pop2}.chr{chr}.bed.gz"
    output:
        tmp_bed = "input_smcpp/merged_masks/{pop1}.{pop2}.chr{chr}.tmp.bed",
        bed = "input_smcpp/merged_masks/{pop1}.{pop2}.chr{chr}.bed.gz"
    shell:
        """
        gunzip -c {input} | sort -k1,1 -k2,2n > {output.tmp_bed}
        bedtools merge -i {output.tmp_bed} | bgzip > {output.bed}
        tabix -f {output.bed}
        """