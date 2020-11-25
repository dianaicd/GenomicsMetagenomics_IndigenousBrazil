configfile: 'multiple_purposes.yaml'

chromosomes = [str(i) for i in range(1, 23)]

input_files = list(config['smcpp'].keys())
pop_inds = {file:config['smcpp'][file]['populations'] for file in input_files}
header_path = config['smcpp_header'] if 'smcpp_header' in config.keys() else 'MN0008.header'
use_docker = config['smcpp_docker'] if 'smcpp_docker' in config.keys() else True
mount_dir = config['smcpp_mount_dir'] if 'smcpp_mount_dir' in config.keys() else "$PWD"

chr_length = {}

with open(header_path, 'r') as header:
    for line in header.readlines():
        if "@SQ\tSN:" in line:
            chr = line.split()[1].replace("SN:", "")
            l = line.split()[2].replace("LN:", "").replace("\n", "")
            chr_length[chr] = l
# chr_length = {'1':249250621}



def smc_command(use_docker, mount_dir):
    if use_docker:
        command = f'docker run --rm -v {mount_dir}:/mnt terhorst/smcpp:latest '
    else:
        command = 'smc++ '
    return command 


rule all:
    input:
        Ne = ['plot/{file}.{population}.csv'.format(file=file, population=population) 
                for file in input_files for population in pop_inds[file].keys()]

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
        tabix {input.vcf}
        '''

rule vcf2smc:
    input:
        docker_img = 'Dockerfile',
        vcf = lambda wildcards: config['smcpp'][wildcards.file]['prefix'] + f'.chr{wildcards.chr}.vcf.gz',
        vcf_index = lambda wildcards: config['smcpp'][wildcards.file]['prefix'] + f'.chr{wildcards.chr}.vcf.gz.tbi'
    output:
        smc = 'smc/{file}.{population}.chr{chr}.smc.gz'
    threads: 2
    params:
        command = smc_command(use_docker, mount_dir = mount_dir),
        individuals = lambda wildcards: pop_inds[wildcards.file][wildcards.population],
        chr_length = lambda wildcards: chr_length[wildcards.chr]
    shell:
        '''
        inds=$( echo {params.individuals} | sed -e 's/[0-9a-zA-Z_-]*/&_&/g ; s/ /,/g ; s/,$//' )
        {params.command} \
            vcf2smc {input.vcf} {output.smc} {wildcards.chr} \
            {wildcards.population}:$inds \
            --length {params.chr_length} --cores {threads}
        '''

rule estimate:
    input:
        smc = expand('smc/{file}.{population}.chr{chr}.smc.gz', chr = chromosomes,
                    file = '{file}', population = '{population}')
    output:
        model = 'analysis/{file}.{population}/model.final.json'
    params:
        command = smc_command(use_docker, mount_dir = mount_dir),
        out_prefix = 'analysis/{file}.{population}',
        mutation_rate = 1.25e-8
    threads: 6
    shell:
        '''
        {params.command} \
            estimate -o {params.out_prefix} \
            {params.mutation_rate} \
            {input.smc} \
            --spline cubic \
            --cores {threads}
        '''

rule plot:
    input:
        model = 'analysis/{file}.{population}/model.final.json'
    output:
        csv = 'plot/{file}.{population}.csv',
        pdf = 'plot/{file}.{population}.pdf'
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