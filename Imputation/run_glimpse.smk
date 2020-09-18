configfile: "multiple_purposes.yaml"
import subprocess

ref_genome = config['ref_genome']
include: "parse_resources.smk"
bamlists = list(config["GLIMPSE"]["bamlists"].keys())

chromosomes = [ str(chr) for chr in range(1, 23) ]

rule all:
    input:
        vcf = expand("vcf/{bamlist}_chr{chr}.vcf.gz",  bamlist = bamlists,
        chr = chromosomes),
        phased_haplotypes = expand( "GLIMPSE_phased/{bamlist}.chr{chr}.phased.bcf",
                                    bamlist = bamlists, chr = chromosomes )


def expand_path( bamlist ):
    paths = list(config["GLIMPSE"]["bamlists"][ bamlist ]["paths"].values()) 
    #full_paths = [os.path.expanduser(p) for p in paths]
    bams = []
    for p in paths:
        myCommand = "shopt -s extglob ; ls " + p 
        files = os.popen(myCommand).read().split('\n')
        files.remove('')
        [bams.append(individual_path) for individual_path in files]

    return(bams)

rule make_bamlist:
    input:
        lambda wildcards: expand_path( wildcards.bamlist )
    output:
        "bamlists/{bamlist}.txt"
    log:
        "logs/{bamlist}_make_bamlist.log"
    run:
        with open(output[0], 'w') as file:
            for line in input:
                file.write(line+"\n")

# extract variable positions
rule extract_positions:
    input:
        bcf = "rm_singleton/{chr}.no_singletons.bcf"
    output:
        sites = "sites/{chr}.sites.vcf.gz",
        tsv = "sites/{chr}.sites.tsv.gz"
    log:
        'logs/extract_pos_{chr}.log'
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps \
        {input.bcf} \
        -Oz -o {output.sites}

        bcftools index -f {output.sites}
        bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' {output.sites} | \
            bgzip -c > {output.tsv}
        tabix -s1 -b2 -e2 {output.tsv}
        """

rule bcftools_mpileup:
    input:
        ref_genome = config['ref_genome'],
        sites = "sites/{chr}.sites.vcf.gz",
        bamlist = "bamlists/{bamlist}.txt",
        tsv = "sites/{chr}.sites.tsv.gz"
    output:
        vcf = "vcf/{bamlist}_chr{chr}.vcf.gz"
    threads:
        2
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mpileup_mem", attempt, 8),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mpileup_time", attempt, 24)
    log:
        'logs/bcftools_mpileup_{bamlist}_{chr}.txt'
    shell:
        """
        bcftools mpileup -f {input.ref_genome} \
            -I -E -a 'FORMAT/DP' \
            -T {input.sites} -r {wildcards.chr} \
            -b {input.bamlist} -Ob --threads {threads} | \
            bcftools call -Aim -C alleles -T {input.tsv} -Oz -o {output.vcf}
        """

def get_num_chunks( chr ):
    myCommand = 'wc -l chunks/chunks.chr' + chr + '.txt'
    proc = subprocess.Popen(myCommand, stdout = subprocess.PIPE, shell = True)
    (out, err) = proc.communicate()
    n_chunks = int(out.decode().split()[0])

    return(n_chunks)

rule index_bcf:
    input:
        bcf = "{file}.bcf"
    output:
        csi = "{file}.bcf.csi"
    log:
        'logs/index_bcf_{file}.txt'
    shell:
        """
        bcftools index -f {input.bcf}
        """

rule index_vcf:
    input:
        vcf = "{file}.vcf.gz"
    output:
        csi = "{file}.vcf.gz.csi"
    log:
        'logs/index_vcf_{file}.vcf.gz'
    shell:
        """
        bcftools index -f {input.vcf}
        """

# # 5. Impute and phase a whole chromosome
rule impute_phase:
    input:
        chunks = "chunks/chunks.chr{chr}.txt",
        vcf_sample = "vcf/{bamlist}_chr{chr}.vcf.gz",
        csi_sample = "vcf/{bamlist}_chr{chr}.vcf.gz.csi",
        vcf_ref = "rm_singleton/{chr}.no_singletons.bcf",
        map = "maps/genetic_maps.b37/chr{chr}.b37.gmap.gz"
    output:
        imputed_bcf = temp( "GLIMPSE_imputed/{bamlist}_chr{chr}.{n}.bcf" ),
        # imputed_csi = temp( "GLIMPSE_imputed/{bamlist}_chr{chr}.{n}.bcf.csi" )
    log:
        'logs/impute_phase_{bamlist}_{chr}_{n}.txt'
    threads:
        8
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("phase_mem", attempt, 16),
        runtime=lambda wildcards, attempt: get_runtime_alloc("phase_time", attempt, 2)
    shell:
        """

        LINE=$( head -n {wildcards.n} {input.chunks} |tail -n1 )
        IRG=$(echo $LINE | cut -d" " -f3)
        ORG=$(echo $LINE | cut -d" " -f4)

        GLIMPSE_phase --input {input.vcf_sample} \
            --reference {input.vcf_ref} \
            --map {input.map} \
            --input-region ${{IRG}} \
            --output-region ${{ORG}} \
            --output {output.imputed_bcf} \
            --thread {threads}

        """


#  # 6. Ligate multiple chunks together
rule ligate_chunks:
    input:
        phased = lambda wildcards: [ "GLIMPSE_imputed/{bamlist}_chr{chr}.{n}.bcf".format(
                                            bamlist = wildcards.bamlist, 
                                            chr = wildcards.chr,
                                            n = n )
                                                for n in range(1, get_num_chunks( wildcards.chr ) ) 
        ],
        index = lambda wildcards: [ "GLIMPSE_imputed/{bamlist}_chr{chr}.{n}.bcf.csi".format(
                                            bamlist = wildcards.bamlist, 
                                            chr = wildcards.chr,
                                            n = n)
                                                for n in range(1, get_num_chunks( wildcards.chr ) ) 
        ]
    output:
        ligated_bcf = "GLIMPSE_ligated/{bamlist}.chr{chr}.merged.bcf",
        list_files = temp("GLIMPSE_ligated/list_{bamlist}_chr{chr}.txt")
    log:
        'logs/ligate_chunks_{bamlist}_{chr}.txt'
    shell:
        """
        for file in {input.phased} ; do echo $file ; done > {output.list_files}
        GLIMPSE_ligate --input {output.list_files} --output {output.ligated_bcf}

        """

#  # 7. Sample haplotypes
rule sample_haplotypes:
    input:
        ligated_bcf = "GLIMPSE_ligated/{bamlist}.chr{chr}.merged.bcf",
        ligated_csi = "GLIMPSE_ligated/{bamlist}.chr{chr}.merged.bcf.csi"
    output:
        phased_bcf = "GLIMPSE_phased/{bamlist}.chr{chr}.phased.bcf"
    log:
        'logs/sample_haplotype_{bamlist}_{chr}.txt'
    shell:
        """
        GLIMPSE_sample --input {input.ligated_bcf} \
        --solve \
        --output {output.phased_bcf}
        """
