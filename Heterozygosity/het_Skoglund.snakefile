configfile: "multiple_purposes.yaml"
import os,glob,itertools
include: "parse_resources.smk"

path_vcf = "~/americas/Panels/SGDP/VCF/*"
metadata_simons = "~/americas/Panels/SGDP/Simons_sample_pop_region_country.txt"
# Don't forget to link required scripts
gitpath="/users/dcruzdav/data/Git/Botocudos-scripts/"
myPaths = ["AlleleCounts/count_and_sample.py",
            "Heterozygosity/het_pi_Skoglund.py",
            "Heterozygosity/remove_nondiallelic.r"]
for p in myPaths:
    myCommand = "ln -s " + gitpath + p + " ."
    os.system(myCommand)

chromosomes = [str(chr) for chr in range(1,23)]

def expand_pop_inds(pop):
    all_ind = list(config["Heterozygosity"]["populations"][pop].keys())
    myCombinations = [pair for pair in itertools.combinations(all_ind, 2)]
    prefix = ["{pop}/{ind1}_{ind2}".format(pop = pop, ind1 = ind1, ind2 = ind2)
    for ind1,ind2 in myCombinations]
    return(prefix)

all_pops = list(config["Heterozygosity"]["populations"].keys())
all_prefix = [prefix for pop in all_pops for prefix in expand_pop_inds(pop)]
all_ind = [ind for pop in all_pops for ind in list(config["Heterozygosity"]["populations"][pop].keys())]
African = config["african_het"]

wildcard_constraints:
    Chr =   "|".join(chromosomes)  ,
    african =   African ,
    ind1 = "(" + "|".join([ind for ind in all_ind])+ ")(?!_" + African + ")" 
#=============================================================================#
rule all:
    input:
        expand("{prefix}_{african}.pi.stats.txt", 
        prefix = all_prefix, african = African)

rule fetch_heterozygous:
    input:
        vcf = "{name}.vcf.gz"
    output:
        "sites_het_{name}.txt"
    shell:
        """
        bcftools view -m2 -M2 -Ou  -v snps {input.vcf}  \
                -e  'GT="hom" || (REF ~ "C" & ALT ~ "T") 
                || (REF ~ "G" & ALT ~ "A") || 
                (REF ~ "T" & ALT ~ "C") || (REF ~ "A" & ALT ~ "G")
                || (ALT ~ "C" & ALT ~ "T") ||(ALT ~ "G"& ALT ~ "A")' |\
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        > {output}
        """

#  rule get_variants_all_pops:
#     input:
#         vcf = "{name}.vcf.gz",
#         vcf_tbi = "{name}.vcf.gz.tbi"
#     output:
#         "NonDiallelic/nondiallelic_{name}.txt"
#     shell:
#         """
#         bcftools view -m3 -Ou -v snps {input.vcf} \
#         | bcftools query -f '%CHROM\t%POS\n' > \
#             {output}
#         """

# rule merge_variants:
#     input:
#         expand("NonDiallelic/nondiallelic_{name}.txt")
#     output:
#         "NonDiallelic/all_nondiallelic_{name}.txt"
#     shell:
#         """
#         cat {input} |sort -nk1,2|uniq > {output}
#         """

# rule rm_non_diallelic:
#     input:
#         "NonDiallelic/all_nondiallelic_{name}.txt"
#     output:
#         "Africa/sites_het_{name}.txt"
#     shell:
#         """
#         Rscript remove_nondiallelic.r {input} {output}
#         """

rule keep_chromosomes_in_bam:
    input:
        sites_afr = "sites_het_{african}.txt"
    output:
        refalt = "{african}.refalt"
    shell:
        """
        chromosomes=({chromosomes})
        # Get sites and reference, alternative alleles

        chromInPanel=($(cut -f1 {input.sites_afr} |sort |uniq))
        for chr in ${{chromInPanel[@]}}
        do
            if [[ " ${{chromosomes[*]}} " == *$chr* ]]
            then
                echo $chr is desired
            else
                echo will remove $chr
                sed -i "/$chr/d" {input.sites_afr}
            fi
        done
        cat {input.sites_afr} |\
            sed 's/\t/_/ ; s/A/0/g ; s/C/1/g ; s/G/2/g ; s/T/3/g' |cut -f1-3 >{output.refalt}

        """

rule get_african_sites:
    input:
        sites_afr = "sites_het_{african}.txt"
    output:
        sites_bed = "{african}_{chr}_sites.bed"
    shell:
        """
        cut -f 1,2 {input} |grep -P "^{wildcards.chr}\t"|\
        awk '{{print($1"\t"$2-1"\t"$2)}}' >{output}
        """

def gimme_bam(ind, population):
    myBam = config["Heterozygosity"]["populations"][population][ind]
    return(myBam)

def expand_path(ind, population):
    path = config["Heterozygosity"]["populations"][population][ind]
    full_path = os.path.expanduser(path) 
    bam = glob.glob(full_path)
    return(bam)

rule make_lists:
    input:
        bam1 = lambda wildcards: expand_path(wildcards.ind1, wildcards.population),
        bam2 = lambda wildcards: expand_path(wildcards.ind2, wildcards.population)
    output:
        bams_list = "{population}/{ind1}_{ind2}.txt"
    run:
        with open(output.bams_list, 'w') as output:
            output.write(input.bam1[0] + "\n")
            output.write(input.bam2[0] + "\n")
    
rule do_mpileup:
    input:
        sites_afr = "sites_het_{african}.txt",
        sites_bed = "{african}_{chr}_sites.bed",
        bams_list = "{population}/{ind1}_{ind2}.txt"
    output:
        mpileup = "{population}/{ind1}_{ind2}_{chr}_{african}.mpileup"
    log:
        "{population}/{ind1}_{ind2}_{chr}_{african}.log"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mpileup_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mpileup_time", attempt, 12)
    shell:
        """
        samtools mpileup -r {wildcards.chr} \
            -Bl {input.sites_bed} \
            -b {input.bams_list} \
            -a -o {output.mpileup} \
            -Q20 \
            > {log}
        """

rule count_and_sample_bam:
    input:
        mpileup = "{population}/{ind1}_{ind2}_{chr}_{african}.mpileup",
        refalt = "{african}.refalt"
    output:
        sampled = "{population}/{ind1}_{ind2}_{chr}_{african}.sampled.gz",
        counts = "{population}/{ind1}_{ind2}_{chr}_{african}.counts.gz"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("count_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("count_time", attempt, 12)
    shell:
        """
        python count_and_sample.py \
            --mpileup {input.mpileup} \
            --counts {output.counts} \
            --sampled {output.sampled} \
            --refalt {input.refalt}
        """

rule merge_counts:
    input:
        sampled = expand("{population}/{ind1}_{ind2}_{chr}_{african}.sampled.gz", chr = chromosomes, 
        ind1 = "{ind1}", ind2 = "{ind2}", 
        population = "{population}", african = African)
    output:
        sampled = "{population}/{ind1}_{ind2}_{african}.sampled.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merge_genos_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merge_genos_time", attempt, 2)
    shell:
        """
        for file in {input.sampled}
        do
            gunzip -c $file 
        done > {output.sampled}
        """

rule het_pi_Skoglund:
    input:
        sampled = "{population}/{ind1}_{ind2}_{african}.sampled.txt",
        sites = "{african}.refalt"
    output:
        pi = "{population}/{ind1}_{ind2}_{african}.pi.stats.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("pi_mem", attempt, 8),
        runtime=lambda wildcards, attempt: get_runtime_alloc("pi_time", attempt, 2)
    shell:
        """
        python3.6 het_pi_Skoglund.py \
            -c {input.sampled} \
            -s {input.sites} \
            -o {wildcards.population}/{wildcards.ind1}_{wildcards.ind2}_{wildcards.african}
    """
#-----------------------------------------------------------------------------#
# Rules below are run for data in vcf.gz
# rule get_pops_SGDP:
#     input:
#         metadata = metadata_simons
#     output:
#         "SGDP_pops.txt"
#     shell:
#         """
#         cut -f4 {input} \
#             |sort |uniq |grep -v Population > {output}
#         """

# rule merge_vcf_pop_chr:
#     input:
#         pop_list = "SGDP_pops.txt",
#         sites_bed = "{african}_{chr}_sites.bed"
#     output:
#         vcf_gz = "Merged/{african}_{pop}_{chr}.vcf.gz"
#     params:
#         path_vcf = path_vcf
#     shell:
#         """
#         samples=($(ls {params.path_vcf}/*{wildcards.pop}*vcf.gz))

#         bcftools merge -0R {input.sites_bed} --force-samples \
#          -Oz {wildcards.african} ${{samples[@]}} > {output.vcf}
#         """

# rule index_vcf:
#     input:
#         "{file}.vcf.gz"
#     output:
#         "{file}.vcf.gz.tbi"
#     shell:
#         """
#         bcftools index -t {input}
#         """

# rule count_from_vcf:
#     input:  
#         diallelic = "Africa/sites_het_{african}_{pop}_{chr}.txt",
#         merged_vcf = "Merged/{african}_{pop}_{chr}.vcf.gz"
#     output:
#         counts_pop = "{pop}/{african}.{pop}.{chr}.counts"
#     params:
#         path_vcf = path_vcf
#     shell:
#         """
#         samples=($(ls {params.path_vcf}/*{pop}*vcf.gz))
#         lastField=$( echo ${{#samples[@]}} +1|bc)

#         bcftools view -R {input.diallelic}  -Ou {input.merged_vcf} | \
#          bcftools query -f '[%GT\t]\n' | \
#             cut -f2-$lastField |sed 's|0/1|0\t1|g ; s|1/1|1\t1|g ; s|0/0|0\t0|g' \
#             > {output.counts_pop}

#         """

# rule merge_counts:
#     input:
#         expand("{pop}/{african}.{pop}.{chr}.counts")
#     output:
#         "{pop}/{african}.{pop}.counts.txt "
#     shell:
#         """
#         cat {input} > {output}
#         """

# rule het_pi_Skoglund_called_genos:
#     input:
#         refalt = "{african}.refalt",
#         counts = "{pop}/{african}.{pop}.counts.txt"
#     output:
#     params:
#         seed = 123
#     shell:
#         """
#         ncol=$(head -n1 {input.counts} |wc -c)
#         if [ $ncol -ge 4 ]
#         then
#             python3.5 het_pi_Skoglund.py --counts {input.counts} \
#                 --sites {input.refalt} \
#                 --output {wildcards.pop}/{wildcards.pop} \
#                 --random_seed {params.seed} \
#                 --called_genotypes
#         else
#             touch {output}
#         fi
#         """