# configfile: "multiple_purposes.yaml"
configfile: "genotype_calling.yaml"
include: "parse_resources.smk"
# Call genotypes on a medium/high coverage sample
# output all possible sites
#-----------------------------------------------------------------------------#
# Variables and functions to begin with
def param_is_defined(name, default_par = False):
    myParameter = config[name] if name in config.keys() else default_par
    return(myParameter)

def expand_path(myBamfile):
    mySample = list(config["geno_calls"]["bamlists"][myBamfile]["paths"].keys())[0]
    path = config["geno_calls"]["bamlists"][myBamfile]["paths"][mySample]
    full_path = os.path.expanduser(path)
    return(full_path)

def get_range(chr = False):
    if(chr):
        ranges = ["_".join([str(chr), 
                    str(lower_chr[str(chr)][i]), 
                    str(upper_chr[str(chr)][i])]) for i in range(0, len(lower_chr[str(chr)]))]
    else:
        ranges = ["_".join([str(chr), 
                    str(lower_chr[str(chr)][i]), 
                    str(upper_chr[str(chr)][i])]) for chr in chromosomes for i in range(0, len(lower_chr[str(chr)]))]
    return(ranges)
#-----------------------------------------------------------------------------#

# 
bamlists = list(config["geno_calls"]["bamlists"].keys())

dict_bamlists = config["geno_calls"]["bamlists"]

myDict = config["geno_calls"]["bamlists"]

ref_genome = param_is_defined(name = "ref_genome", 
                            default_par = "/scratch/axiom/FAC/FBM/DBC/amalaspi/popgen/reference_human/hs.build37.1/hs.build37.1.fa"
)

blockSize = param_is_defined(name = "blockSize", default_par = 5e6)

with open(ref_genome + ".fai", 'r') as index:
    chr_size = {}
    for line in index.readlines():
        chr,size = line.split()[0:2]
        chr_size[chr] = size
    chromosomes = list(chr_size.keys())
    chromosomes = [str(x) for x in range(1,23)] 
    lower_chr = {}
    upper_chr = {}
    for chr in chromosomes:
        chr = str(chr)
        lower_chr[chr] = [i for i in range(1, int(chr_size[chr]), int(float(blockSize)))]
        upper_chr[chr] = [i - 1 for i in lower_chr[chr][1:]]
        upper_chr[chr].append(chr_size[chr])

chromosomes = [str(i) for i in range(1, 23)]

wildcard_constraints:
    chr = "(" + "|".join([i for i in chromosomes]) + ")"
#-----------------------------------------------------------------------------#
# Get chromosome size and break it in blocks
# include: "break_blocks.smk"
# include: "make_bamlist.smk"
#-----------------------------------------------------------------------------#

rule all:
    input:
        expand(
            "Raw/{bamfile}_chr{chr}.bcf",
            bamfile = bamlists,
            chr = chromosomes
            )
        # bed_ind = expand("Filtered/{bamlist}_chr{chr}_depth_filter_{rmTrans}.bed",
                        # bamlist = bamlists, chr = chromosomes, rmTrans = ["all", "rmTrans"])

rule print_positions:
    output:
        positions = "Raw/{bamlist}_positions.txt"
    run:
        with open(output.positions, "w") as out:
            [out.write(line.replace("_", "\t")+"\n") for chr in chromosomes for line in get_range(chr)]

rule call_genos:
    input:
        bamlist = "{bamlist}/{bamlist}.txt"#lambda wildcards: expand_path(wildcards.bamfile),
    output:
        raw_genos = temp("Raw/{bamlist}_{chr}_{start}_{end}.bcf")
    log:
        "logs/{bamlist}_{chr}_{start}_{end}.log"
    params:
        minMapQ=30,
        minBaseQ=20,
        ref =param_is_defined(name = "ref_genome"),
        threads = 4,
        moreno2019 = param_is_defined("Genos_Moreno2019", "Yes")
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("call_genos_mem", attempt, 4),
        runtime= 5#lambda wildcards, attempt: get_runtime_alloc("call_genos_time", attempt, 2)
    shell:
        """
            if [ {params.moreno2019} == "Yes" ]
            then
               samtools mpileup -q {params.minMapQ} \
                    -t DP -C50 -uf {params.ref} \
                    --region {wildcards.chr}:{wildcards.start}-{wildcards.end} \
                -b {input.bamlist} | bcftools call -f GQ -c \
                    --threads {params.threads} \
                    -Ob -o {output.raw_genos} - \
                2>{log}
            else
                bcftools mpileup -C 50 -q {params.minMapQ} -Q {params.minBaseQ} -a FMT/DP,SP \
                        -f {params.ref} --threads {params.threads} \
                        --regions {wildcards.chr}:{wildcards.start}-{wildcards.end} \
                        -Ou  {input.bamlist} | \
                    bcftools annotate -c RPB | \
                    bcftools call --threads {params.threads} -c -V indels \
                        -Ob -o {output.raw_genos} - \
                2>{log}
            fi
        """

rule concat_genos_chr:
    input:
        positions = "Raw/{bamfile}_positions.txt",
        bcf = lambda wildcards: expand("Raw/{bamfile}_{range}.bcf", 
                    bamfile = "{bamfile}",
                    range =  get_range(wildcards.chr))
    wildcard_constraints:
        chr = "|".join([str(x) for x in chromosomes])
    output:
        genos_list = "Raw/{bamfile}_{chr}.list",
        genos_all_depths = "Raw/{bamfile}_chr{chr}.bcf"
    params:  
        threads = 4
    log:
    run:
        def write_line(line, myChr):
            chr,start,end = line.replace("\n", "").split()
            myNewLine = "Raw/{bamfile}_{chr}_{start}_{end}.bcf\n".format(bamfile = wildcards.bamfile, 
                                                                                                    chr = chr, start = start, 
                                                                                                    end = end)
            if chr == myChr:
                genos_list.write(myNewLine)
            #return(myNewLine)

        with open(input.positions, 'r') as positions, open(output.genos_list, 'w') as genos_list:
            [write_line(line, wildcards.chr) for line in positions.readlines()]

        myCommand = "bcftools concat -f {bcf_list} -Ob -o {out} --threads {thr}".format(bcf_list = output.genos_list, 
        out = output.genos_all_depths, thr = params.threads)
        os.system(myCommand)

rule filter_genos:
    input:
        raw_genos = "Raw/{bamfile}_chr{chr}.bcf"
    output:
        filtered_genos = temp("Filtered/{bamfile}_chr{chr}_allDepths_all.bcf")
    params:
        threads = 20,
        genoqual = 30,
        allelic_imbalance = 0.2,
        read_position_bias = 0,
        variant_distance_bias = 0,
        moreno2019 = param_is_defined("Genos_Moreno2019", "Yes")

    log:
        "logs/{bamfile}_chr{chr}.log"
    shell:
        """
        if [ {params.moreno2019} == "Yes" ]
        then

            bcftools filter --threads {params.threads} \
            --SnpGap 5 --IndelGap 5 \
             --exclude "QUAL<{params.genoqual} | 
             PV4[0]<0.0001 | PV4[1]=0 | PV4[2]=0 | PV4[3]<0.0001 | 
             RPB=0 | 
             (GT='het' && 
                        ( (DP4[0]+DP4[1])/SUM(DP4) <= {params.allelic_imbalance} ||
                        (DP4[2]+DP4[3])/SUM(DP4) <= {params.allelic_imbalance}  )
                    )"\
            {input.raw_genos} | grep -v INDEL |\
            bcftools view -Ob -o {output.filtered_genos}

        else

        bcftools filter --threads {params.threads} --exclude \
            "(SP > 40) 
                || (GT='het' && 
                        ( (DP4[0]+DP4[1])/SUM(DP4) <= {params.allelic_imbalance} ||
                        (DP4[2]+DP4[3])/SUM(DP4) <= {params.allelic_imbalance}  )
                    )
                || (REF='C' & ALT='T') || (REF=='T' & ALT == 'C')
                || (REF='G' & ALT='A') || (REF=='A' & ALT == 'G')
                || N_ALT > 1 || QUAL <= {params.genoqual} 
                || VDB == {params.variant_distance_bias} 
                || RPB == {params.read_position_bias}" \
            -o {output.filtered_genos} {input.raw_genos} 2>> {log}
        fi    
        """

rule filter_transitions:
    input:
        genos = "Filtered/{bamfile}_chr{chr}_allDepths_all.bcf"
    output:
        genos = temp("Filtered/{bamfile}_chr{chr}_allDepths_rmTrans.bcf")
    params:
        threads = 4
    shell:
        """
        bcftools filter --threads {params.threads} \
            --exclude "(REF='C' & ALT='T') || (REF=='T' & ALT == 'C')|| (REF='G' & ALT='A') || (REF=='A' & ALT == 'G')" \
            -Ob -o {output.genos} {input.genos} 
        """


rule filter_depth_chr:
    input:
        genos_all_depths = "Filtered/{bamfile}_chr{chr}_allDepths_{rmTrans}.bcf"
    output:
        genos_depth_filter = temp("Filtered/{bamfile}_chr{chr}_depth_filter_{rmTrans}.bcf"),
        # stats = "Filtered/{bamfile}_{chr}_stats_{rmTrans}.txt",
        # new_name = "Filtered/{bamfile}_{chr}_new_sample_name_{rmTrans}.txt"
    params:
    # Minimum and maximum depth relative to average depth of coverage 
    # on filtered genotypes in this sample
        min_depth = "1/3",
        max_depth = "2",
        # sample = lambda wildcards: list(config["geno_calls"]["bamlists"][wildcards.bamfile]["paths"].keys())[0]
    log:
        "logs/filter_depth_{bamfile}_{chr}_{rmTrans}.log"
    threads: 20
    shell:
        """
        set +e
        # sample=$(bcftools query -f '[%SAMPLE]\n' {input.genos_all_depths} \
        #             |head -n1) >{log} 2>&1

        avgdp=$(bcftools stats -s $sample {input.genos_all_depths} \
            |grep -P "^PSC" |cut -f 10); echo "average depth: $avgdp " >>{log} 2>&1

        mindp=$avgdp*{params.min_depth} >>{log} 2>&1
        maxdp=$avgdp*{params.max_depth} >>{log} 2>&1

        # echo "$sample {params.sample}" > {output.new_name}  2>>{log}
        bcftools filter \
            --threads {threads} \
            -Ob \
            --exclude "(SUM(DP4)< $mindp| SUM(DP4) > $maxdp )" \
            -o {output.genos_depth_filter} \
            {input.genos_all_depths} \
                >>{log} 2>&1
        
            # bcftools reheader \
            #     -s {output.new_name} \

        # bcftools stats -s {params.sample} {output.genos_depth_filter} > {output.stats} 2>>{log}
        exitcode=$? 
        if [ $exitcode -eq 1 ] ; then exit 1; else exit 0 ; fi
        """

rule bcf_to_vcf:
    input:
        bcf = "Filtered/{bamfile}_chr{chr}_depth_filter_{rmTrans}.bcf"
    output:
        vcf = "Filtered/{bamfile}_chr{chr}_depth_filter_{rmTrans}.vcf"
    shell:
        """
        bcftools view -Ov -o {output.vcf} {input.bcf}        
        """
#-----------------------------------------------------------------------------#
rule recode_plink_chr:
    input:
        vcf = "Filtered/{bamfile}_chr{chr}_depth_filter_{rmTrans}.vcf"
    output:
        #vcf = "Filtered/{bamfile}_{chr}_depth_filter.vcf",
        bed_ind = "Filtered/{bamfile}_chr{chr}_depth_filter_{rmTrans}.bed",
        bim_ind = "Filtered/{bamfile}_chr{chr}_depth_filter_{rmTrans}.bim"
    shell:
        """        
        #bcftools view -Ov -o  
        plink --recode --vcf {input.vcf} \
            --out Filtered/{wildcards.bamfile}_{wildcards.chr}_depth_filter_{wildcards.rmTrans} \
            --make-bed --double-id --set-missing-var-ids @:#
        """

#=====================================================
