
configfile: "genotype_calling.yaml"
# Call genotypes on a medium/high coverage sample
# and merge to a panel
def param_is_defined(name, default_par = False):
    myParameter = config[name] if name in config.keys() else default_par
    return(myParameter)
#-----------------------------------------------------------------------------#
# Variables and functions to begin with
chromosomes = [str(x) for x in range(1,23)] 
chromosomes.append("X")

samples_names = [list(myPaths)[0] for mySamples in config["bamlists"].values() for myPaths in mySamples.values()]
bamfile = list(config["bamlists"].keys())
panels = list(config["panels"].keys()) 

# Awful regular expresions to deal with files bearing bamfile and panel name
wildcard_constraints:
    extension = "(bed|vcf|ped)",
    # Chr =   "|".join(chromosomes)  ,
    panel =  "(" + "|".join([p for p in panels])+ ")(?!_" + "|_".join([b for b in bamfile]) + ")" ,
    bamfile = ")".join(["(?<!"+p for p in panels])  + "_)" + "|".join([b for b in bamfile])


def expand_path(myBamfile):
    mySample = list(config["bamlists"][myBamfile]["paths"].keys())[0]
    path = config["bamlists"][myBamfile]["paths"][mySample]
    full_path = os.path.expanduser(path)
    #bams = [f for p in full_paths for f in glob.glob(p)]
    return(full_path)
#-----------------------------------------------------------------------------#
rule all:
    input:
        vcf = expand("Filtered/{bamfile}_{panel}_depth_filter.vcf",
                        bamfile = bamfile, panel = panels),
        bed_ind = expand("Filtered/{bamfile}_{panel}_depth_filter.bed",
                        bamfile = bamfile, panel = panels),
        # bed_ind = expand("Filtered/{bamfile}_{panel}_depth_filter.bed",
        #                     bamfile = bamlists, panel = panels),
        bed_merged = expand("Merged/{panel}_merged.bed",#"Merged/{panel}_{bamfile}.bed",
                            panel = panels, bamfile = bamfile)

rule get_positions:
    input:
        panel = "{panel}.bim"
    output:
        positions = "Raw/{panel}_{chr}_positions.txt"
    shell:
        """
        grep -P "^{wildcards.chr}\t" {input.panel} |cut -f1,4 > {output.positions}
        """

rule call_genos:
    input:
        bamfile = lambda wildcards: expand_path(wildcards.bamfile),
        positions = "Raw/{panel}_{chr}_positions.txt"
    output:
        raw_genos = "Raw/{bamfile}_{panel}_{chr}.bcf"
    log:
        "logs/{bamfile}_{panel}_{chr}.log"
    params:
        minMapQ=30,
        minBaseQ=20,
        ref = "/scratch/axiom/FAC/FBM/DBC/amalaspi/popgen/reference_human/hs.build37.1/hs.build37.1.fa",
        threads = 4,
        moreno2019 = param_is_defined("Genos_Moreno2019", "Yes")
    shell:
        """
            if [ {params.moreno2019} == "Yes" ]
            then
                samtools mpileup -q {params.minMapQ} \
                -t DP -C50 -uf {params.ref} \
                -r {wildcards.chr} -l {input.positions} \
                {input.bamfile} | bcftools call -f GQ -c \
                --threads {params.threads} \
                -Ob -o {output.raw_genos}  \
                2>{log}
            else
                bcftools mpileup -C 50 -q {params.minMapQ} -Q {params.minBaseQ} -a FMT/DP,SP \
                        -f {params.ref} --threads {params.threads} -r {wildcards.chr} \
                        -R {input.positions} -Ou  {input.bamfile} | \
                    bcftools annotate -c RPB | \
                    bcftools call --threads {params.threads} -c -V indels \
                        -Ob -o {output.raw_genos} \
                2>{log}
            fi
        """

rule filter_genos:
    input:
        raw_genos = "Raw/{bamfile}_{panel}_{chr}.bcf"
    output:
        filtered_genos = "Filtered/{bamfile}_{panel}_{chr}_allDepths.bcf"
    params:
        threads = 20,
        genoqual = 30,
        allelic_imbalance = 0.2,
        read_position_bias = 0,
        variant_distance_bias = 0,
        moreno2019 = param_is_defined("Genos_Moreno2019", "Yes")
    log:
        "logs/{bamfile}_{panel}_{chr}.log"
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


rule concat_genos:
    input:
        expand("Filtered/{bamfile}_{panel}_{chr}_allDepths.bcf", bamfile = "{bamfile}", chr = chromosomes, panel = "{panel}")
    output:
        genos_all_depths = "Filtered/{bamfile}_{panel}_allDepths.bcf"
    params:
        threads = 20
    log:
    shell:
        """
        bcftools concat {input} -Ob -o {output.genos_all_depths} --threads {params.threads}
        """

rule filter_depth:
    input:
        genos_all_depths = "Filtered/{bamfile}_{panel}_allDepths.bcf"
    output:
        genos_depth_filter = "Filtered/{bamfile}_{panel}_depth_filter.bcf",
        stats = "Filtered/{bamfile}_{panel}_stats.txt",
        new_name = "Filtered/{bamfile}_{panel}_new_sample_name.txt"
    params:
    # Minimum and maximum depth relative to average depth of coverage 
    # on filtered genotypes in this sample
        min_depth = "1/3",
        max_depth = "2",
        sample = lambda wildcards: list(config["bamlists"][wildcards.bamfile]["paths"].keys())[0]
    log:
        "logs/filter_depth_{bamfile}_{panel}.log"
    threads: 20
    shell:
        """
        set +e
        sample=$(bcftools query -f '[%SAMPLE]\n' {input.genos_all_depths} \
                    |head -n1) >{log} 2>&1

        avgdp=$(bcftools stats -s $sample {input.genos_all_depths} \
            |grep -P "^PSC" |cut -f 10); echo "average depth: $avgdp " >>{log} 2>&1

        mindp=$avgdp*{params.min_depth} >>{log} 2>&1
        maxdp=$avgdp*{params.max_depth} >>{log} 2>&1

        echo "$sample {params.sample}" > {output.new_name}  2>>{log}
        bcftools filter --threads {threads} --exclude \
            "(SUM(DP4)< $mindp| SUM(DP4) > $maxdp )" {input.genos_all_depths} \
            -Ob |bcftools reheader -s {output.new_name} -o {output.genos_depth_filter} \
            >>{log} 2>&1
        
        bcftools stats -s {params.sample} {output.genos_depth_filter} > {output.stats} 2>>{log}
        exitcode=$? 
        if [ $exitcode -eq 1 ] ; then exit 1; else exit 0 ; fi
        """

#-----------------------------------------------------------------------------#
rule recode_plink:
    input:
        bcf = "Filtered/{bamfile}_{panel}_depth_filter.bcf"
    output:
        vcf = "Filtered/{bamfile}_{panel}_depth_filter.vcf",
        bed_ind = "Filtered/{bamfile}_{panel}_depth_filter.bed",
        bim_ind = "Filtered/{bamfile}_{panel}_depth_filter.bim"
    shell:
        """        
        bcftools view -Ov -o {output.vcf} {input.bcf}
        plink --recode --vcf {output.vcf} \
            --out Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter \
            --make-bed --double-id #--set-missing-var-ids @:#
        """

rule update_alleles:
    input:
        sample_bim = "Filtered/{bamfile}_{panel}_depth_filter.bim",
        panel_bim = "{panel}.bim"
    output:
        alleles_to_change = "Filtered/{bamfile}_{panel}_change_alleles.txt",
        alleles_to_remove = "Filtered/{bamfile}_{panel}_rm_alleles.txt",
        new_bed = "Filtered/{bamfile}_{panel}_depth_filter_refalt.bed",
        old_bim = "Filtered/{bamfile}_{panel}_depth_filter_old.bim"
    shell:
        """
        mv {input.sample_bim} {output.old_bim}
        python ~/data/Git/Botocudos-scripts/Genotypes/match_alleles_to_change.py \
            --sample_bim {output.old_bim} \
            --panel_bim {input.panel_bim} \
            --output_alleles {output.alleles_to_change} \
            --output_bim {input.sample_bim}

        if [ $(grep -cP "0\t0$" {output.alleles_to_change}) -eq 0 ]
        then
            echo "no alleles to remove"
            touch {output.alleles_to_remove}
            
            plink --bfile Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter \
            --update-alleles {output.alleles_to_change} \
            --make-bed --out Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter_refalt
        else
            echo "will remove $(grep -cP "0\t0$" {output.alleles_to_change}) alleles"
            grep -P "0\t0$" {output.alleles_to_change} > {output.alleles_to_remove}

            plink \
                --bfile Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter \
                --exclude {output.alleles_to_remove} \
                --make-bed \
                --out Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter_diallelic

            grep -vP "0\t0$" {output.alleles_to_change} > {output.alleles_to_change}.tmp
            mv {output.alleles_to_change}.tmp {output.alleles_to_change}

            plink --bfile Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter_diallelic \
            --update-alleles {output.alleles_to_change} \
            --make-bed --out Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter_refalt
        fi


        """

rule make_merge_list:
    input:
        bed_ind = expand("Filtered/{bamfile}_{panel}_depth_filter_refalt.bed",
                         panel = "{panel}", bamfile = bamfile) 
    output:
        list_to_merge = "Filtered/{panel}_to_merge.txt"
    run:
        with open(output[0], 'w') as myList:
            for line in input:
                myList.write(line.replace(".bed"," ") + "\n")

rule merge_plink:
    input:
        list_to_merge = "Filtered/{panel}_to_merge.txt", #"Filtered/{bamfile}_{panel}_depth_filter_refalt.bed",
        panel = "{panel}.bed"
    output:
        bed_merged = "Merged/{panel}_merged.bed" #"Merged/{panel}_{bamfile}.bed"
    log:
    params:
    wildcard_constraints:
        panel = ".*",
        bamfile = ".*"
    shell:
        """
        plink --bfile {wildcards.panel} \
            --merge-list {input.list_to_merge} \
            --out Merged/{wildcards.panel}_merged
        """
            # Filtered/{wildcards.bamfile}_{wildcards.panel}_depth_filter_refalt \
 #{wildcards.bamfile}
#=====================================================
