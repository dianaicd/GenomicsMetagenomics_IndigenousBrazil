# Snakefile version to merge datasets
configfile: "multiple_purposes.yaml"
import os, glob, math
include: "parse_resources.smk"

bamlists = list(config['merge_sample']["bamlists"].keys())
panels = list(config['merge_sample']["panels"].keys()) 
git_path = "/users/dcruzdav/data/Git/Botocudos-scripts/"

chromosomes = [str(x) for x in range(1,23)] 
mind = [str(round(x*0.01, 2)) for x in range(95, 100, 5)]

wildcard_constraints:
    extension = "(bed|vcf|ped)",
    Chr =   "|".join(chromosomes)  ,
    panel =  "(" + "|".join([p for p in panels])+ ")(?!_" + "|_".join([b for b in bamlists]) + ")" ,
    bamlist = ")".join(["(?<!"+p for p in panels])  + "_)" + "|".join([b for b in bamlists])


rule all:
    input:
        haploid_ped = expand("{panel}/{bamlist}_{panel}.haploid.ped", bamlist = bamlists,
                            panel = panels),
        dist = expand("{panel}/{bamlist}_{panel}.haploid.dist.gz", bamlist = bamlists,
                             panel = panels),
        dist_mind = expand("{panel}/{bamlist}_{panel}.mind{mind}.haploid.dist.gz", bamlist = bamlists,
                             panel = panels, mind = mind)#,


rule link_programs:
    input:
        count_and_sample = git_path + "AlleleCounts/count_and_sample.py",
        sample_ped = git_path + "MDS/sample_ped.pl",
        genolike = git_path + "GenoLike/workflow_genolike.sh"
    output:
        count_and_sample ="count_and_sample.py", 
        sample_ped = "sample_ped.pl", 
        genolike = "workflow_genolike.sh"
    shell:
        "ln -s {input.count_and_sample} ./ ;"
        "ln -s {input.sample_ped} ./ ;"
        "ln -s {input.genolike} ./"

def expand_path(bamlist):
    paths = list(config["merge_sample"]["bamlists"][bamlist].values())
    full_paths = [os.path.expanduser(p) for p in paths]
    bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

# def expand_bamlists(bamlist):  

# rule make_bamlists:
#     input:
#         lambda wildcards: expand_path(wildcards.bamlist)
#     output:
#         "{bamlist}.txt"
#     log:
#         "logs/{bamlist}_make_bamlist.log"
#     run:
#         with open(output[0], 'w') as file:
#             for line in input:
#                 file.write(line+"\n")

bamlist_nFiles = {}
bamlist_nInd = {}

for bamlist in bamlists:
    input = expand_path(bamlist)
    output = bamlist + ".txt"
    bamlist_nFiles[ bamlist ] = math.ceil( len(input) / 10 )
    bamlist_nInd[ bamlist ] = len( input )
    n_file = 1

    with open(output, 'w') as file:
        for line in input:
            file.write( line + "\n" )
    split_command = "split -l 10 --numeric-suffixes=1 --additional-suffix='.txt' " + output + " " + output.replace(".txt", "_") 
    # print(split_command)
    os.system(split_command)

def return_bai(bamlist):
    bai = [line + ".bai" for line in expand_path(bamlist)]
    return(bai)

rule samtools_index:
    input:
        bam = "{file}.bam"
    output:
        bai = "{file}.bam.bai"
    log:
        "logs/index_samtools_{file}.log"
    shell:
        """
        samtools index {input.bam} &>{log}
        """

rule mpileup:
    input:
        bai = lambda wildcards: return_bai(wildcards.bamlist),
        bamlist = "{bamlist}_{nGroup}.txt",
        panel = lambda wildcards: config['merge_sample']["panels"][wildcards.panel]["path"] ,
        sites = "{panel}/{Chr}_sites.bed"
    output:
        "{panel}/{bamlist}_{nGroup}_{Chr}.mpileup"
    wildcard_constraints:
        bamlist = ")".join(["(?<!"+p for p in panels])  + "_)" + "|".join([b for b in bamlists])
    params:
        baseQ = config["BaseQuality"]
    log:
        "{panel}/logs/{bamlist}_{nGroup}_{Chr}.log"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mpileup_mem", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mpileup_time", attempt, 24), 
    shell:
        """
        samtools mpileup -r {wildcards.Chr} -Bl {input.sites} \
           -b {input.bamlist} -a -o {output} -Q {params.baseQ} \
           &> {log}  
        """

# Bamlists are split into groups of 10 individuals or less
rule merge_mpileup_samples:
    input:
        mpileup = lambda wildcards: 
                            ["{panel}/{bamlist}_{nGroup}_{Chr}.mpileup".format(
                                panel = wildcards.panel,
                                Chr = wildcards.Chr,
                                bamlist = wildcards.bamlist,
                                nGroup = str(i) if i >= 10 else "0" + str(i)
                                ) 
                                for i in range(1,
                                                bamlist_nFiles[ wildcards.bamlist ] + 1
                                                )
                                             
                            ]
                            
    output:
        tmp_mpileup = temp("{panel}/{bamlist}_{Chr}_tmp.mpileup"),
        mpileup = "{panel}/{bamlist}_{Chr}.mpileup",
        samples_order = "{panel}/{bamlist}_{Chr}_samples_order.txt"
    run:
        def write_clean_pileup( line, mpileup ):

            pileup_data = "\t".join( [ item 
                                        for i in range(0, n_files ) 
                                        for item in line.split("\t")[ start_points[i] : end_points[i] ]
                                        ] )

            pileup_data = "\t".join( line.split( "\t" )[0:3] ) + "\t" + pileup_data 
            
            mpileup.write( pileup_data )

        n_ind = bamlist_nInd[ wildcards.bamlist ]
        n_files = bamlist_nFiles[ wildcards.bamlist ]
        start_points = [ ( i*33 ) + 3 for i in range(0, n_files) ]
        end_points = [i - 3 for i in start_points[1:] ] + [ 3*(n_ind + n_files) ]

        command_paste = "paste " + " ".join([mpileup for mpileup in input.mpileup]) + " > " +output.tmp_mpileup
        os.system(command_paste)

        with open( output.tmp_mpileup, "r") as tmp_mpileup, open( output.mpileup, "w" ) as mpileup:
            [ write_clean_pileup( line, mpileup ) for line in tmp_mpileup ]

        with open( output.samples_order, "w") as samples:
            [ samples.write(line+"\n") for line in input.mpileup ]


panelExtension={"vcf":"vcf", "bed":"bim", "bim":"bim", "ped":"map"}
columns = {"vcf":"1,2,4,5", "bim":"1,4,5,6", "map":"1,4,5,6"}

rule refalt:
    input:
        panel = lambda wildcards: config['merge_sample']["panels"][wildcards.panel]["path"]
    params:
        extension = lambda wildcards: panelExtension[config['merge_sample']["panels"][wildcards.panel]["type"]],
        col = lambda wildcards: columns[panelExtension[config['merge_sample']["panels"][wildcards.panel]["type"]]]
    output:
        refalt = "{panel}/sites.refalt"
    log:
        "{panel}/logs/refalt.log"
    shell:
        'cut -f {params.col} {wildcards.panel}.{params.extension} '
        '   | sed "s/\\t/_/ ;s/A/0/g ; s/C/1/g; s/G/2/g ; s/T/3/g" > {output.refalt} 2>{log}'

rule refalt_chr:
    input:
        refalt = "{panel}/sites.refalt"
    output:
        refalt = "{panel}/{Chr}.refalt",
        sites = "{panel}/{Chr}_sites.bed"
    log:
        "{panel}/logs/refalt_{Chr}.log"
    shell:
        'cut -f 1,2 {input.refalt} |grep -P "^{wildcards.Chr}\_"| sed "s/_/\\t/"| '
        '   awk \'{{print($1"\\t"$2-1"\\t"$2)}}\' > {output.sites} ; '
        '   grep -P "^{wildcards.Chr}\_" {input.refalt} > {output.refalt} '
        '   2>{log}'

rule count_and_sample:
    input:
        mpileup = "{panel}/{bamlist}_{Chr}.mpileup",
        refalt = "{panel}/{Chr}.refalt",
        program = "count_and_sample.py"
    output:
        counts = "{panel}/{bamlist}_{Chr}.counts.gz",
        sampled = "{panel}/{bamlist}_{Chr}.sampled.gz"
    params:
        allmutations = config['merge_sample']["allmutations"]
    log:
        "{panel}/logs/{bamlist}_{Chr}_count_and_sample.log"
    shell:
        "python count_and_sample.py --mpileup {input.mpileup} "
        " --counts {output.counts} "
        " --sampled {output.sampled} "
        " --refalt {input.refalt} --ped {params.allmutations}"
        " 2>{log}"

def expand_chrs(wildcards):
    myInput = expand("{panel}/{bamlist}_{Chr}.sampled.gz", 
                    bamlist = wildcards.bamlist, Chr = chromosomes, panel = wildcards.panel)
    return(myInput)

rule merge_sampled:
    input:
        sampled = expand_chrs
    output:
        merged = "{panel}/{bamlist}.counts.sampled"
    log:
        "{panel}/logs/{bamlist}.log"
    shell:
        "zcat {input.sampled} > {output.merged} 2>{log}"

rule panel_to_tped:
    input:
        panel = lambda wildcards: config['merge_sample']["panels"][wildcards.panel]["path"]
    output:
        panel = "{panel}/{panel}.tped",
        tfam = "{panel}/{panel}.tfam"
    log:
        "{panel}/logs/panel_to_tped.log"
    shell:
        "extension=$(echo {input.panel}|rev |cut -f1 -d. |rev) ;"
        "if [ $extension == 'bed' ] ; then " 
        "   plink --recode transpose --bfile {wildcards.panel} --out {wildcards.panel}/{wildcards.panel} 2>{log} ;"
        "elif [ $extension == 'vcf' ] ; then"
        "   plink --recode --vcf {wildcards.panel} --out {wildcards.panel}/{wildcards.panel} 2>{log} ;"
        "elif [ $extension == 'ped' ] ; then "
        "   plink --recode transpose --file {wildcards.panel} --out {wildcards.panel}/{wildcards.panel} 2>{log} ;"
        "fi ;"

rule merge_bam_to_tped:
    input:
        panel_tped = "{panel}/{panel}.tped",
        panel_tfam = "{panel}/{panel}.tfam",
        counts = "{panel}/{bamlist}.counts.sampled",
        bamlist = "{bamlist}.txt"
    output:
        merged_tped = "{panel}/{bamlist}_{panel}.tped",
        merged_tfam = "{panel}/{bamlist}_{panel}.tfam",
        temp_fam = temp("{panel}/{bamlist}_{panel}.temp.fam"),
        merged_fam = "{panel}/{bamlist}_{panel}.fam"
    log:
        "{panel}/logs/{bamlist}.log"
    shell:
        """
        paste {input.panel_tped} {input.counts} -d ' ' > {output.merged_tped}  2>>{log} 
        cp {input.panel_tfam} {output.merged_tfam}     2>> {log} 

        while read line 
        do 
            name=$(basename $line .bam) 
            echo \"$name $name 0 0 0 -9\" >> {output.merged_tfam} 
        done < {input.bamlist}             2>> {log} 

        plink --recode --make-bed --tfile {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel} \
           --out {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel}         2>>{log};

        awk '{{print $1,$2,$3,$4,$5,1}}' \
           {output.merged_fam} > {output.temp_fam}         2>>{log}

        cp {output.temp_fam} {output.merged_fam}          &>>{log}
        """

rule tped_to_ped:
    input:
        tped = "{panel}/{bamlist}_{panel}.tped"
    output:
        ped = "{panel}/{bamlist}_{panel}.ped",
        map = "{panel}/{bamlist}_{panel}.map"
    log:
        "{panel}/logs/{bamlist}.log"
    shell:
        """
        plink --recode --tfile {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel} \
            --out {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel}       &>>{log}
        """


rule ped_to_haploid:
    input:
        ped = "{panel}/{bamlist}_{panel}.ped",
        map = "{panel}/{bamlist}_{panel}.map",
        program = "sample_ped.pl"
    output:
        haploid = "{panel}/{bamlist}_{panel}.haploid.ped",
        dist = "{panel}/{bamlist}_{panel}.haploid.dist.gz",
        map = "{panel}/{bamlist}_{panel}.haploid.map"
    log:
        "{panel}/logs/{bamlist}.log"
    shell:
        """
        perl sample_ped.pl -ped {input.ped} -out {output.haploid} &>>{log} ;
        cp {input.map} {output.map} &>> {log} ;
        plink --distance square gz 'flat-missing' \
            --file {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel}.haploid \
            --out {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel}.haploid \
                   &>>{log}
        """

rule calc_dist:
    input:
        ped = "{panel}/{bamlist}_{panel}.haploid.ped"
    output:
        dist = "{panel}/{bamlist}_{panel}.mind{mind}.haploid.dist.gz"
    wildcard_constraints:
    log:
        "{panel}/logs/{bamlist}_mind{mind}.log"
    shell:
        """
        plink --distance square gz 'flat-missing' \
            --file {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel}.haploid \
            --mind {wildcards.mind} \
            --out {wildcards.panel}/{wildcards.bamlist}_{wildcards.panel}.mind{wildcards.mind}.haploid \
                   &>>{log}
        """