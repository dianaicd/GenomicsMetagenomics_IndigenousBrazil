configfile: "multiple_purposes.yaml"
import os, glob, math
include: "parse_resources.smk"


panels = list( config['roh_hapROH']['panels'].keys() ) 
bamlists = { panel:list( config['roh_hapROH']['panels'][panel]['bamlists'].keys() ) for panel in panels }
ind_pop_files = { bamlist:config['roh_hapROH']['panels'][panel]['bamlists'][bamlist]['ind_pop'] for panel in panels for bamlist in bamlists[ panel ]}

wildcard_constraints:
    bamlist = ")".join(["(?<!"+p for p in panels])  + "_)" + "|".join([ bamlist for panel in panels for bamlist in bamlists[ panel ]])


rule all:
    input:
        geno = ["{panel}/{bamlist}.geno".format( panel = panel,
                                                            bamlist = bamlist)
                for panel in panels for bamlist in bamlists[ panel ]
        ]


rule sample_alleles:
    input:
        snakefile = "source4merging.snakefile",
        # config = "multiple_purposes.yaml",
        parse_resources = "parse_resources.smk"
    output:
        merged = "{panel}/{bamlist}.counts.sampled"
    shell:
        """
        snakemake -prs source4merging.snakefile --until {output} --profile axiom
        """


# Make a .map file with four fields:
rule make_map:
    input:
        map = "{panel}.map"
    output:
        map = "{panel}/{bamlist}.map"
    shell:
        """
        cat {input.map} | sed 's/ \+/ /g' |cut -f1-4 -d ' ' > {output.map}
        """

# Must generate a tped file by pastinf counts.sampled
# and four additional fields (those that should be found in a .map file)
rule make_tped:
    input:
        sampled_alleles = "{panel}/{bamlist}.counts.sampled",
        map = "{panel}/{bamlist}.map"
    output:
        tped = "{panel}/{bamlist}.tped",
    shell:
        """
        paste {input.map} {input.sampled_alleles} -d ' ' > {output.tped} 
        """

rule make_fam:
    input:
        ids = "{bamlist}.ids"
    output:
        fam = "{panel}/{bamlist}.fam",
        tfam = "{panel}/{bamlist}.tfam",
        pedind =  "{panel}/{bamlist}.pedind"
    shell:
        """
        awk '{{print $1,$1,0,0,0,1}}' {input.ids} > {output.fam}
        cp {output.fam} {output.pedind}
        cp {output.fam} {output.tfam}
        """

rule recode_tped_bed:
    input:
        tped = "{panel}/{bamlist}.tped",
        tfam = "{panel}/{bamlist}.tfam",
        map = "{panel}/{bamlist}.map"
    output:
        bed = "{panel}/{bamlist}.bed"
    params:
        prefix = "{panel}/{bamlist}"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("recode_tped_mem", attempt, 8),
        runtime=lambda wildcards, attempt: get_runtime_alloc("recode_tped_time", attempt, 4), 
    shell:
        """
        plink --recode --make-bed --tfile {params.prefix} --out {params.prefix}
        """

# There will be an error when running convertf
# if the IDs exceed certain length; switch them to simple numbers
def match_ind_pop( bamlist):
    ind_pop = ind_pop_files[ bamlist ]
    population = {}

    with open(ind_pop, 'r') as file:
        for line in file.readlines():
            myId,myPop = line.replace("\n", "").split("\t")
            population[myId] = myPop 
    return(population)
    

rule make_short_ids:
    input:
        ind = "{panel}/{bamlist}.fam"
    output:
        ind = "{panel}/{bamlist}_newID.pedind",
        link_table = "{panel}/{bamlist}.link_ids"
    run:
        myPops = match_ind_pop( wildcards.bamlist )

        with open(input.ind, 'r') as fam, open(output.ind, 'w') as pedind, open(output.link_table, 'w') as link:
            for i,line in enumerate(fam.readlines()):
                fam,ind,dad,mom,sex,pheno = line.replace("\n","").split()
                pop = myPops[ind]
                pedind.write(" ".join([pop,str(i+1),dad,mom,sex,"1"]) + "\n")
                link.write("\t".join([fam,ind,pop,str(i+1)]) + "\n")

rule update_ids:
    input:
        link_table = "{file}.link_ids",
        bed = "{file}.bed"
    output:
        new_bed = "{file}_newID.bed",
        new_ped = "{file}_newID.ped",
        new_fam = "{file}_newID.fam",
        bim = "{file}_newID.bim",
        tmp_fam = temp("{file}_newID.tmp.fam")
    shell:
        """
        plink --make-bed --update-ids {input.link_table} --bfile {wildcards.file} --out {wildcards.file}_newID
        cp {output.new_fam} {output.tmp_fam}
        cat {output.tmp_fam} | awk '{{print $1,$2,$3,$4,$5,1}}' > {output.new_fam}
        plink --recode --bfile {wildcards.file}_newID --out {wildcards.file}_newID

        # plink --recode --make-bed --bfile {wildcards.file}_newID --out {wildcards.file}_newID
        """

rule par_ped_to_eigenstrat:
    input:
        ped = "{file}_newID.ped",
        pedind =  "{file}_newID.pedind",
        map = "{file}.map"
    output:
        par = temp("{file}_ped2packedancestry.par")
    shell:
        """
        argument=(genotype snpname indivname outputformat genotypeoutname 
                    snpoutname indivoutname familynames) 

        options=({input.ped} {input.map} {input.ped} PACKEDANCESTRYMAP {wildcards.file}.geno 
                {wildcards.file}.snp {wildcards.file}.ind YES) 
        
        for i in $(seq 0 7) 
        do 
          echo "${{argument[$i]}}: ${{options[$i]}}" >> {output.par} 
        done
        """

rule convert_eigenstrat:
    input:
        par = "{file}_ped2packedancestry.par"
    output:
        geno = "{file}.geno",
        snp = "{file}.snp",
        ind = "{file}.ind"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("convertf_mem", attempt, 16),
        runtime=lambda wildcards, attempt: get_runtime_alloc("convertf_time", attempt, 8), 
    shell:
        """
        convertf -p {input.par}
        """

