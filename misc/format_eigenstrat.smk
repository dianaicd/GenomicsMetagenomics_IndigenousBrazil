# Snakemake to change a bed panel to eigenstrat format
configfile: "qp3pop_test.yaml"

# This file should have a column with indId and
# one column with population ID
ind_pop = config["ind_pop"]

panel = [config["panels"][d]["path"].split(".")[0] for d in list(config["panels"].keys())][0]

#=============================================================================#
# Mess, mess, mess

def match_ind_pop(panel = ind_pop):
    population = {}
    with open(ind_pop, 'r') as file:
        for line in file.readlines():
            myId,myPop = line.replace("\n", "").split("\t")
            population[myId] = myPop 
    return(population)
    

#=============================================================================#
# You must decide on one outgroup population (e.g., Yoruba)
# and three lists of populations for which gene flow will be tested
# I usually have a test of the form (O, NatAm1; NatAm2, Australasian)

def read_pops(fileName):
    with open(fileName, 'r') as myFile:
        myPops = [line.replace("\n", "") for line in myFile.readlines()]
    return(myPops)

#=============================================================================#
# Define rules
rule all:
    input:
        par_conversion = "{panel}_ped2eigenstrat.par".format(panel = panel)
#=============================================================================#
# Clean data and get things in the rigth format
# The genotypes should be converted to eigenstrat format

# There will be an error when running convertf
# if the IDs exceed certain length; switch them to simple numbers
rule make_short_ids:
    input:
        ind = "{panel}.fam"
    output:
        ind = "{panel}.pedind",
        link_table = "{panel}.link_ids"
    run:
        myPops = match_ind_pop()

        with open(input.ind, 'r') as fam, open(output.ind, 'w') as pedind, open(output.link_table, 'w') as link:
            for i,line in enumerate(fam.readlines()):
                fam,ind,dad,mom,sex,pheno = line.replace("\n","").split("\t")
                pop = myPops[ind]
                pedind.write(" ".join([str(i+1),str(i+1),dad,mom,sex,"1"]) + "\n")
                link.write("\t".join([fam,ind,pop,str(i+1)]) + "\n")

rule update_ids:
    input:
        link_table = "{panel}.link_ids",
        bed = "{panel}.bed"
    output:
        new_ped = "{panel}_newID.ped",
        new_bed = "{panel}_newID.bed",
        new_bim = "{panel}_newID.bim",
        new_fam = "{panel}_newID.fam",
        tmp_fam = temp("{panel}_newID.tmp.fam")
    shell:
        """
        plink --make-bed --update-ids {input.link_table} --bfile {wildcards.panel} --out {wildcards.panel}_newID
        cp {output.new_fam} {output.tmp_fam}
        cat {output.tmp_fam} | awk '{{print $1,$2,$3,$4,$5,1}}' > {output.new_fam}
        plink --recode --bfile {wildcards.panel}_newID --out {wildcards.panel}_newID
        """

rule par_ped_to_eigenstrat:
    input:
        ped = "{panel}_newID.ped",
        fam = "{panel}_newID.fam",
        bim = "{panel}_newID.bim"
    output:
        par = "{panel}_ped2eigenstrat.par"
    shell:
        """
        argument=(genotype snpname indivname outputformat genotypeoutname 
                    snpoutname indivoutname familynames) 

        options=({input.ped} {input.bim} {input.ped} EIGENSTRAT {wildcards.panel}.eigenstratgeno 
                {wildcards.panel}.snp {wildcards.panel}.ind YES) 
        
        for i in $(seq 0 7) 
        do 
          echo "${{argument[$i]}}: ${{options[$i]}}" >> {output.par} 
        done
        """

rule convert_eigenstrat:
    input:
        par = "{panel}_ped2eigenstrat.par"
    output:
        geno = "{panel}.eigenstratgeno",
        snp = "{panel}.snp",
        ind = "{panel}.ind"
    shell:
        """
        convertf -p {input.par}
        """

# For some reason we need this
rule replace_third_col:
    input:
        ind = "{panel}.ind"
    output:
        ind = "{panel}.mod.ind",
        pops = "{panel}.pops"
    run:
        myPops = {}
        with open(input.ind, 'r') as OldInd, open(output.ind, 'w') as NewInd, open(output.pops, 'w') as PopFile:
            for line in OldInd.readlines():
                pop,ind = line.split()[0].split(":")
                NewInd.write("\t".join([pop+":"+ind, "U", pop])+"\n")
                if pop not in myPops:
                    myPops[pop] = 1
            [PopFile.write(pop+"\n") for pop in myPops.keys()]
