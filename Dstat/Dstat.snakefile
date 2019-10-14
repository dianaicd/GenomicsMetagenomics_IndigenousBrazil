configfile: "Dstat.yaml"
import itertools
import re
# Rund qpDstat from AdmixTools

# Unfortunately, there is a lot of cleaning to do 
pops_h1 = config["H1"]
pops_h2 = config["H2"]
pops_h3 = config["H3"]
pops_outgroup = config["Outgroup"]
# This file should have a column with indId and
# one column with population ID
ind_pop = config["ind_pop"]

max_files_per_dir = config["max_files_per_dir"] if "max_files_per_dir" in config.keys() else 100
max_lines = config["max_lines"] if "max_lines" in config.keys() else 10

panel = [config["panels"][d]["path"].split(".")[0] for d in list(config["panels"].keys())][0]

#=============================================================================#
# Mess, mess, mess
def expand_range(panel):
    nLines =  len(open(panel + ".combinations").readlines())
    lower_bound = [i for i in range(1, nLines, max_lines)]
    upper_bound = [i-1 for i  in lower_bound[1:]]
    upper_bound.append(nLines)
    nRanges = range(0, len(lower_bound))
    myRanges = [str(lower_bound[i]) + "_" + str(upper_bound[i]) for i in nRanges]
    return(myRanges)

def expand_range_dir(panel):
    nLines =  len(open(panel + ".combinations").readlines())
    N = max_files_per_dir*max_lines
    lower_bound = [i for i in range(1, nLines, N)]
    upper_bound = [i-1 for i  in lower_bound[1:]]
    upper_bound.append(nLines)
    nRanges = range(0, len(lower_bound))
    myRanges = [str(lower_bound[i]) + "_" + str(upper_bound[i]) for i in nRanges]
    return(myRanges)

def expand_dstat(panel):
    nLines =  len(open(panel + ".combinations").readlines())
    sizeDir = max_files_per_dir*max_lines
    allDirs = []
    lower_bound = [i for i in range(1, nLines, sizeDir)]
    upper_bound = [i-1 for i  in lower_bound[1:]]
    upper_bound.append(nLines)
    def gimme_my_dir(i):
        l = [x for x in range(lower_bound[i], upper_bound[i], max_lines)]
        u = [x-1 for x  in l[1:]]
        u.append(upper_bound[i])
        myDirs = [ str(lower_bound[i]) + "_" + 
        str(upper_bound[i]) + "/" + 
        panel + "_" + str(l[j]) + "_" + 
        str(u[j]) for j in range(0, len(l))]
        return(myDirs)

    allDirs =[item for i in range(0, len(lower_bound)) for item in gimme_my_dir(i)]

    return(allDirs)



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

def combine_pops(pops_h1, pops_h2, pops_h3, outgroup, output):
    h1 = read_pops(pops_h1)
    h2 = read_pops(pops_h2)
    h3 = read_pops(pops_h3)
    h1_h2 = [subset for subset in itertools.product(h1,h2) if subset[0] != subset[1] ]

    with open(output, 'w') as outFile:
        [outFile.write("\t".join([combination[0], combination[1], pop3, outgroup])+"\n") 
            for combination in h1_h2 
            for pop3 in h3 if pop3 not in combination]
combine_pops(pops_h1, pops_h2, pops_h3, pops_outgroup, panel + ".combinations")

wildcard_constraints:
    myRange =  "(" + "|".join([r for r in expand_range(panel)])+ ")"

#=============================================================================#
# Define rules
rule all:
    input:
        par_conversion = "{panel}_ped2eigenstrat.par".format(panel = panel),
        par_Dstat = "{panel}.qpDstat.par".format(panel = panel),
        combinations = "{panel}.combinations".format(panel = panel),
        big_result = "Results/{panel}_qpDstat.results".format(panel = panel)

#=============================================================================#
# Clean data and get things in the rigth format
# The genotypes should be converted to eigenstrat format
# rule clean_fam:
# # last column should be 1
#     input:
#         "{panel}.fam"
#     output:
#         temp("{panel}.tmp.fam")
#     shell:
#         """
#         mv {input} {output}
#         cat {output} | awk 'BEGIN {{OFS="\t"}} {{print $1,$2,$3,$4,$5,1}}' > {input}
#         """

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


rule par_dstat:
    input:
        geno = "{panel}.eigenstratgeno",
        snp = "{panel}.snp",
        ind = "{panel}.mod.ind",
        pop_name = "{panel}.pops",
        pop_file = "{panel}.combinations"
    output:
        par = "{panel}.qpDstat.par"
    params:
        f4 = config["f4mode"] if "f4mode" in config.keys() else "NO"
    shell:
        """
        argument=(genotypename snpname indivname poplistname popfilename 
                    f4mode) 

        options=({input.geno} {input.snp} {input.ind} {input.pop_name} {input.pop_file} 
                {params.f4}) 
        
        for i in $(seq 0 5) 
        do 
          echo "${{argument[$i]}}: ${{options[$i]}}" >> {output.par} 
        done
        """

#=============================================================================#
# qpDstat may require a lot of memory, and it is handy to run only 
# a few combinations each time, in parallel
# Avoid a huge number of files in your directory

rule run_qpDstat:
    input:
        par = "{panel}.qpDstat.par",
        pop_file = "{panel}.combinations"
    output:
        result = "Results/{myFolder}/{panel}_{myRange}.qpDstat.results" # expand_range("{panel}")
    wildcard_constraints:
        myFolder = "(\d+)_(\d+)",
        myRange = "(\d+)_(\d+)"
    shell:
        """
        low=$(echo {wildcards.myRange} |cut -f1 -d '_')
        high=$(echo {wildcards.myRange} |cut -f2 -d '_')
        qpDstat -p {input.par} -l $low -h $high > {output.result} 
        """
 

rule wrap_qpDstat:
    input:
        dstat = expand("Results/{file}.qpDstat.results",#"Results/" + file + ".qpDstat.results" for file in expand_dstat(panel)]
                        #panel = panel, 
                        #myRange = expand_range(panel),
                        file = expand_dstat(panel))
    output:
        result = "Results/"+ panel + "_qpDstat.results"

    shell:
        """
            cat {input.dstat} |  grep -P "result:" | sed 's/ \+/\t/g' | cut -f2- | sort -nrk6 > {output.result}
        """
