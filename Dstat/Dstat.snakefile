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
    h1_h2 = [subset for subset in itertools.combinations(set(h1+h2),2) if subset[0] != subset[1] and subset[0] != outgroup and subset[1] != outgroup]

    with open(output, 'w') as outFile:
        [outFile.write("\t".join([combination[0], combination[1], pop3, outgroup])+"\n") 
            for combination in h1_h2 
            for pop3 in h3 if pop3 not in combination and pop3 != outgroup]

if not os.path.exists(panel+".combinations"):
    combine_pops(pops_h1, pops_h2, pops_h3, pops_outgroup, panel + ".combinations")

wildcard_constraints:
    myRange =  "(" + "|".join([r for r in expand_range(panel)])+ ")"

#=============================================================================#
# Define rules
rule all:
    input:
        par_Dstat = "{panel}.qpDstat.par".format(panel = panel),
        combinations = "{panel}.combinations".format(panel = panel),
        big_result = "Results/{panel}_qpDstat.results".format(panel = panel)

#=============================================================================#
# Clean data and get things in the rigth format
# The genotypes should be converted to eigenstrat format

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
            #myFiles=({input.dstat})
            cat Results/*/{panel}_*.qpDstat.results |  grep -P "result:" | sed 's/ \+/\\t/g' | cut -f2- | sort -nrk6 > {output.result}
        """
