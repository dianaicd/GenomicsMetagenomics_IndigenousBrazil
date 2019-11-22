# Subworkflow to make lists with bamfiles

# A config file should have a structure like:
# bamlists:
    # Boto_Posth_M2014_fV:
    #     paths:
    #         # Remember to add "/*bam" at the end
    #         Botocudos: "~/Project/Botocudos/BAM/*bam"
    #         Posth: "~/Project/Americas/Posth/BAM/*bam"
    #         Malaspinas2014: "~/Project/Botocudos/Malaspinas2014/*bam"
    #         fromVictor_Ancient: "~/Project/Americas/fromVictor/Ancient/*bam"
    #         fromVictor_Modern : "~/Project/Americas/fromVictor/Modern/*bam"
    # MN00010:
    #     paths:
    #         MN00010: '~/Project/Botocudos/BAM/MN00010.bam'

# etc...

# configfile: "samples_panel.yaml"

import glob,os

#=============================================================================#

def expand_path(wildcards):
    paths = list(config["bamlists"][wildcards.bamlist]["paths"].values())
    full_paths = [os.path.expanduser(p) for p in paths]
    bams = [f for p in full_paths for f in glob.glob(p)]
    return(bams)

def expand_path_groups(path): #bamlist, group):
    #path = config["bamlists"][bamlist]["paths"][group]
    full_paths = os.path.expanduser(path)
    bams = glob.glob(full_paths)
    return(bams)

rule make_bamlist:
    input:
        expand_path 
    output:
        "{bamlist}.txt"
    run:
        with open(output[0], 'w') as file:
            for line in input:
                file.write(line+"\n")
                
paths = {}
bamlists = list(config["bamlists"].keys())
def add_key(group,b):
    paths[group] = config["bamlists"][b]["paths"][group]


[add_key(group,b) for b in bamlists 
for group in list(config["bamlists"][b]["paths"].keys())]

rule make_bamlist_groups:
    input:
        bams = lambda wildcards: expand_path_groups(paths[wildcards.group])#lambda wildcards: expand_path_groups(wildcards.group)
    output:
        "Bams_group/{group}.txt"
    run:
        with open(output[0], 'w') as file:
            for line in input:
                file.write(line+"\n")
