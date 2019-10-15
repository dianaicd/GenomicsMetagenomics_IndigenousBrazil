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

config = "samples_panel.yaml"
import glob,os

rule all:
    input:
        bamlist = expand("{bamlist}.txt", bamlist = bamlists)

def expand_path(wildcards):
    paths = list(config["bamlists"][wildcards.bamlist]["paths"].values())
    full_paths = [os.path.expanduser(p) for p in paths]
    bams = [f for p in full_paths for f in glob.glob(p)]
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
