# Snakemake for Error analysis in ANGSD
myGit = "~/Projects/Botocudos/Scripts/"
config_name = "~/Projects/Botocudos/Files/test/"
def get_git_path(myGit, folder = ""):
    return(myGit + folder)

subworkflow make_bamlist:
    workdir:
        get_git_path(myGit, folder = "misc/")
    snakefile:
        get_git_path(myGit, folder = "misc/make_bamlist.smk")
    configfile:
        "dummy_samples.yaml"

