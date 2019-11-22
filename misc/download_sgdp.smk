import os
with open("download-links-2.txt", "r") as file:
    def parse_line(line):
        if ".gz.tbi" in line:
            sufix = ".gz.tbi"
        else:
            sufix = ".gz"
        parsed = line.split("/")[-1].split(".gz")[0] + sufix
        return(parsed)
    myInput = [line for line in file.readlines()]
    myOutput = [parse_line(line) for line in myInput]
    myFiles = dict(zip(myOutput, myInput))

wildcard_constraints:
    file = "|".join([file for file in myOutput])
    
rule all:
    input:
        myFiles = myOutput

def gimme_input(file):
    return(myFiles[file].replace("&", "\&"))

rule download_things:
    output:
        file = "{file}"
    run:
        myCommand = "wget -c --no-use-server-timestamps -O "+ output.file + " " + gimme_input(output.file)
        os.system(myCommand)