configfile: "merge_plink_files.yaml"
import os

wildcard_constraints:
    file_1 = "(" + "|".join([f for f in list(config["files"].keys())])+ ")(?!diallelic)",
    file_2 =  "(" + "|".join([f for f in list(config["files"].keys())])+ ")(?!diallelic)"

rule all:
    input:
        bed = config["output"]
ruleorder: keep_autosomes > exclude_duplicates > exclude_triallelic > merge_bfiles
# include autosomes only
rule keep_autosomes:
    input:
        bed = lambda wildcards: config["files"][wildcards.file]["prefix"] + ".bed"
    output:
        bed = temp("{file}_autosomes.bed"),
        bim = temp("{file}_autosomes.bim"),
        fam = temp("{file}_autosomes.fam")
    log:
        "{file}_autosomes.log"
    params:
        input_basename = lambda wildcards: config["files"][wildcards.file]["prefix"],
        output_basename = "{file}_autosomes"
    shell:
        """
        plink --bfile {params.input_basename} \
            --make-bed --autosome \
            --out {params.output_basename} --alleleACGT
        """

# remove variants with duplicated positions
rule exclude_duplicates:
    input: 
        bim = "{file}.bim",
        bed = "{file}.bed",
        fam = "{file}.fam"
    output:
        dup_ids = "{file}_dup_ids.txt",
        vars_to_exclude = "{file}_dups.txt",
        bed = temp("{file}_unique.bed"),
        bim = temp("{file}_unique.bim"),
        fam = temp("{file}_unique.fam")
    log:
        "{file}_unique.log"
    params:
        basename = "{file}",
        out_basename = "{file}_unique"
    run:
        with open(input.bim, "r") as bim_file:
            vars = {line.split()[1]: [line.split()[0], line.split()[3]]
            for line in bim_file.readlines()}
        with open(input.bim, "r") as bim_file:
            all_ids = [line.split()[1]  for line in bim_file.readlines()]

        duplicates = {}
        all_pos = {}
        print("finding duplicated position")
        for ID in vars:
            pos = vars[ID][0] + "_" + vars[ID][1]
            if pos in all_pos:
                all_pos[pos].append(ID)
                duplicates[pos] = all_pos[pos]
            else:
                all_pos[pos] = [ID]
        print("finding duplicated ids")
        find_dup_id_command = f"""\
        cat {input.bim} |cut -f2 |\
            sort |uniq -c |grep -v "1 "|\
                sed 's/ \+/ /g'|\
                    cut -f3 -d " " > {output.dup_ids}
        """
        os.system(find_dup_id_command)

        print("writing duplicated positions")
        with open(output.vars_to_exclude, "w") as vars_to_exclude:
            vars_to_exclude.writelines([
                ID + "\n"
                for  pos in duplicates
                for ID in duplicates[pos]
            ])
        print("writing duplicated ids")
        os.system(f"cat  {output.dup_ids} >> {output.vars_to_exclude}")
        
        plink_command = f"""plink --bfile {params.basename} \
            --out {params.out_basename} --make-bed \
            --exclude {output.vars_to_exclude} \
          """
        os.system(plink_command)

# update IDs to chr_pos
rule update_ids:
    input:
        bim = "{file}_autosomes_unique.bim",
        bed = "{file}_autosomes_unique.bed",
        fam = "{file}_autosomes_unique.fam",
    output:
        bed = temp("{file}_newID.bed"),
        bim = temp("{file}_newID.bim"),
        fam = temp("{file}_newID.fam"),
        new_ids = "{file}_newID.txt"
    params:
        basename = "{file}_autosomes_unique",
        out_basename = "{file}_newID",
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{print $2,$1"_"$4}}' {input.bim} > {output.new_ids}
        plink --update-name {output.new_ids} \
            --make-bed --out {params.out_basename} \
            --bfile {params.basename}
        """
# exclude triallelic variants
# and reduce the sets to the intersection of variants
rule exclude_triallelic:
    input:
        bim_1 = "{file_1}_newID.bim",
        bed_1 = "{file_1}_newID.bed",
        fam_1 = "{file_1}_newID.fam",
        bim_2 = "{file_2}_newID.bim",
        bed_2 = "{file_2}_newID.bed",
        fam_2 = "{file_2}_newID.fam"
    output:
        vars_to_extract_1 = "{file_1}_{file_2}_inter_diallelic.txt",
        vars_to_extract_2 = "{file_2}_{file_1}_inter_diallelic.txt",
        bed_1 = temp("{file_1}_{file_2}_diallelic.bed"),
        bim_1 = temp("{file_1}_{file_2}_diallelic.bim"),
        fam_1 = temp("{file_1}_{file_2}_diallelic.fam"),
        bed_2 = temp("{file_2}_{file_1}_diallelic.bed"),
        bim_2 = temp("{file_2}_{file_1}_diallelic.bim"),
        fam_2 = temp("{file_2}_{file_1}_diallelic.fam"),
        intersection = "intersection_{file_1}_{file_2}.txt"
    params:
        basename_1 = "{file_1}_newID",
        basename_2 = "{file_2}_newID",
        out_basename_1 = "{file_1}_{file_2}_diallelic",
        out_basename_2 = "{file_2}_{file_1}_diallelic"
    run:
        with open(input.bim_1, 'r') as bim_file:
            vars_1 = {line.split()[1]: [line.split()[4], line.split()[5], line.split()[1]]
                    for line in bim_file.readlines()}

        with open(input.bim_2, 'r') as bim_file:
            vars_2 = {line.split()[1]: [line.split()[4], line.split()[5], line.split()[1]]
                    for line in bim_file.readlines()}
        #---------------------------------------------------------
        # find triallelic sites
        triallelic_1 = {}
        triallelic_2 = {}
        with open(output.intersection, "w") as inter:
            for ID in vars_1:
                if ID in vars_2:

                    inter.write("\t".join([ID,
                         vars_1[ID][0], vars_1[ID][1], vars_1[ID][2] ,
                          vars_2[ID][0], vars_2[ID][1], vars_2[ID][2]]) + "\n")

                    vars = set(vars_1[ID][0:2] + vars_2[ID][0:2])
                    vars.discard('0')
                    if len(vars) > 2:
                        triallelic_1[vars_1[ID][2]] = 1
                        triallelic_2[vars_2[ID][2]] = 1
        #---------------------------------------------------------
        # get intersection and discard triallelics
        var1_var2 = set([vars_1[ID][2] for ID in set([i for i in vars_1]).intersection(set([i for i in vars_2]))])
        [var1_var2.discard(i) for i in set([i for i in triallelic_1])]

        var2_var1 = set([vars_2[ID][2] for ID in set([i for i in vars_2]).intersection(set([i for i in vars_1]))])
        [var2_var1.discard(i) for i in set([i for i in triallelic_2])]

        with open(output.vars_to_extract_1, 'w') as f:
            f.writelines([
                ID + "\n" for ID in var1_var2
            ])

        with open(output.vars_to_extract_2, 'w') as f:
            f.writelines([
                ID + "\n" for ID in var2_var1
            ])    

        plink_command = f"""plink --bfile {params.basename_1} \
            --out {params.out_basename_1} --make-bed \
            --extract {output.vars_to_extract_1} \
            """
        os.system(plink_command)

        plink_command = f"""plink --bfile {params.basename_2} \
            --out {params.out_basename_2} --make-bed \
            --extract {output.vars_to_extract_2} \
            """
        os.system(plink_command)

# try to merge
rule merge_bfiles:
    input:
        bed_1 = "{file_1}_{file_2}_diallelic.bed",
        bim_1 = "{file_1}_{file_2}_diallelic.bim",
        fam_1 = "{file_1}_{file_2}_diallelic.fam",
        bed_2 = "{file_2}_{file_1}_diallelic.bed",
        bim_2 = "{file_2}_{file_1}_diallelic.bim",
        fam_2 = "{file_2}_{file_1}_diallelic.fam"
    output:
        bed = "{file_1}_{file_2}_merged.bed"
    params:
        basename_1 = "{file_1}_{file_2}_diallelic",
        basename_2 = "{file_2}_{file_1}_diallelic",
        out_basename = "{file_1}_{file_2}_merged"
    shell:
        """
        plink --bfile {params.basename_1} \
            --bmerge {params.basename_2} \
            --make-bed --out {params.out_basename} \
            --merge-equal-pos
        """