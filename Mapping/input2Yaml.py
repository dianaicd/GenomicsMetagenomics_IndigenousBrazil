#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:35:40 2019

@author: dcruz
"""
# %%
myInput = "/Users/dcruz/Projects/Botocudos/Files/test/24samples.txt"
column = "LB"
# %%
def get_index_column(input, columnName, sep = " "):
    with open(input, 'r') as file:
        header = file.readline().replace("\n", "")
        try:
            columnIndex = [index for index,element in enumerate(header.split(sep)) if(element == columnName)][0]
        except:
            columnIndex = False
    return(columnIndex)

# %%
def is_paired(line, indexColumn, sep = " "):
    if line.split(sep)[indexColumn] == "NULL":
        return("Single-end")
    else:
        return("Paired-end")

# %%
def get_values_column(input, column, sample = False, lib = False, sep = " "):
    # Function to return unique values per column
    myValues = {}
    with open(input, 'r') as file:
        header = file.readline().replace("\n", "")
        if column in header.split(sep):
            indexColumn = [index for index,element in enumerate(header.split(sep))
                    if(column == element)][0]
            
            colSM = get_index_column(input, "SM", sep = sep) #[index for index,element in enumerate(header.split(sep)) if(element == "SM\n")][0]
            colLB = get_index_column(input, "LB", sep = sep) #[index for index,element in enumerate(header.split(sep)) if(element == "LB")][0]

            for line in file.readlines():
                line = line.replace("\n", "")
                if sample:
                    if sample ==  line.split(sep)[colSM]:
                        if lib:
                            if lib == line.split(sep)[colLB]:
                                if column == "Data2":
                                    myValues[is_paired(line, indexColumn, sep)] = 1
                                else:
                                    myValues[line.split(sep)[indexColumn]] = 1
                        else:
                            myValues[line.split(sep)[indexColumn]] = 1
                else:
                    myValues[line.split(sep)[indexColumn]] = 1
                # Ask for column Data2, but it is not in the input file; libraries are then single-end
        elif column == "Data2":
            myValues["Single-end"] = 1
        else:
            myValues["NULL"] = 1
        return(list(myValues.keys()))

# %%
get_values_column(myInput, "ID")
depth = "ID"
# %%
def expand_input(input, depth = "ID", sep = " "):
    # returns prefix for an input
    myPrefix = {}
    with open(input, 'r') as file:
        header = file.readline().replace("\n", "")
        smIndex = get_index_column(input, "SM", sep)
        lbIndex = get_index_column(input, "LB", sep)
        idIndex = get_index_column(input, "ID", sep)
        for line in file.readlines():
            sm = line.split(sep)[smIndex].replace("\n", "")
            lb = line.split(sep)[lbIndex]
            id = line.split(sep)[idIndex]
            if depth == "SM":
                myPrefix[sm+"/"+sm] = 1
            elif depth == "LB":
                myPrefix[sm+"/"+lb] = 1
            elif depth == "ID":
                myPrefix[sm+"/"+lb+"/"+id+"/"+id] = 1
    return(list(myPrefix.keys()))
#%%
expand_input(myInput, depth = "LB")

# %%
def ident_me(text, nIdent, colon = True):
    ident = "  " #ident with two spaces
    myString = ""
    for i in range(0, nIdent):
        myString = myString + ident
    myString = myString + text
    if colon:
        myString = myString + ":\n"
    else:
        myString = myString + "\n"
    return(myString)
# %%
def input2yaml(input, sep = " "):
#%%
    yaml = input.replace(".txt", ".yaml")
    SMs = get_values_column(input, "SM", sep = sep)
    #has_paired_data = get_index_column(myInput, "Data2", sep = sep)
    #if has_paired_data:
    
    with open(yaml, 'w') as YAML:
        YAML.write(ident_me("Samples", 0))
        for sm in SMs:
            YAML.write(ident_me(sm, 1))
            LBs = get_values_column(input = input, column = "LB",
                                   sample = sm, sep = sep)
            YAML.write(ident_me("Libraries", 2))
            for lb in LBs:
                YAML.write(ident_me(lb, 3))
                lib_type = get_values_column(input = input, column = "Data2",
                                   sample = sm, lib = lb, sep = sep)
                YAML.write(ident_me("Type", 4))  
                YAML.write(ident_me(lib_type[0], 5, False))               
                IDs = get_values_column(input = input, column = "ID",
                                   sample = sm, lib = lb, sep = sep)
                YAML.write(ident_me("IDs", 4))
                for id in IDs:
                    YAML.write(ident_me(id, 5, False))
        YAML.close()