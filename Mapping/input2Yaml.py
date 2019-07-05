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
def get_values_column(input, column, sample = False, lib = False):
    # Function to return unique values per column
    myValues = {}
    with open(input, 'r') as file:
        header = file.readline()
        indexColumn = [index for index,element in enumerate(header.split())
                if(column == element)][0]
        for line in file.readlines():
            if sample:
                if sample ==  line.split()[5]:
                    if lib:
                        if lib == line.split()[3]:
                            myValues[line.split()[indexColumn]] = 1
                    else:
                        myValues[line.split()[indexColumn]] = 1
            else:
                myValues[line.split()[indexColumn]] = 1

    return(list(myValues.keys()))


# %%

get_values_column(myInput, "ID")
depth = "ID"
# %%
def expand_input(input, depth = "ID"):
    # returns prefix for an input
    myPrefix = {}
    with open(input, 'r') as file:
        header = file.readline()
        for line in file.readlines():
            sm = line.split()[5]
            lb = line.split()[3]
            id = line.split()[0]

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
def identMe(text, nIdent, colon = True):
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
def input2yaml(input):
#%%
    yaml = input.replace(".txt", ".yaml")
    SMs = get_values_column(input, "SM")
    with open(yaml, 'w') as YAML:
        YAML.write(identMe("Samples", 0))
        for sm in SMs:
            YAML.write(identMe(sm, 1))
            LBs = get_values_column(input = input, column = "LB",
                                   sample = sm )
            YAML.write(identMe("Libraries", 2))
            for lb in LBs:
                YAML.write(identMe(lb, 3))
                IDs = get_values_column(input = input, column = "ID",
                                   sample = sm, lib = lb)
                YAML.write(identMe("IDs", 4))
                for id in IDs:
                    YAML.write(identMe(id, 5, False))
        YAML.close()