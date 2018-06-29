# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 12:39:34 2018

@author: Angela
"""

#make "chromsome painting" input from LAMP output
    #was too slow in R
    
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--ancestry", type = str, action = "store", dest = "ancestry", help = "Path to ancestry.txt")
parser.add_argument("--pos", type = str, action = "store", dest = "pos", help = "Path to pos.txt")
parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", help = "Prefix for output .csv for R input")
args = parser.parse_args()

print("Starting conversion of LAMP output to ggplot-able input.\nLoading ancestry input.")
ancestry = pd.read_table(args.ancestry, sep = "\t", header = None).transpose().dropna()
print("Loaded ancestry file.\nLoading position file.")
pos = pd.read_table(args.pos, header = None).dropna()
print("Loaded position file.")

pop_ancs = pd.DataFrame(index = range(len(ancestry) - 1), columns = range(len(ancestry.columns)/3))
start = 0
print("Starting parsing of individual ancestries.")
for prob in range(0, len(ancestry.columns)/3): 
    person = ancestry.iloc[:, start:(start + 3)] #split ancestry file into df's of 3 cols each (3 ref. pops)
    person_ID = person.iloc[0, 0].replace(":", "") #make id w/o :
    person = person.iloc[1:] #remove id line
    personal_anc = []
    for SNP in range(0, len(person)): #most likely ancestry by SNP
        if person.iloc[SNP, 0] > (person.iloc[SNP, 1] + person.iloc[SNP, 2]):
            personal_anc.append("IBS")
        elif person.iloc[SNP, 1] > (person.iloc[SNP, 0] + person.iloc[SNP, 2]):
            personal_anc.append("PEL")
        elif person.iloc[SNP, 2] > (person.iloc[SNP, 0] + person.iloc[SNP, 1]):
            personal_anc.append("YRI")
        else:
            personal_anc.append("NA")
    pop_ancs[prob] = personal_anc #add to data frame of all
    pop_ancs.columns.values[prob] = person_ID
    start = start + 3
    print("Finished with person " + str(prob + 1) + " out of " + str(len(ancestry.columns)/3) + ".")
pop_ancs.to_csv(args.output_prefix + ".csv", na_rep = "NA", header = True, index = False)
print("Done with conversion. Exiting program. Have a nice day :).")

