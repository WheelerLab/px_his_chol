# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 10:09:04 2018

@author: Angela
"""
import os
#extract YRI, CEU, and PEL from 1000G to use as anchors in Hispanic PCA
open_dir = "/home/angela/"
#open_dir = "Z:/"
super_pop = ["AFR", "AMR", "EUR"]
sub_pop = ["YRI", "PEL", "IBS"]

anchor_pops = open(open_dir + "px_his_chol/ADMIXTURE/anchor_pops.txt", "w")
for select_pop in range(0, len(super_pop)):
    fam = open(open_dir + "1000G/fams/" + super_pop[select_pop] + "_fam.txt", "r")
    for line in fam:
        arr = line.strip().split()
        (FID, IID, pop) = arr[0:3]
        if pop == sub_pop[select_pop]:
            anchor_pops.write(FID + "\t" + IID + "\n")

os.system("plink --bfile " + open_dir + "1000G/1000G --keep " + open_dir + "px_his_chol/ADMIXTURE/anchor_pops.txt --make-bed --out " + open_dir + "px_his_chol/ADMIXTURE/anchor_pops")
