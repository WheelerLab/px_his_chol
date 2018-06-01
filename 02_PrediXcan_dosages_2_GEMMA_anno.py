#Takes information from PrediXcan-style dosages to make SNP annotation files for GEMMA
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 09:31:31 2018

@author: Angela
"""

import gzip
dosage_path = "/home/angela/px_his_chol/Imputation/UMich/UMich_dosages/"
anno_path = "/home/angela/px_his_chol/GEMMA/anno/"

for i in range(1,23):
    anno = open(anno_path + "anno" + str(i) + ".txt", "w")
    for line in gzip.open(dosage_path + "chr" + str(i) + ".maf0.01.r20.8.dosage.txt.gz", "rb"):
        arr = line.strip().split()
        (chr, rs, bp) = arr[0:3]
        anno.write(rs + "\t" + bp + "\t" + chr + "\n")
    anno.close()
