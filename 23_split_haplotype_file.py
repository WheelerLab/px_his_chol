# -*- coding: utf-8 -*-
"""
Created on Mon Oct 01 20:48:02 2018

@author: Angela
"""

import numpy as np
import pandas as pd

chr = 22
phind = pd.read_table("HCHS_chr22.phind", header = None, delim_whitespace = True)
all_inds = list(phind[0])
genofile = pd.read_table("test_1000_genofile.22", header = None)
    #no built in separator so gotta do it myself

#split into each character because no delimiters
haplotypes = []
for genofile_row in genofile.itertuples():
    haplotypes.append(list(genofile_row[1]))
haplotypes = pd.DataFrame(haplotypes)
haplotypes.columns = all_inds

#okay so now we have a dataframe of haplotypes
split_inds = np.array_split(all_inds, 12) #is about 2k reasonable?
for ind_chunk in range(0, len(split_inds)):
    #subset to within chunks
    ind_chunk_include = split_inds[ind_chunk].tolist()
    ind_chunk_genofile = haplotypes[ind_chunk_include]
    ind_chunk_genofile = ind_chunk_genofile.transpose()
    ind_chunk_genofile = pd.DataFrame(ind_chunk_genofile.astype(str).apply(''.join)) #since to_csv() doesn't like running without a delimiter, we cheat a bit
    ind_chunk_genofile_l = ind_chunk_genofile.values.tolist()
    ind_chunk_phind = phind.loc[phind[0].isin(ind_chunk_include)]

    #write to chunks    
    ind_chunk_genofile.to_csv("test_" + str(ind_chunk) + "genofile." + str(chr), sep = ",", na_rep = "NA", header = False, index = False)
    ind_chunk_phind.to_csv("test_" + str(ind_chunk) + "_" + str(chr) + ".phind", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f')
   
    
