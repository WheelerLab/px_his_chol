#MOSAIC input with all haplotypes is too large, so this program chunks the ~24k haplotypes into 12 equal chunks for an easier time on the memory
import argparse
import numpy as np
import pandas as pd
parser = argparse.ArgumentParser()

parser.add_argument("--input_path", type = str, action = "store", dest = "input_path", required = False, default = "", help = "Path to folder containing input for MOSAIC")
parser.add_argument("--phind_prefix", type = str, action = "store", dest = "phind_prefix", required = False, default = "HCHS_chr", help = "(Optional) prefix of MOSAIC phind input")
parser.add_argument("--genofile_prefix", type = str, action = "store", dest = "genofile_prefix", required = False, default = "HCHSgenofile", help = "(Optional) prefix of MOSAIC genofile input")
parser.add_argument("--output_path", type = str, action = "store", dest = "output_path", required = False, default = "output/", help = "(Optional) folder to output to")
args = parser.parse_args()

input_path = args.input_path
phind_prefix = args.phind_prefix
genofile_prefix = args.genofile_prefix
output_path = args.output_path

'''
input_path = ""
phind_prefix = "HCHS_chr"
genofile_prefix = "HCHSgenofile"
output_path = "output/"
chr = 1
'''

for chr in range(1, 23):
    phind = pd.read_table(input_path + phind_prefix + str(chr) + ".phind", header = None, delim_whitespace = True)
    all_inds = list(phind[0])
    genofile = pd.read_table(input_path + genofile_prefix + "." + str(chr), header = None)
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
        ind_chunk_genofile.to_csv(output_path + genofile_prefix + str(ind_chunk) + "." + str(chr), sep = ",", na_rep = "NA", header = False, index = False)
        ind_chunk_phind.to_csv(output_path + phind_prefix + str(chr) + "_" + str(ind_chunk) + ".phind", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f')
        
