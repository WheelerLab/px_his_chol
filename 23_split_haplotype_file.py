#MOSAIC input with all haplotypes is too large, so this program chunks the ~24k haplotypes into 12 equal chunks for an easier time on the memory
import argparse
import numpy as np
import pandas as pd
import os
parser = argparse.ArgumentParser()

parser.add_argument("--input_path", type = str, action = "store", dest = "input_path", required = False, default = "", help = "Path to folder containing input for MOSAIC")
parser.add_argument("--phind_prefix", type = str, action = "store", dest = "phind_prefix", required = False, default = "HCHS_chr", help = "Prefix of MOSAIC phind input")
parser.add_argument("--genofile_prefix", type = str, action = "store", dest = "genofile_prefix", required = False, default = "HCHSgenofile", help = "Rrefix of MOSAIC genofile input")
parser.add_argument("--output_path", type = str, action = "store", dest = "output_path", required = False, default = "", help = "Folder to output to")
parser.add_argument("--num_splits", type = str, action = "store", dest = "num_splits", required = False, default = "12", help = "Number of splits to make")
parser.add_argument("--pop_name", type = str, action = "store", dest = "pop_name", required = False, default = "HCHS", help = "Name of studied population")
args = parser.parse_args()

print("Reading in input.")
input_path = args.input_path
phind_prefix = args.phind_prefix
genofile_prefix = args.genofile_prefix
output_path = args.output_path
num_splits = int(args.num_splits)
pop_name = args.pop_name

'''
input_path = ""
phind_prefix = "HCHS_chr"
genofile_prefix = "HCHSgenofile"
output_path = "output/"
pop_name = "HCHS"
chr = 1
num_splits = 12
ind_chunk = 0
'''
#for writing sample.names for each chunk
sample_names = np.loadtxt(input_path + "sample.names", dtype = str).tolist()
ref_samples = good = [sample for sample in sample_names if sample != pop_name]

for chr in range(1, 23):
    print("Starting processes on chromosome " + str(chr) + ".")
    phind = pd.read_table(input_path + phind_prefix + str(chr) + ".phind", header = None, delim_whitespace = True)
    all_inds = list(phind[0])
    genofile_chunks = []
    for chunk in pd.read_table(input_path + genofile_prefix + "." + str(chr), header = None, chunksize = 20000):
        genofile_chunks.append(chunk) #when your data too thicc
    #no built in separator so gotta do it myself
    genofile = pd.concat(genofile_chunks, axis = 0)
    del genofile_chunks
    
    #split into each character because no delimiters
    haplotypes = []
    for genofile_row in genofile.itertuples():
        haplotypes.append(list(genofile_row[1]))
    haplotypes = pd.DataFrame(haplotypes)
    haplotypes.columns = all_inds

    #okay so now we have a dataframe of haplotypes
    split_inds = np.array_split(all_inds, num_splits) #is about 2k reasonable?
    for ind_chunk in range(0, len(split_inds)):
        #subset to within chunks
        ind_chunk_include = split_inds[ind_chunk].tolist()
        
        ind_chunk_genofile = haplotypes[ind_chunk_include]
        ind_chunk_genofile = ind_chunk_genofile.transpose()
        ind_chunk_genofile = pd.DataFrame(ind_chunk_genofile.astype(str).apply(''.join)) #since to_csv() doesn't like running without a delimiter, we cheat a bit
        ind_chunk_genofile_l = ind_chunk_genofile.values.tolist()
        ind_chunk_phind = phind.loc[phind[0].isin(ind_chunk_include)]
        ind_chunk_sample_names = ref_samples
        ind_chunk_sample_names.extend([pop_name] * len(ind_chunk_include))
        ind_chunk_sample_names = pd.DataFrame(ind_chunk_sample_names)

        #make chunk directory and copy snpfile and rate files
        os.system("mkdir -p " + output_path + "chr" + str(chr) + "_" + str(ind_chunk))
        os.system("cp " + input_path + "snpfile." + str(chr) + " " + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        os.system("cp " + input_path + "rates." + str(chr) + " " + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        
        #and need all ancestries files to run MOSAIC smoothly
        os.system("cp " + input_path + "IBS_chr" + str(chr) + "*" + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        os.system("cp " + input_path + "NAT_chr" + str(chr) + "*" + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        os.system("cp " + input_path + "YRI_chr" + str(chr) + "*" + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        os.system("cp " + input_path + "IBSgenofile." + str(chr) + "*" + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        os.system("cp " + input_path + "NATgenofile." + str(chr) + "*" + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        os.system("cp " + input_path + "YRIgenofile" + str(chr) + "*" + output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/")
        
        #write to chunks    
        ind_chunk_genofile.to_csv(output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/" + genofile_prefix + "." + str(chr), sep = ",", na_rep = "NA", header = False, index = False)
        ind_chunk_phind.to_csv(output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/" + phind_prefix + str(chr) + ".phind", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f')
        ind_chunk_sample_names.to_csv(output_path + "chr" + str(chr) + "_" + str(ind_chunk) + "/sample.names", sep = ",", na_rep = "NA", header = False, index = False)
        print("Complete with chunk " + str(ind_chunk + 1) + ".")
   
