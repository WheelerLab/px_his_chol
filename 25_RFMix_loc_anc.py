#Translates HAPI-UR and RFMix output into GEMMA-style input dosages to be used in 22_GEMMA_wrapper_admixture_mapping.py
#And yes this is just a subset of 19_MOSAIC_loc_anc.py but RFMix's output is a lot easier to parse
#test input: angela@wheelerlab3:~/px_his_chol/local_anc_GEMMA/RFMix_output$ python 25_RFMix_loc_anc.py --phind HCHS_chr22.phind --phsnp test.phsnp --Viterbi test.Viterbi.txt --output_prefix RFMix_for_GEMMA_test
import argparse
import numpy as np
import pandas as pd
import os

#input arguments
parser = argparse.ArgumentParser()
parser.add_argument("--admixed_hap", type = str, action = "store", dest = "admixed_hap", required = True, help = "Path to file of haplotype names including only admixed and no reference.")
parser.add_argument("--phsnp", type = str, action = "store", dest = "phsnp", required = True, help = "Path to .phsnp output by HAPI-UR")
parser.add_argument("--Viterbi", type = str, action = "store", dest = "Viterbi", required = True, help = "Path to .Viterbi. output by RFMix")
parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", required = False, default = "RFMix_for_GEMMA", help = "Prefix for GEMMA-input ancestry file")
args = parser.parse_args()

print("Reading input files.")
Viterbi_chunks = []
for chunk in pd.read_table(args.Viterbi, header = None, delim_whitespace = True, chunksize = 20000):
    Viterbi_chunks.append(chunk) #when your data too thicc
Viterbi = pd.concat(Viterbi_chunks, axis = 0).transpose()
del Viterbi_chunks
del chunk
haps = np.loadtxt(args.admixed_hap, dtype = 'string')
phsnp = pd.read_table(args.phsnp, header = None, delim_whitespace = True)
output_prefix = args.output_prefix

''' 
#test data
Viterbi_chunks = []
for chunk in pd.read_table("test.Viterbi.txt", header = None, delim_whitespace = True, chunksize = 20000):
    Viterbi_chunks.append(chunk) #when your data too thicc
Viterbi = pd.concat(Viterbi_chunks, axis = 0).transpose()
del Viterbi_chunks
del chunk
haps = np.loadtxt("haps.txt", dtype = 'string')
phsnp = pd.read_table("test.phsnp", header = None, delim_whitespace = True)
output_prefix = "test"
'''

#Convert input to suitable format for the loop
phsnp.columns = ['rs', 'chr', 'cM', 'bp', 'A1', 'A2']
SNPs = phsnp['rs'].tolist()
Viterbi.index = haps
Viterbi.columns = SNPs
ind_list = []
for hap in haps:
    ind_list.append(hap[:-2])
ind_list = sorted(set(ind_list), key = ind_list.index) #remove duplicates but preserve order
ind_file = open(output_prefix + "_ind.txt", "w") #to be used in subsetting GEMMA input
ind_file.write("\n".join(ind_list))
ind_file.close()

print("Starting to make dosage file.")
anc_dosage_write = open(output_prefix + ".csv", "a+")
anc_dosage_write.write("IID," + ",".join(SNPs) + "\n")
progress_landmarks_ind = np.linspace(0, len(ind_list), 21, dtype = int).tolist()
num_ind = 0
SNP_file = open(output_prefix + "_SNPs.txt", "w") #to be used in subsetting GEMMA input
SNP_file.write("\n".join(SNPs))
SNP_file.close()

for ind in ind_list:
    #Extract haplotypes
    ind_haplotype_A = []
    hap_A = ind + "_A"
    ind_haplotype_A = pd.DataFrame(Viterbi.loc[hap_A])
    ind_haplotype_B = []
    hap_B = ind + "_B"
    ind_haplotype_B = pd.DataFrame(Viterbi.loc[hap_B])
    loc_anc_dosage = []
    ind_anc = ind_haplotype_A.join(ind_haplotype_B, lsuffix = "_A", rsuffix = "_B") #marry hap A and hap B
    
    #first col is NAT, second is IBS, and third is YRI
    anc_dosage = []
    for ind_anc_row in ind_anc.itertuples(): #translate from RFMix codes to ancestry dosages
        if ind_anc_row[1] == 1 and ind_anc_row[2] == 1:
            anc_dosage.append([ind_anc_row[0], "0.0\t2.0\t0.0"]) #both IBS
        elif ind_anc_row[1] == 2 and ind_anc_row[2] == 2:
            anc_dosage.append([ind_anc_row[0], "2.0\t0.0\t0.0"]) #both NAT
        elif ind_anc_row[1] == 3 and ind_anc_row[2] == 3:
            anc_dosage.append([ind_anc_row[0], "0.0\t0.0\t2.0"]) #both YRI
        elif (ind_anc_row[1] == 1 and ind_anc_row[2] == 2) or (ind_anc_row[1] == 2 and ind_anc_row[2] == 1):
            anc_dosage.append([ind_anc_row[0], "1.0\t1.0\t0.0"]) #one NAT, one IBS
        elif (ind_anc_row[1] == 2 and ind_anc_row[2] == 3) or (ind_anc_row[1] == 3 and ind_anc_row[2] == 2): 
            anc_dosage.append([ind_anc_row[0], "1.0\t0.0\t1.0"]) #one NAT, one YRI
        elif (ind_anc_row[1] == 3 and ind_anc_row[2] == 1) or (ind_anc_row[1] == 1 and ind_anc_row[2] == 3):
            anc_dosage.append([ind_anc_row[0], "0.0\t1.0\t1.0"]) #one IBS, one YRI  
        else: #it shouldn't have to execute this ever
            anc_dosage.append([ind_anc_row[0], "NA\tNA\tNA"]) #who knows
    anc_dosage_df = pd.DataFrame(anc_dosage)
    anc_dosage_df.columns = ["rs", ind]
    anc_dosage_df = anc_dosage_df.drop_duplicates()
    anc_dosage_list = anc_dosage_df[ind].tolist()
    anc_dosage_write.write(ind + "," + ",".join(anc_dosage_list) + "\n") #write to dosage file
    
    '''
    If someone happens to read this far
    The reason I append to a file and not a data frame is b/c the input is >12k people
    Appending to a file is a constant speed
    Appending to a data frame has increasing speed as the data frame grows
    '''

    num_ind = num_ind + 1
    if num_ind in set(progress_landmarks_ind): #print progress by 5% increments
      progress = progress_landmarks_ind.index(num_ind)
      print("SNP ancestry dosage conversion is " + str(progress * 5) + "% complete.")
anc_dosage_write.close() #yay we're done!
print("Completed writing individual, SNP, and SNP ancestry dosage file to " + output_prefix + "_ind.txt, " + output_prefix + "_SNPs.txt and " + output_prefix + ".csv. Have a nice day!")
   
