#uses local ancestry as a dosage to run an ancestry-by-ancestry level admixture mapping analysis.
#all input except for genotypes must be already in GEMMA format, and if some individuals and/or SNPs got removed, use 02_PrediXcan_dosages_to_GEMMA.py and 21_make_GEMMA_input.R to make input proper again
#cd /home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/
#python 22_GEMMA_wrapper_admixture_mapping.py --BIMBAM BIMBAM/chr22.txt.gz --anno anno/anno22.txt --relatedness relatedness.txt --pheno pheno.txt --covariates covariates.txt --snptable RFMix_for_GEMMA_22.csv --output chr22
import argparse
import numpy as np
import pandas as pd
import os
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
#novel-ish part of using GEMMA
parser.add_argument("--snptable", type = str, action = "store", dest = "snptable", required = True, help = "Path to file containing the .csv output of 25_RFMix_loc_anc.py")

#established part of using GEMMA
parser.add_argument("--relatedness", type = str, action = "store", dest = "relatedness", required = True, help = "Path to file containing relatedness matrix w/o IIDs for only individuals in analysis.")
parser.add_argument("--BIMBAM", type = str, action = "store", dest = "BIMBAM", required = True, help = "Path to file with BIMBAM-formatted genotypes.")
parser.add_argument("--pheno", type = str, action = "store", dest = "pheno", required = True, help = "Path to file containing phenotypic information w/o IIDs for only individuals in analysis.")
parser.add_argument("--covariates", type = str, action = "store", dest = "covariates", required = False, help = "Path to file containing covariates w/o IIDs for only individuals in analysis.")
parser.add_argument("--anno", type = str, action = "store", dest = "anno", required = False, help = "Path to file containing the annotations.")
parser.add_argument("--output", type = str, action = "store", dest = "output", required = False, default = "", help = "Name of output file")
args = parser.parse_args()

print("Reading input files.")
loc_anc_cov = pd.read_csv(args.snptable, delimiter=',', encoding="utf-8-sig")
if args.BIMBAM.endswith(".gz"):
    BIMBAM_chunks = []
    for chunk in pd.read_table(args.BIMBAM, compression='gzip', sep='\t', header = None, index_col = 0, chunksize = 20000):
        BIMBAM_chunks.append(chunk) #when your data too thicc
    BIMBAM = pd.concat(BIMBAM_chunks, axis = 0)
    del BIMBAM_chunks
    del chunk
else:
    BIMBAM_chunks = []
    for chunk in pd.read_table(args.BIMBAM, sep='\t', header = None, index_col = 0, chunksize = 20000):
        BIMBAM_chunks.append(chunk) #when your data too thicc
    BIMBAM = pd.concat(BIMBAM_chunks, axis = 0)
    del BIMBAM_chunks
    del chunk

#following are just to be used in GEMMA input
covariates_file = args.covariates
if args.anno is None:
    anno = " "
else:
    anno = " -a " + args.anno + " "
if args.covariates is None:
    covariates_file = " "
else:
    covariates_file = " -c " + args.covariates + " "
pheno_file = args.pheno
relatedness = args.relatedness
output = args.output

'''
BIMBAM_file = "BIMBAM/test.txt.gz"
covariates_file = "covariates.txt"
anno = " -a anno/test.txt "
pheno_file = "pheno.txt"
relatedness = "relatedness.txt"
output = "test"

SNPs = np.loadtxt("test_SNPs.txt", dtype = 'string')#, engine='python')
loc_anc_cov = pd.read_csv("test.csv", delimiter=',', encoding="utf-8-sig")
BIMBAM_chunks = []
for chunk in pd.read_table(BIMBAM_file, compression='gzip', sep='\t', header = None, index_col = 0, chunksize = 20000):
    BIMBAM_chunks.append(chunk) #when your data too thicc
BIMBAM = pd.concat(BIMBAM_chunks, axis = 0)
del BIMBAM_chunks
del chunk

pheno_num = 1
pheno_name_rank = "CHOL_rank"
ind = inds[0]
'''

print("Formatting input for processing.")
loc_anc_cov['IID'] = loc_anc_cov['IID'].str.replace(r'.*:', '')
inds = loc_anc_cov['IID'].tolist()
loc_anc_cov = loc_anc_cov.set_index('IID').transpose()
SNPs = BIMBAM[[1, 2]]
SNPs = SNPs.reset_index()
SNPs.columns = ['rs', 'A1', 'A0']
SNP_list = SNPs['rs'].tolist()

#format of phenotype file
pheno = range(1,5)
pheno_name = ["CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank"]

print("Creating ancestry-specific dosage files.")

IBS_file = open("BIMBAM/IBS" + output + ".txt", "a+")
NAT_file = open("BIMBAM/NAT" + output + ".txt", "a+")
YRI_file = open("BIMBAM/YRI" + output + ".txt", "a+")
    
progress_landmarks_ind = np.linspace(0, len(inds), 21, dtype = int).tolist()
num_ind = 0
    
'''
Okay so for whoever bothers to read this
I tried making a large data frame and lists of lists with the dosages
However, this made the addition time for each individual raise incrementally
And that's not good for 12,000 people
So, I write each dosage as a line into a file individually so the addition time is constant
And then read that all back in and add the BIMBAM information
So I swear I tried to do it a fancier way but it's so slow with so many people
'''
    
for ind in inds:
    #iterate through cols
    ind = str(ind)
    ind_df = loc_anc_cov[[ind]] 
    ind_df['NAT'], ind_df['IBS'], ind_df['YRI'] = ind_df[ind].str.split('\t', 2).str #split each individual's column into 3
    ind_df = ind_df.drop(ind, axis = 1).transpose().applymap(str)
    #pull from one ancestry each
        #assemble a BIMBAM file except it's local ancestries
    IBS_list = ind_df.loc['IBS'].tolist()
    NAT_list = ind_df.loc['NAT'].tolist()
    YRI_list = ind_df.loc['YRI'].tolist()
       
    IBS_file.write("\t".join(IBS_list) + "\n")
    NAT_file.write("\t".join(NAT_list) + "\n")
    YRI_file.write("\t".join(YRI_list) + "\n")
        
    num_ind = num_ind + 1
    if num_ind in set(progress_landmarks_ind): #print progress by 5% increments
        progress = progress_landmarks_ind.index(num_ind)
        print("Individual ancestry dosage conversion is " + str(progress * 5) + "% complete.")
    
#so now you have dosages for each ancestry
IBS_file.close()
NAT_file.close()
YRI_file.close()
    
#make into BIMBAM
IBS_BIMBAM = pd.read_table("BIMBAM/IBS" + output + ".txt", sep='\t', header = None).transpose()
NAT_BIMBAM = pd.read_table("BIMBAM/NAT" + output + ".txt", sep='\t', header = None).transpose()
YRI_BIMBAM = pd.read_table("BIMBAM/YRI" + output + ".txt", sep='\t', header = None).transpose()

#add SNP info
IBS_BIMBAM = pd.concat([SNPs, IBS_BIMBAM], axis=1)
NAT_BIMBAM = pd.concat([SNPs, NAT_BIMBAM], axis=1)
YRI_BIMBAM = pd.concat([SNPs, YRI_BIMBAM], axis=1)
    
#write to file
IBS_BIMBAM.to_csv("BIMBAM/IBS" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
NAT_BIMBAM.to_csv("BIMBAM/NAT" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
YRI_BIMBAM.to_csv("BIMBAM/YRI" + output + ".txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
    
#phenotype loop    
for pheno_num, pheno_name_rank in zip(pheno, pheno_name):
    print("Starting analyses on " + pheno_name_rank + ".")
    for pop in ['NAT', 'IBS', 'YRI']:
        GEMMA_command = "gemma -g BIMBAM/" + pop + output + ".txt.gz -p " + pheno_file + " -n " + str(pheno_num) + anno + " -k " + relatedness + covariates_file + " -lmm 4 -notsnp -o " + output + "_" + pheno_name_rank + "_" + pop
        os.system(GEMMA_command + " >> GEMMA_log.txt")
        print("Completed with " + pop + " for " + pheno_name_rank + ".")
    print("Ending analyses on " + pheno_name_rank + ".")

print("Removing extra files.")
os.system("rm -f BIMBAM/IBS" + output + ".txt")
os.system("rm -f BIMBAM/NAT" + output + ".txt")
os.system("rm -f BIMBAM/YRI" + output + ".txt")
os.system("rm -f BIMBAM/IBS" + output + ".txt.gz")
os.system("rm -f BIMBAM/NAT" + output + ".txt.gz")
os.system("rm -f BIMBAM/YRI" + output + ".txt.gz")
print("Analyses in all phenotypes is complete. Have a nice day :)!")
