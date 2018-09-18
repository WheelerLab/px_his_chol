#wrapper for admixture mapping in GEMMA
import argparse
import numpy as np
import pandas as pd
import os
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
#novel-ish part of using GEMMA
parser.add_argument("--snplist", type = str, action = "store", dest = "snplist", required = True, help = "Path to file containing a list of SNPs to be included in analysis.")
parser.add_argument("--snptable", type = str, action = "store", dest = "snptable", required = True, help = "Path to file containing the .csv output of 19_loc_anc_to_covar.py")
parser.add_argument("--region", type = str, action = "store", dest = "region", required = True, help = "Path to file containing region information.")

#established part of using GEMMA
parser.add_argument("--BIMBAM", type = str, action = "store", dest = "BIMBAM", required = True, help = "Path to file with BIMBAM-formatted genotypes.")
parser.add_argument("--anno", type = str, action = "store", dest = "anno", required = False, help = "Path to file containing the annotations.")
parser.add_argument("--pheno", type = str, action = "store", dest = "pheno", required = True, help = "Path to file containing phenotypic information.")
parser.add_argument("--covariates", type = str, action = "store", dest = "covariates", required = True, help = "Path to file containing covariates.")
parser.add_argument("--relatedness", type = str, action = "store", dest = "relatedness", required = True, help = "Path to file containing relatedness matrix w/o IIDs.")
parser.add_argument("--output", type = str, action = "store", dest = "output", required = False, help = "Name of output file")
args = parser.parse_args()

print("Reading input files.")
SNPs = np.loadtxt(args.snplist, dtype = 'string')#, engine='python')
loc_anc_cov = pd.read_csv(args.snptable, delimiter=',', encoding="utf-8-sig")
region = pd.read_table(args.region, delim_whitespace = True, dtype = {'region':object})
BIMBAM = pd.read_table(args.BIMBAM, header = None, index_col = 0)

#following are just to be used in GEMMA input
BIMBAM_file = args.BIMBAM
covariates_file = args.covariates
if args.anno is None:
    anno = " "
else:
    anno = " -a " + args.anno + " "
pheno_file = args.pheno
relatedness = args.relatedness

'''
#testing files
BIMBAM_file = "BIMBAM/chr22.txt"
anno =  " -a anno/anno22.txt "
pheno_file = "pheno_woIID.txt"
relatedness = "relatedness_matrix_woIID_noNeg.txt"
covariates_file = "cov_region_5_PCs_100_ind.txt"

print("Reading in input files.")
SNPs = list(np.loadtxt("MOSAIC_for_GEMMA_22_snps.txt", dtype = 'string'))#, engine='python')
loc_anc_cov = pd.read_csv("MOSAIC_for_GEMMA_22.csv", delimiter=',', encoding="utf-8-sig")
region = pd.read_table("region.txt", delim_whitespace = True, dtype = {'region':object})
BIMBAM = pd.read_table(BIMBAM_file, header = None)
pheno_num = 1
pheno_name_rank = "CHOL_rank"
pop = "NAT"
ind = "SoL100865"
'''

print("Formatting input for processing.")
BIMBAM = pd.read_table(BIMBAM_file, header = None)
loc_anc_cov['IID'] = loc_anc_cov['IID'].str.replace(r'.*:', '')
inds = loc_anc_cov['IID'].tolist()
loc_anc_cov = loc_anc_cov.set_index('IID').transpose()
SNPs = BIMBAM[[0, 1, 2]]
SNPs.columns = ['rs', 'A1', 'A0']

#format of phenotype file
pheno = range(1,5)
pheno_name = ["CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank"]

#phenotype loop
for pheno_num, pheno_name_rank in zip(pheno, pheno_name):
    print("Starting analyses on " + pheno_name_rank + ".")
    
    IBS = SNPs
    NAT = SNPs
    YRI = SNPs
    #start individual loop
    for ind in inds:
        #iterate through cols
        ind_df = loc_anc_cov[[ind]]
        ind_df['NAT'], ind_df['IBS'], ind_df['YRI'] = ind_df[ind].str.split('\t', 2).str
        ind_df = ind_df.drop(ind, axis = 1)
        #pull from one ancestry each
            #assemble a BIMBAM file except it's local ancestries
        IBS = IBS.set_index('rs').join(ind_df[['IBS']])
        IBS = IBS.rename(columns = {'IBS':ind})
        IBS.index.name = 'rs'
        IBS.reset_index(inplace = True)
    
        NAT = NAT.set_index('rs').join(ind_df[['NAT']])
        NAT = NAT.rename(columns = {'NAT':ind})
        NAT.index.name = 'rs'
        NAT.reset_index(inplace = True)
    
        YRI = YRI.set_index('rs').join(ind_df[['YRI']])
        YRI = YRI.rename(columns = {'YRI':ind})
        YRI.index.name = 'rs'
        YRI.reset_index(inplace = True)
    
    #so now you have BIMBAMs for each ancestry
    IBS.to_csv("BIMBAM/IBS.txt", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f')
    NAT.to_csv("BIMBAM/NAT.txt", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f')
    YRI.to_csv("BIMBAM/YRI.txt", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f')
    
    for pop in ['NAT', 'IBS', 'YRI']:
        GEMMA_command = "gemma -g BIMBAM/" + pop + ".txt -p " + pheno_file + " -n " + str(pheno_num) + anno + " -k " + relatedness + " -c " + covariates_file + " -lmm 4 -o " + pheno_name_rank + "_output_" + pop
        os.system(GEMMA_command + " >> GEMMA_log.txt")
    
    print("Ending analyses on " + pheno_name_rank + ".")
    
print("Analyses in all phenotypes is complete. Have a nice day :)!")
