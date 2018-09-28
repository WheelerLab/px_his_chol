#uses local ancestry as a dosage to run an ancestry-by-ancestry level admixture mapping analysis.
#cd /home/angela/px_his_chol/local_anc_GEMMA/MOSAIC_RESULTS/
#python 21_GEMMA_wrapper_admixture_mapping.py --snplist MOSAIC_for_GEMMA_1_snps.txt --snptable MOSAIC_for_GEMMA_1.csv --ind_list MOSAIC_for_GEMMA_1_ind.txt --BIMBAM BIMBAM/chr1.txt.gz --anno anno/anno1.txt --pheno pheno_chr1.txt --relatedness relatedness_chr1.txt --covariate covariates_chr1.txt --output chr1
import argparse
import numpy as np
import pandas as pd
import os
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
#novel-ish part of using GEMMA
parser.add_argument("--snplist", type = str, action = "store", dest = "snplist", required = True, help = "Path to file containing a list of SNPs to be included in analysis (output of 19_loc_anc.py).")
parser.add_argument("--snptable", type = str, action = "store", dest = "snptable", required = True, help = "Path to file containing the .csv output of 19_loc_anc.py")
parser.add_argument("--ind_list", type = str, action = "store", dest = "ind_list", required = True, help = "Path to file containing individuals in the analysis (output of 19_loc_anc.py).")

#established part of using GEMMA
parser.add_argument("--relatedness", type = str, action = "store", dest = "relatedness", required = True, help = "Path to file containing relatedness matrix w/o IIDs for only individuals in analysis.")
parser.add_argument("--BIMBAM", type = str, action = "store", dest = "BIMBAM", required = True, help = "Path to file with BIMBAM-formatted genotypes.")
parser.add_argument("--pheno", type = str, action = "store", dest = "pheno", required = True, help = "Path to file containing phenotypic information w/o IIDs for only individuals in analysis.")
parser.add_argument("--covariates", type = str, action = "store", dest = "covariates", required = False, help = "Path to file containing covariates w/o IIDs for only individuals in analysis.")
parser.add_argument("--anno", type = str, action = "store", dest = "anno", required = False, help = "Path to file containing the annotations.")
parser.add_argument("--output", type = str, action = "store", dest = "output", required = False, help = "Name of output file")
args = parser.parse_args()

print("Reading input files.")
SNPs = np.loadtxt(args.snplist, dtype = 'string')#, engine='python')
loc_anc_cov = pd.read_csv(args.snptable, delimiter=',', encoding="utf-8-sig")
#region = pd.read_table(args.region, delim_whitespace = True, dtype = {'region':object})
if args.BIMBAM.endswith(".gz"):
    BIMBAM = pd.read_table(args.BIMBAM, compression='gzip', sep='\t', header = None, index_col = 0)
else:
    BIMBAM = pd.read_table(args.BIMBAM, header = None, index_col = 0)

#following are just to be used in GEMMA input
BIMBAM_file = args.BIMBAM
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

'''
BIMBAM_file = "BIMBAM/chr1.txt.gz"
covariates_file = "covariates_chr1.txt"
anno = " -a anno/anno1.txt "
pheno_file = "pheno_chr1.txt"
relatedness = "relatedness_chr1.txt"

SNPs = np.loadtxt("MOSAIC_for_GEMMA_1_snps.txt", dtype = 'string')#, engine='python')
loc_anc_cov = pd.read_csv("MOSAIC_for_GEMMA_1.csv", delimiter=',', encoding="utf-8-sig")
BIMBAM = pd.read_table(BIMBAM_file, compression='gzip', sep='\t', header = None, index_col = 0)

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
    
    progress_landmarks_ind = np.linspace(0, len(inds), 21, dtype = int).tolist()
    num_ind = 0
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
        
        num_ind = num_ind + 1
        if num_ind in set(progress_landmarks_ind): #print progress by 5% increments
            progress = progress_landmarks_ind.index(num_ind)
            print("Individual ancestry dosage conversion is " + str(progress * 5) + "% complete.")
    
    #so now you have BIMBAMs for each ancestry
    IBS.to_csv("BIMBAM/IBS.txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
    NAT.to_csv("BIMBAM/NAT.txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
    YRI.to_csv("BIMBAM/YRI.txt.gz", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format='%12f', compression = "gzip")
    
    for pop in ['NAT', 'IBS', 'YRI']:
        if args.output is not None:
            GEMMA_command = "gemma -g BIMBAM/" + pop + ".txt -p " + pheno_file + " -n " + str(pheno_num) + anno + " -k " + relatedness + covariates_file + " -lmm 4 -o " + pheno_name_rank + "_" + pop
            os.system(GEMMA_command + " >> GEMMA_log.txt")
        else:
            GEMMA_command = "gemma -g BIMBAM/" + pop + ".txt -p " + pheno_file + " -n " + str(pheno_num) + anno + " -k " + relatedness + covariates_file + " -lmm 4 -o " + args.output + "_" + pheno_name_rank + "_" + pop
            os.system(GEMMA_command + " >> GEMMA_log.txt")
        
    
    print("Ending analyses on " + pheno_name_rank + ".")
print("Analyses in all phenotypes is complete. Have a nice day :)!")
