#wrapper for admixture mapping analysis in GEMMA
#NOT UPDATED TO ADMIXTURE MAPPING AS OF 9/17/18
import argparse
import gzip
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
if args.anno is None:
    anno = ""
else:
    anno = " -a " + args.anno + " "
pheno_file = args.pheno
relatedness = args.relatedness

'''
#testing files
BIMBAM_file = "BIMBAM/chr22.txt"
anno = "anno/anno22.txt"
pheno_file = "pheno_woIID.txt"
relatedness = "relatedness_matrix_woIID_noNeg.txt"

print("Reading in input files.")
SNPs = list(np.loadtxt("MOSAIC_for_GEMMA_22_snps.txt", dtype = 'string'))#, engine='python')
loc_anc_cov = pd.read_csv("MOSAIC_for_GEMMA_22.csv", delimiter=',', encoding="utf-8-sig")
region = pd.read_table("region.txt", delim_whitespace = True, dtype = {'region':object})
BIMBAM = pd.read_table(BIMBAM_file, header = None)
'''

print("Formatting input for processing.")
#remove FID from loc_anc_cov
loc_anc_cov['IID'] = loc_anc_cov['IID'].str.replace(r'.*:', '')

#clean up region and add intercept
region.columns = ['IID', 'region']
region['intercept'] = 1
region = region[['IID', 'intercept', 'region']]

#merge all covariates to be selectively pulled later
covariates = region.set_index('IID').join(loc_anc_cov.set_index('IID'))

#for progress updates
progress_landmarks = np.linspace(0, len(loc_anc_cov.columns), 11, dtype = int).tolist()

#format of phenotype file
pheno = range(1,5)
pheno_name = ["CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank"]

#phenotype loop
for pheno_num, pheno_name_rank in zip(pheno, pheno_name):
    print("Starting analyses on " + pheno_name_rank + ".")
    pheno_results = pd.DataFrame(columns=['chr', 'rs', 'ps', 'n_miss', 'allele1', 'allele0', 'af', 'beta', 'se', 'l_remle', 'l_mle', 'p_wald', 'p_lrt', 'p_score', 'anc'])
    
    #start SNP loop
    SNP_num = 0
    for SNP in SNPs:
        #rebuild covariate file for each SNP
            #should be the same length as the number of people
        #pull column of imputations from 19_loc_anc_to_covar
        SNP_cov = covariates[['intercept', 'region', SNP]]
        
        #estimates need to be split b/c .to_csv doesn't like when they're not
        SNP_cov['NAT'], SNP_cov['IBS'], SNP_cov['YRI'] = SNP_cov[SNP].str.split('\t', 2).str
        SNP_cov = SNP_cov.drop(SNP, axis = 1)
        
        SNP_cov['NAT'] = SNP_cov['NAT'].astype(float)
        SNP_cov['IBS'] = SNP_cov['IBS'].astype(float)
        SNP_cov['YRI'] = SNP_cov['YRI'].astype(float)
        #intercept region NAT IBS YRI

        for pop in ['NAT', 'IBS', 'YRI']: #because it doesn't work all at once
            #run one ancestry at a time b/c they're collinear
            single_pop_cov = SNP_cov[['intercept', 'region', pop]]

            #tmp files
            SNP_name = open("tmp_SNP.txt", "w")
            SNP_name.write(SNP) #take argument -SNPs to specify SNPs to study
            SNP_name.close()
        
            #change float when necessary
            single_pop_cov.to_csv("tmp_SNP_cov.txt", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format = '%.12g')
        
            #run GEMMA
            GEMMA_command = "gemma -g " + BIMBAM_file + " -p " + pheno_file + " -n " + str(pheno_num) + anno + " -k " + relatedness + " -snps tmp_SNP.txt -c tmp_SNP_cov.txt -lmm 4 -notsnp -o tmp_output_" + pop
            os.system(GEMMA_command + " >> GEMMA_log.txt")
        
            #add to output file
            GEMMA_output = pd.read_table("output/tmp_output_" + pop + ".assoc.txt", delim_whitespace = True)
            GEMMA_output['anc'] = pop
            pheno_results = pheno_results.append(GEMMA_output)
            
        SNP_num = SNP_num + 1
        if SNP_num in set(progress_landmarks): #print progress by 10% increments
            progress = progress_landmarks.index(SNP_num)
            print("Analysis in " + pheno_name_rank + " is " + str(progress * 10) + "% complete.")
        
    #write to results
    print("Ending analyses on " + pheno_name_rank + ".")
    pheno_results.to_csv(pheno_name_rank + "_results.txt", sep = "\t", na_rep = "NA", header = True, index = False, quoting = 3, float_format='%g')

print("Analyses in all phenotypes is complete. Have a nice day :)!")
