#wrapper for SNP-by-SNP level GEMMA
import argparse
import gzip
import numpy as np
import pandas as pd
import subprocess
pd.options.mode.chained_assignment = None

mount = ""
#mount = "/home/angela"

'''
parser = argparse.ArgumentParser()
#novel-ish part of using GEMMA
parser.add_argument("--snplist", type = str, action = "store", dest = "snplist", required = True, help = "Path to file containing a list of SNPs to be included in analysis.")
parser.add_argument("--snptable", type = str, action = "store", dest = "snptable", required = True, help = "Path to file containing the .csv output of 19_loc_anc_to_covar.py")
parser.add_argument("--region", type = str, action = "store", dest = "region", required = True, help = "Path to file containing region information.")
#established part of using GEMMA
parser.add_argument("--BIMBAM", type = str, action = "store", dest = "BIMBAM", required = True, help = "Path to file with BIMBAM-formatted genotypes.")
parser.add_argument("--anno", type = str, action = "store", dest = "anno", required = True, help = "Path to file containing the annotations.")
parser.add_argument("--pheno", type = str, action = "store", dest = "pheno", required = True, help = "Path to file containing phenotypic information.")
parser.add_argument("--relatedness", type = str, action = "store", dest = "relatedness", required = True, help = "Path to file containing relatedness matrix w/o IIDs.")
parser.add_argument("--output", type = str, action = "store", dest = "output", required = False, help = "Name of output file")
args = parser.parse_args()
print("Reading input files.")
SNPs = np.loadtxt(args.snplist, dtype = 'string')#, engine='python')
loc_anc_cov = pd.read_csv(args.snptable, delimiter=',', encoding="utf-8-sig")
region = pd.read_table(args.region, delim_whitespace = True, dtype = {'region':object})
#following are just to be used in GEMMA input
BIMBAM = args.BIMBAM
anno = args.anno
pheno = args.pheno
relatedness = args.relatedness
'''

BIMBAM_file = "BIMBAM/chr22.txt"
anno = "anno/anno22.txt"
pheno_file = "pheno_woIID.txt"
relatedness = "relatedness_matrix_woIID_noNeg.txt"

#testing files
SNPs = list(np.loadtxt("MOSAIC_for_GEMMA_22_snps.txt", dtype = 'string'))#, engine='python')
loc_anc_cov = pd.read_csv("MOSAIC_for_GEMMA_22.csv", delimiter=',', encoding="utf-8-sig")
region = pd.read_table("region.txt", delim_whitespace = True, dtype = {'region':object})
BIMBAM = pd.read_table(BIMBAM_file, header = None)
BIMBAM = BIMBAM.set_index(0)

#remove FID from loc_anc_cov
loc_anc_cov['IID'] = loc_anc_cov['IID'].str.replace(r'.*:', '')

#clean up region and add intercept
region.columns = ['IID', 'region']
region['intercept'] = 1
#region[]
region = region[['IID', 'intercept', 'region']]

#merge all covariates to be selectively pulled later
covariates = region.set_index('IID').join(loc_anc_cov.set_index('IID'))

#format of phenotype file
pheno = range(1,5)
pheno_name = ["CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank"]

#phenotype loop
for pheno_num, pheno_name_rank in zip(pheno, pheno_name):
    #make file for new GEMMA files to append to
    
    #testing
    pheno_num = pheno[0]
    pheno_name_rank = pheno_name[0]
    
    pheno_results = pd.DataFrame(columns=['chr', 'rs', 'ps', 'n_miss', 'allele1', 'allele0', 'af', 'beta', 'se', 'l_remle', 'l_mle', 'p_wald', 'p_lrt', 'p_score', 'anc'])
    
    #start SNP loop
    for SNP in SNPs:
        #testing
        SNP = SNPs[0]
        
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
            #try keeping one column
            single_pop_cov = SNP_cov[['intercept', 'region', pop]]

            #tmp files
            SNP_name = open("tmp_SNP.txt", "w")
            SNP_name.write(SNP) #take argument -SNPs to specify SNPs to study
            SNP_name.close()
        
            #change float when necessary
            single_pop_cov.to_csv("tmp_SNP_cov.txt", sep = "\t", na_rep = "NA", header = False, index = False, quoting = 3, float_format = '%.12g')
        
            #-snp doesn't seem to be working. Pull SNP from BIMBAM file?
            BIMBAM_to_write = BIMBAM.loc[[SNP], :]
            BIMBAM_to_write.to_csv("tmp_BIMBAM.txt", sep = "\t", na_rep = "NA", header = False, index = True, float_format = '%.12g', quoting = 3)
        
            #run GEMMA
            GEMMA_command = "/usr/local/bin/gemma -g tmp_BIMBAM.txt -p " + pheno_file + " -n " + str(pheno_num) + " -a " + anno + " -k " + relatedness + " -c tmp_SNP_cov.txt -lmm 1 -notsnp -o tmp_output_" + pop
        
            #learn how to use this
            subprocess.call([
                'gemma', '-g', BIMBAM_file, '-p', pheno_file, '-n', str(pheno_num), '-a', anno, '-k', relatedness, '-c', 'tmp_SNP_cov.txt', '-snps', 'tmp_SNP.txt', '-lmm', str(4), '-o', 'tmp_output'
            ])
    
            #add to output file
            GEMMA_output = pd.read_table("output/tmp_output_" + pop + ".assoc.txt", delim_whitespace = True)
            GEMMA_output['anc'] = pop
            pheno_results = pheno_results.append(GEMMA_output)

    #write to results
    pheno_results.to_csv(pheno_name_rank + "_results.txt", sep = "\t", na_rep = "NA", header = True, index = False, quoting = 3)
