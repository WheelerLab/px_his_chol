#"imputes" local ancestry between markers to use local ancestry as a covariate in GEMMA on a SNP-by-SNP basis
#example input - python 19_loc_anc_to_covar.py --loc_anc local_anc_22_HCHS.csv --snpfile /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/snpfile.22

import argparse
import numpy as np
import pandas as pd

#unnote when out of testing
parser = argparse.ArgumentParser()
parser.add_argument("--loc_anc", type = str, action = "store", dest = "loc_anc", required = True, help = "Path to file containing local ancestry converted from MOSAIC (output of 18_convert_MOSAIC_output.R)")
parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", required = False, default = "MOSAIC_for_GEMMA", help = "Prefix for 'imputed' ancestry file")
parser.add_argument("--snpfile", type = str, action = "store", dest = "snpfile", required = True, help = "Path to snpfile used in HAPI-UR")
parser.add_argument("--sig_genes", type = str, action = "store", dest = "sig_genes", required = False, help = "List of genes to prune the list of SNPs around (1 Mb before and after)")
args = parser.parse_args()

print("Reading input files.")
local_anc = pd.read_csv(args.loc_anc, dtype={'bp':float})#, engine='python')
snpfile = pd.read_table(args.snpfile, delim_whitespace = True, header = None)
output_prefix = args.output_prefix
chr = int(args.snpfile[-2:].replace(".", ""))
if args.sig_genes is None:
    sig_genes = [] #keep all SNPs. SIGNIFICANTLY SLOWER.
    print("No significant gene file called, so program will be keeping all SNPs and will be significantly slowed.")
else:
    sig_genes = np.loadtxt(args.sig_genes, dtype = 'string')
gene_start_end = pd.read_csv("/home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/gene_start_end.csv")


#testing files
#in /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/100_ind/
'''
local_anc = pd.read_csv("local_anc_1_HCHS_100000.csv", dtype={'bp':float}, engine='python')
snpfile = pd.read_table("snpfile.1", delim_whitespace = True, header = None)
gene_start_end = pd.read_csv("gene_start_end.csv")
output_prefix = "TEST_"
chr = 1
sig_genes = ["PABPC4", "TRIT1"]
'''

local_anc = local_anc[['bp', 'haplotype', 'anc']]
snpfile.columns = ['rs', 'chr', 'cM', 'bp', 'A1', 'A2']
snpfile = snpfile[['rs', 'bp']]
hap_list = local_anc['haplotype'].unique() #keep one of each hap

#or input fam file used to make GEMMA and go from there?
ind_list = []
for hap in hap_list:
    ind_list.append(hap[:-2])
ind_list = sorted(set(ind_list), key = ind_list.index) #remove duplicates but preserve order

gene_start_end_chr = gene_start_end.loc[gene_start_end['chr'] == chr] #subset genes to only relevant chr

#prune gene_start_end_chr to just sig_genes
keep_gene_start_end = []
if len(sig_genes) == 0:
    for gene_start_end_chr_row in gene_start_end_chr.itertuples():
        keep_gene_start_end.append(list(gene_start_end_chr_row))
        keep_gene_start_end = pd.DataFrame(keep_gene_start_end)
    print("Kept all genes on chr. " + str(chr) + ".")
    keep_gene_start_end.columns = ['index', 'chr', 'gene', 'gene_name', 'start', 'end']
    keep_gene_start_end = keep_gene_start_end.drop('index', axis=1).drop('gene', axis=1).drop('chr', axis=1)
else:
    for gene_start_end_chr_row in gene_start_end_chr.itertuples():
        if gene_start_end_chr_row[3] in sig_genes:
            keep_gene_start_end.append(list(gene_start_end_chr_row))
    keep_gene_start_end = pd.DataFrame(keep_gene_start_end)
    print("Kept the " + str(len(keep_gene_start_end.index)) + " significant genes on chr. " + str(chr) + ".")
    keep_gene_start_end.columns = ['index', 'chr', 'gene', 'gene_name', 'start', 'end']
    keep_gene_start_end = keep_gene_start_end.drop('index', axis=1).drop('gene', axis=1).drop('chr', axis=1)
    
#keep SNPs in snpfile only if they are w/in 1 Mb of start/end of sig genes
keep_SNP = []
test_SNPs_snpfile = set(snpfile['bp'].tolist()) #is this faster than having nested loops for a dataframe?
if len(sig_genes) == 0:
    for SNP in iter(test_SNPs_snpfile):
        for keep_gene_start_end_row in keep_gene_start_end.itertuples():
            keep_SNP.append(SNP)
    keep_SNP = pd.DataFrame(keep_SNP)
    keep_SNP.columns = ['bp']
    keep_SNP = keep_SNP.merge(snpfile, on = 'bp', how = 'left')
    keep_SNP = keep_SNP.sort_values('bp').drop_duplicates()
    print("Kept all SNPs from an original " + str(len(snpfile)) + " SNPs from the snpfile.") 
else:
    for SNP in iter(test_SNPs_snpfile):
        for keep_gene_start_end_row in keep_gene_start_end.itertuples():
            if (SNP > (keep_gene_start_end_row[2] - 1000000)) and (SNP < (keep_gene_start_end_row[3] + 1000000)): #if base pair of known SNP is w/in 1 Mb of start or end of gene
                keep_SNP.append(SNP)
    keep_SNP = pd.DataFrame(keep_SNP)
    keep_SNP.columns = ['bp']
    keep_SNP = keep_SNP.merge(snpfile, on = 'bp', how = 'left')
    keep_SNP = keep_SNP.sort_values('bp').drop_duplicates()
    print("Kept " + str(len(keep_SNP)) + " SNPs from an original " + str(len(snpfile)) + " SNPs from the snpfile.") 

#prune local_anc to relevant SNPs (this wasn't very useful in the example data)
#make local anc row bps into a set for speed
keep_local_anc = []
test_local_anc_SNPs = set(local_anc['bp'].tolist()) #for speed purposes only
if len(sig_genes) == 0:
    for SNP in iter(test_local_anc_SNPs):
        for keep_gene_start_end_row in keep_gene_start_end.itertuples():
            keep_local_anc.append(SNP)
else:
    for SNP in iter(test_local_anc_SNPs):
        for keep_gene_start_end_row in keep_gene_start_end.itertuples():
            if (SNP > (keep_gene_start_end_row[2] - 1000000)) and (SNP < (keep_gene_start_end_row[3] + 1000000)): #if base pair of MOSAIC SNP is w/in 1 Mb of start or end of gene
                keep_local_anc.append(SNP)
keep_local_anc = pd.DataFrame(keep_local_anc)
keep_local_anc.columns = ['bp']
keep_local_anc = keep_local_anc.merge(local_anc, on = 'bp', how = 'left')
keep_local_anc['bp'] = keep_local_anc['bp'].astype(float)
keep_local_anc['anc'] = keep_local_anc['anc'].astype('category')

#impute SNPs using .ffill()
#RESTRUCTURE SO YOU JUST APPEND INSTEAD OF MAKING A LARGE DATA FRAME
print("Starting SNP ancestry imputation and dosage file.")
anc_dosage_write = open(output_prefix + "_" + str(chr) + ".csv", "a+")
anc_dosage_write.write("IID," + ",".join(keep_SNP['rs'].tolist()) + "\n")
keep_SNP.rs.to_csv(output_prefix + "_" + str(chr) + "_snps.txt", index = False, header = False)
progress_landmarks_ind = np.linspace(0, len(ind_list), 21, dtype = int).tolist()
num_ind = 0

for ind in ind_list: #what part in here takes so long?
    imputed_haplotypes = pd.DataFrame(columns=['anc', 'bp', 'haplotype', 'rs'])
    #FIRST HAPLOTYPE
    ind_haplotype_A = []
    hap_A = ind + "_A"
    ind_haplotype_A = keep_local_anc.loc[keep_local_anc['haplotype'] == hap_A]
    if ind_haplotype_A.empty: #remove from ind_list to prevent further issues
        print(hap_A + " did not have sufficient coverage. Skipping haplotype and removing individual from further analyses.")
        ind_list.remove(ind) 
        num_ind = num_ind + 1
        continue

    #impute ancestry for SNPs
    hap_SNP_A = pd.concat([ind_haplotype_A, keep_SNP]).sort_values('bp')
    hap_SNP_A['haplotype'] = hap_A
    hap_SNP_A['anc'] = hap_SNP_A['anc'].ffill().bfill() #fills NA with the closest value before and after
        #try to keep it NA in between flips? but I don't know of a method that does that
            #unless I write one I guess

    #remove non-SNPs
    hap_SNP_A = hap_SNP_A.dropna(how = 'any', axis = 0)
    imputed_haplotypes = imputed_haplotypes.append(hap_SNP_A)
    
    #SECOND HAPLOTYPE
    ind_haplotype_B = []
    hap_B = ind + "_B"
    ind_haplotype_B = keep_local_anc.loc[keep_local_anc['haplotype'] == hap_B]
    if ind_haplotype_B.empty: #remove from ind_list to prevent further issues
        print(hap_B + " did not have sufficient coverage. Skipping haplotype and removing individual from further analyses.")
        ind_list.remove(ind)  
        num_ind = num_ind + 1
        continue

    #impute ancestry for SNPs
    hap_SNP_B = pd.concat([ind_haplotype_B, keep_SNP]).sort_values('bp')
    hap_SNP_B['haplotype'] = hap_B
    hap_SNP_B['anc'] = hap_SNP_B['anc'].ffill().bfill() #fills NA with the closest value before and after
        #try to keep it NA in between flips? but I don't know of a method that does that
            #unless I write one I guess

    #remove non-SNPs
    hap_SNP_B = hap_SNP_B.dropna(how = 'any', axis = 0)
    imputed_haplotypes = imputed_haplotypes.append(hap_SNP_B)
    ind_anc = hap_SNP_A.set_index('rs').join(hap_SNP_B.set_index('rs'), lsuffix = "_A", rsuffix = "_B")
    ind_anc = ind_anc[["anc_A", "anc_B"]]
    
    #now set ancestry "dosages"
    #first col is NAT, second is IBS, and third is YRI
    anc_dosage = []
    for ind_anc_row in ind_anc.itertuples():
        if ind_anc_row[1] == "NAT" and ind_anc_row[2] == "NAT":
            anc_dosage.append([ind_anc_row[0], "2.0\t0.0\t0.0"]) #both NAT
        elif ind_anc_row[1] == "IBS" and ind_anc_row[2] == "IBS":
            anc_dosage.append([ind_anc_row[0], "0.0\t2.0\t0.0"]) #both IBS
        elif ind_anc_row[1] == "YRI" and ind_anc_row[2] == "YRI":
            anc_dosage.append([ind_anc_row[0], "0.0\t0.0\t2.0"]) #both YRI
        elif (ind_anc_row[1] == "NAT" and ind_anc_row[2] == "IBS") or (ind_anc_row[1] == "IBS" and ind_anc_row[2] == "NAT"):
            anc_dosage.append([ind_anc_row[0], "1.0\t1.0\t0.0"]) #one NAT, one IBS
        elif (ind_anc_row[1] == "NAT" and ind_anc_row[2] == "YRI") or (ind_anc_row[1] == "YRI" and ind_anc_row[2] == "NAT"): 
            anc_dosage.append([ind_anc_row[0], "1.0\t0.0\t1.0"]) #one NAT, one YRI
        elif (ind_anc_row[1] == "YRI" and ind_anc_row[2] == "IBS") or (ind_anc_row[1] == "IBS" and ind_anc_row[2] == "YRI"):
            anc_dosage.append([ind_anc_row[0], "0.0\t1.0\t1.0"]) #one IBS, one YRI  
        else:
            anc_dosage.append([ind_anc_row[0], "NA\tNA\tNA"]) #who knows
    anc_dosage_df = pd.DataFrame(anc_dosage)
    anc_dosage_df.columns = ["rs", ind]
    anc_dosage_df = anc_dosage_df.drop_duplicates()
    anc_dosage_list = anc_dosage_df[ind].tolist()
    
    anc_dosage_write.write(ind + "," + ",".join(anc_dosage_list) + "\n")
    
    num_ind = num_ind + 1
    if num_ind in set(progress_landmarks_ind): #print progress by 5% increments
      progress = progress_landmarks_ind.index(num_ind)
      print("SNP ancestry covariate conversion is " + str(progress * 5) + "% complete.")

#write list of SNPs to use in GEMMA (-snps)
anc_dosage_write.close()
print("Completed writing SNP and SNP ancestry covariate file to " + output_prefix + "_" + str(chr) + "_snps.txt and " + output_prefix + "_" + str(chr) + ".csv. Have a nice day!")
   
