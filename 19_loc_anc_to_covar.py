#"imputes" local ancestry between markers to use local ancestry as a covariate in GEMMA on a SNP-by-SNP basis
#probably run once per chromosome?

import argparse
import numpy as np
import pandas as pd

#unnote when out of testing
parser = argparse.ArgumentParser()
parser.add_argument("--loc_anc", type = str, action = "store", dest = "loc_anc", required = True, help = "Path to file containing local ancestry converted from MOSAIC")
parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", required = False, default = "MOSAIC_for_GEMMA", help = "Prefix for 'imputed' ancestry file")
parser.add_argument("--snpfile", type = str, action = "store", dest = "snpfile", required = True, help = "Path to snpfile used in HAPI-UR")
parser.add_argument("--chr", type = int, action = "store", dest = "chr", required = True, help = "Chromosome under analysis")
parser.add_argument("--sig_genes", type = str, action = "store", dest = "sig_genes", required = False, help = "List of genes to prune the list of SNPs around (1 Mb before and after)")
args = parser.parse_args()

print("Reading input files.")
local_anc = pd.read_csv(args.loc_anc, dtype={'bp':float})#, engine='python')
snpfile = pd.read_table(args.snpfile, delim_whitespace = True, header = None)
output_prefix = args.output_prefix
chr = args.chr
if args.sig_genes is None:
    sig_genes = [] #keep all SNPs. SIGNIFICANTLY SLOWER.
    print("No significant gene file called, so program will be keeping all SNPs and will be significantly slowed.")
else:
    sig_genes = np.loadtxt(args.sig_genes, dtype = 'string')

#testing files
#local_anc = pd.read_csv("loc_anc.csv", dtype={'bp':float}, engine='python')
#snpfile = pd.read_table("snpfile.22", delim_whitespace = True, header = None)
#chr = 22
#output_prefix = "test_100_ind"
#sig_genes = ["CECR1", "SNAP29", "GNAZ", "TEF", "PRR5"]  

#local_anc = local_anc.loc[local_anc['prob'] > 0.9] #keep ancestry probabilities > 0.9
local_anc = local_anc.loc[local_anc['p'] > 0.9] #keep ancestry probabilities > 0.9
local_anc = local_anc[['bp', 'haplotype', 'anc']]
snpfile.columns = ['rs', 'chr', 'cM', 'bp', 'A1', 'A2']
snpfile = snpfile[['rs', 'bp']]
hap_list = local_anc['haplotype'].unique() #keep one of each hap

#or input fam file used to make GEMMA and go from there?
ind_list = []
for hap in hap_list:
    ind_list.append(hap[:-2])
ind_list = sorted(set(ind_list), key=ind_list.index) #remove duplicates but preserve order

gene_start_end = pd.read_csv("gene_start_end.csv")
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
    print("Kept all SNPs from an original " + str(len(snpfile)) + " SNPs from the snpfile.") 
else:
    for SNP in iter(test_SNPs_snpfile):
        for keep_gene_start_end_row in keep_gene_start_end.itertuples():
            if (SNP > (keep_gene_start_end_row[2] - 1000000)) and (SNP < (keep_gene_start_end_row[3] + 1000000)): #if base pair of known SNP is w/in 1 Mb of start or end of gene
                keep_SNP.append(SNP)
    keep_SNP = pd.DataFrame(keep_SNP)
    keep_SNP.columns = ['bp']
    keep_SNP = keep_SNP.merge(snpfile, on = 'bp', how = 'left')
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
print("Starting SNP ancestry imputation.")
imputed_haplotypes = pd.DataFrame(columns=['anc', 'bp', 'haplotype', 'rs'])
for hap in hap_list:
    ind_haplotype = []
    for keep_local_anc_row in keep_local_anc.itertuples():
        if hap == keep_local_anc_row[2]:
            ind_haplotype.append(keep_local_anc_row)
    hap_df = pd.DataFrame(ind_haplotype)
    if hap_df.empty: #remove from ind_list to prevent further issues
        print(hap + " is empty. Skipping haplotype and removing individual from further analyses.")
        ind = hap[:-2]
        ind_list.remove(ind)        
        continue
    hap_df = hap_df.drop('Index', axis = 1)
    
    #impute ancestry for SNPs
    hap_SNP = pd.concat([hap_df, keep_SNP]).sort_values('bp')
    hap_SNP['haplotype'] = hap
    hap_SNP['anc'] = hap_SNP['anc'].ffill().bfill() #fills NA with the closest value before and after
        #try to keep it NA in between flips? but I don't know of a method that does that
            #unless I write one I guess

    #remove non-SNPs
    hap_SNP = hap_SNP.dropna(how = 'any', axis = 0)
    imputed_haplotypes = imputed_haplotypes.append(hap_SNP)
print("Haplotype ancestry imputation complete.")

#make covariate file SNP by SNP
print("Starting SNP ancestry covariate file.")
study_SNPs = pd.DataFrame(imputed_haplotypes['rs'].unique())
study_SNPs.columns = ['rs']

#collapse haplotype into one individual
for ind in ind_list:
    ind_anc_A = [] #haplotype 1
    ind_anc_B = [] #haplotype 2
    for imputed_haplotypes_row in imputed_haplotypes.itertuples():
        if imputed_haplotypes_row[3] == (ind + "_A"):
            ind_anc_A.append(imputed_haplotypes_row)
        elif imputed_haplotypes_row[3] == (ind + "_B"):
            ind_anc_B.append(imputed_haplotypes_row)
    ind_anc_A = pd.DataFrame(ind_anc_A)
    ind_anc_A = ind_anc_A.drop('Index', axis = 1).drop('bp', axis = 1).drop('haplotype', axis = 1)
    ind_anc_A.columns = ['anc_A', 'rs']
    ind_anc_B = pd.DataFrame(ind_anc_B)
    ind_anc_B = ind_anc_B.drop('Index', axis = 1).drop('bp', axis = 1).drop('haplotype', axis = 1)
    ind_anc_B.columns = ['anc_B', 'rs']
    ind_anc = ind_anc_A.set_index('rs').join(ind_anc_B.set_index('rs'))
    
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
    study_SNPs = study_SNPs.set_index('rs').join(anc_dosage_df.set_index('rs'))
    
    #restore SNPs column
    study_SNPs.index.name = 'rs'
    study_SNPs.reset_index(inplace = True)

#write list of SNPs to use in GEMMA (-snps)
study_SNPs.rs.to_csv(output_prefix + "_" + str(chr) + "_snps.txt", index = False, header = False)

#write file
study_SNPs_t = study_SNPs.transpose()
study_SNPs_t.index.name = 'rs'
study_SNPs_t.reset_index(inplace = True)
study_SNPs_t = study_SNPs_t.set_value(0, 'rs', 'IID')
study_SNPs_t.to_csv(output_prefix + "_" + str(chr) + ".csv", sep = ",", na_rep = "NA\tNA\tNA", index = False, header = False)
print("Completed writing SNP and SNP ancestry covariate file to " + output_prefix + "_" + str(chr) + "_snps.txt and " + output_prefix + "_" + str(chr) + ".csv. Have a nice day!")

#from here, parse on a SNP-by-SNP basis for each GEMMA run (see 20_GEMMA_wrapper.py)
    
