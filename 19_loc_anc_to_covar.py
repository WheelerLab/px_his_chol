#VERY VERY MESSED UP FORMATTING, WORK IN PROGRESS
#"imputes" local ancestry between markers
#probably run once per chromosome?

import pandas as pd
import argparse
#should I make a three dimensional array or is that too ambitious

#unnote when out of testing
#parser = argparse.ArgumentParser()
#parser.add_argument("--loc_anc", type = str, action = "store", dest = "loc_anc", required = True, help = "Path to file containing local ancestry converted from MOSAIC")
#parser.add_argument("--output_prefix", type = str, action = "store", dest = "output_prefix", required = True, help = "Prefix for 'imputed' ancestry file")
#parser.add_argument("--snpfile", type = str, action = "store", dest = "snpfile", required = True, help = "Path to snpfile used in HAPI-UR")
#parser.add_argument("--chr", type = int, action = "store", dest = "chr", required = True, help = "Chromosome under analysis")
#parser.add_argument("--sig_gene_SNPs", type = str, action = "store", dest = "sig_SNP_genes", required = True, help = "List of genes to prune the list of SNPs around (1 Mb before and after)")
#args = parser.parse_args()

print("Reading input file.")
#local_anc = pd.read_csv(args.loc_anc)
#snpfile = pd.read_table(args.snpfile)
#chr = args.chr

#testing files
local_anc = pd.read_csv("loc_anc.csv", dtype={'bp':str}, engine='python')
snpfile = pd.read_table("snpfile.22", delim_whitespace = True, header = None)
chr = 22

snpfile.columns = ['rs', 'chr', 'cM', 'bp', 'A1', 'A2']
snpfile = snpfile[['rs', 'bp']]
hap_list = local_anc['haplotype'].unique() #keep one of each hap
gene_start_end = pd.read_csv("gene_start_end.csv")
gene_start_end_chr = gene_start_end.loc[gene_start_end['chr'] == chr] #subset genes to only relevant chr

#PRUNING TO RELEVENT SNPs
#with open(args.sig_SNP_genes) as f:
#    sig_genes = f.readlines()
        
#FOR TESTING
sig_genes = ["CECR1", "SNAP29", "GNAZ", "TEF", "PRR5"]        

#prune gene_start_end_chr to just sig_genes

for gene_start_end_chr_row in gene_start_end_chr.itertuples():
    keep_gene_start_end = []
    if gene_start_end_chr_row[3] in sig_genes:
        keep_gene_start_end.append(list(gene_start_end_chr_row))
    keep_gene_start_end = pd.DataFrame(keep_gene_start_end)
    if keep_gene_start_end.empty:
        print("There are no significant genes on chr. " + str(chr) + ". Exiting program now.")
        exit()
    keep_gene_start_end.columns = ['index', 'chr', 'gene', 'gene_name', 'start', 'end']
    keep_gene_start_end = keep_gene_start_end.drop('index', axis=1).drop('gene', axis=1).drop('chr', axis=1)
    
    #keep SNPs in hap_SNPs only if they are w/in 1 Mb of start/end of sig genes
    keep_SNP = []
    for keep_gene_start_end_row in keep_gene_start_end.itertuples():
        for snpfile_row in snpfile.itertuples():
            if (snpfile_row[2] > (keep_gene_start_end_row[2] - 1000000)) and (snpfile_row[2] < (keep_gene_start_end_row[3] + 1000000)): #if base pair of SNP is w/in 1 Mb of start or end of gene
                keep_SNP.append(list(snpfile_row))
    keep_SNP = pd.DataFrame(keep_SNP)
    keep_SNP.columns = ['index', 'anc', 'bp', 'haplotype', 'rs']
    keep_SNP = keep_SNP.drop('index', axis=1)
        





#ugh I should've done SNPs as the outer loop shouldn't I
for hap in hap_list:
    hap = hap_list[15] #remove when testing is done
    print("Creating most likely ancestry file for " + hap + ".")
    hap_LoL = []
    
    #get the hap SNPs from the dataframe
    for local_anc_row in local_anc.itertuples():
        if local_anc_row[3] == hap and local_anc_row[1] > 0.9: #keep matching hap and determine which ancestry is most likely at that BP
            hap_LoL.append(list(local_anc_row)) #add rows to new hap-specific dataframe
    hap_df = pd.DataFrame(hap_LoL)
    hap_df.columns = ['index', 'p', 'bp', 'haplotype', 'rs', 'cM', 'anc']
    hap_df = hap_df.drop('index', axis=1).drop('p', axis=1).drop('cM', axis=1)
    
    #merge w/ SNP file
    hap_SNP = pd.concat([hap_df, snpfile])
    hap_SNP['haplotype'] = hap
    
    #order ascendingly
    hap_SNP['bp'] = hap_SNP['bp'].astype(float) #because MOSAIC has decimal base pairs for some reason
    hap_SNP['anc'] = hap_SNP['anc'].astype('category')
    hap_SNP = hap_SNP.sort_values('bp') #wait do I order by base pair or centimorgans?
    
    #impute ancestry for SNPs (this method was easier than I thought)
    hap_SNP['anc'] = hap_SNP['anc'].interpolate(method='pad') #fills NA with the closest value before
        #occasionally "cannot reindex from a duplicate axis"

    #remove non- and duplicated SNPs
    hap_SNP = hap_SNP.dropna()
    
    #OPTIONAL: prune SNPs to those w/in a Mb of sig. genes
    if args.sig_SNP_genes is not None:
        print("Pruning SNPs to those only within a megabase of significant genes.")
        keep_SNPs = []
        
        #with open(args.sig_SNP_genes) as f:
        #    sig_genes = f.readlines()
        
        #FOR TESTING
        sig_genes = ["CECR1", "SNAP29", "GNAZ", "TEF", "PRR5"]        
        
        #prune gene_start_end_chr to just sig_genes
        #find a way to not make this run into an error if there are no sig genes on the chr
        keep_gene_start_end = []
        for gene_start_end_chr_row in gene_start_end_chr.itertuples():
            if gene_start_end_chr_row[3] in sig_genes:
                keep_gene_start_end.append(list(gene_start_end_chr_row))
        keep_gene_start_end = pd.DataFrame(keep_gene_start_end)
        if keep_gene_start_end.empty:
            break #will this skip the args.sig_SNP_genes? I can't google one properly
        keep_gene_start_end.columns = ['index', 'chr', 'gene', 'gene_name', 'start', 'end']
        keep_gene_start_end = keep_gene_start_end.drop('index', axis=1).drop('gene', axis=1).drop('chr', axis=1)
    
        #keep SNPs in hap_SNPs only if they are w/in 1 Mb of start/end of sig genes
        keep_hap_SNP = []
        for keep_gene_start_end_row in keep_gene_start_end.itertuples():
            for hap_SNP_row in hap_SNP.itertuples():
                if (hap_SNP_row[2] > (keep_gene_start_end_row[2] - 1000000)) and (hap_SNP_row[2] < (keep_gene_start_end_row[3] + 1000000)): #if base pair of SNP is w/in 1 Mb of start or end of gene
                    keep_hap_SNP.append(list(hap_SNP_row))
        keep_hap_SNP = pd.DataFrame(keep_hap_SNP)
        keep_hap_SNP.columns = ['index', 'anc', 'bp', 'haplotype', 'rs']
        hap_SNP = keep_hap_SNP.drop('index', axis=1)
        #hap_SNP has now been pruned to only those around sig genes
        
    #NEXT STEPS
    #format: hap, NAT (0 or 1), IBS (0 or 1), YRI (0 or 1)
        #on a per-SNP basis, so probably transpose
    #then combine haps so data is on an individual basis?
