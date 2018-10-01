#Takes information from PrediXcan-style dosages to make SNP annotation and BIMBAM files for GEMMA
#cd /home/angela/px_his_chol/local_anc_GEMMA/MOSAIC_RESULTS/
#python 02_PrediXcan_dosages_to_GEMMA.py --dosage_path /home/angela/px_his_chol/Imputation/UMich/UMich_dosages/ --local_anc_samples MOSAIC_for_GEMMA_1_ind.txt --local_anc_SNPs MOSAIC_for_GEMMA_1_snps.txt --chr 1
import argparse
import gzip
import numpy as np
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("--dosage_path", type = str, action = "store", dest = "dosage_path", required = True, help = "Path to folder containing dosages and samples.txt")
#parser.add_argument("--BIMBAM_path", type = str, action = "store", dest = "BIMBAM_path", required = False, default = "BIMBAM/", help = "Name of output folder for BIMBAMs.")
#parser.add_argument("--anno_path", type = str, action = "store", dest = "anno_path", required = False, default = "anno/", help = "Name of output folder for annotations.")
parser.add_argument("--local_anc_samples", type = str, action = "store", dest = "local_anc_samples", required = False, help = "Path to file containing samples to include (output of 19_loc_anc.py).")
parser.add_argument("--local_anc_SNPs", type = str, action = "store", dest = "local_anc_SNPs", required = False, help = "Path to file containing SNPs to include (output of 19_loc_anc.py).")
parser.add_argument("--chr", type = str, action = "store", dest = "chr", required = False, help = "Path to chromosome to analyze. If no input, analyzes all 22 pairs of chromosomes.")
args = parser.parse_args()

print("Reading input files.")
dosage_path = args.dosage_path
BIMBAM_path = "BIMBAM/"
anno_path = "anno/"

os.system("mkdir -p BIMBAM")
os.system("mkdir -p anno")

if args.local_anc_samples is not None:
    pd.read_csv(args.local_anc_samples, sep = ":", header = None)
if args.chr is None:
    chrs_to_test = range(1, 23)
if args.chr is not None:
    chrs_to_test = range(int(args.chr), (int(args.chr) + 1))
else:
    chrs_to_test = range(1, 23)
if args.local_anc_SNPs is not None:
    local_anc_SNPs = set(np.loadtxt(args.local_anc_SNPs, dtype = 'string'))

'''
mount = ""
dosage_path = mount + "/px_his_chol/Imputation/UMich/UMich_dosages/"
BIMBAM_path = mount + "/px_his_chol/local_anc_GEMMA/MOSAIC_RESULTS/BIMBAM/"
anno_path = mount + "/px_his_chol/local_anc_GEMMA/MOSAIC_RESULTS/anno/"
local_anc_samples = pd.read_csv(mount + "/px_his_chol/local_anc_GEMMA/MOSAIC_RESULTS/MOSAIC_for_GEMMA_1_ind.txt", sep = ":", header = None)
i = 1
'''

#PX_SNPs = set(numpy.loadtxt("/home/angela/px_his_chol/SNPs_in_PrediXcan_models.txt", dtype = 'string'))
dosage_samples = pd.read_csv(dosage_path + "samples.txt", sep = " ", header = None)
dosage_samples = dosage_samples[[1]]

if args.local_anc_samples is not None:
    local_anc_samples = pd.read_csv(args.local_anc_samples, sep = ":", header = None)
    local_anc_samples = local_anc_samples[[1]]
    local_anc_samples = local_anc_samples.set_index(1)

print("Starting conversion from PLINK dosage to GEMMA input BIMBAM and anno.")
for i in chrs_to_test:
    anno = open(anno_path + "anno" + str(i) + ".txt", "w")
    BIMBAM = gzip.open(BIMBAM_path + "chr" + str(i) + ".txt.gz", "wb")
    for line in gzip.open(dosage_path + "chr" + str(i) + ".maf0.01.r20.8.dosage.txt.gz", "rb"):
    #line = gzip.open(dosage_path + "chr" + str(i) + ".maf0.01.r20.8.dosage.txt.gz", "rb").readline()    
        arr = line.strip().split()
        (chr, rs, bp, A1, A2, MAF) = arr[0:6]
        #if len(A1) < 2 and len(A2) < 2 and rs in PX_SNPs
        if args.local_anc_SNPs is not None: #abbreviates to only SNPs to be tested around genes
            if len(A1) < 2 and len(A2) < 2 and rs in local_anc_SNPs:
                #force to be in order of local anc samples
                if args.local_anc_samples is not None:
                    dosages = pd.DataFrame(arr[6:])
                    dosages = pd.concat([dosages, dosage_samples], axis = 1)
                    dosages = dosages.set_index(1)
                    dosages = local_anc_samples.join(dosages)
                    dosages = dosages[0].tolist()
                else:
                    dosages = arr[6:]
            
                dosages_str = '\t'.join(dosages)
                BIMBAM_format = (rs + "\t" + A1 + "\t" + A2 + "\t" + dosages_str + "\n")
                BIMBAM.write(BIMBAM_format)
                anno.write(rs + "\t" + bp + "\t" + str(i) + "\n")
        else:
            if len(A1) < 2 and len(A2) < 2:
                #force to be in order of local anc samples
                if args.local_anc_samples is not None:
                    dosages = pd.DataFrame(arr[6:])
                    dosages = pd.concat([dosages, dosage_samples], axis = 1)
                    dosages = dosages.set_index(1)
                    dosages = local_anc_samples.join(dosages)
                    dosages = dosages[0].tolist()
            else:
                dosages = arr[6:]
                
                dosages_str = '\t'.join(dosages)
                BIMBAM_format = (rs + "\t" + A1 + "\t" + A2 + "\t" + dosages_str + "\n")
                BIMBAM.write(BIMBAM_format)
                anno.write(rs + "\t" + bp + "\t" + str(i) + "\n")
    BIMBAM.close()
    anno.close()
    print("Completed with chr " + str(i) + ".")
print("Conversion from PLINK dosage to GEMMA input has been completed. Have a nice day!")
