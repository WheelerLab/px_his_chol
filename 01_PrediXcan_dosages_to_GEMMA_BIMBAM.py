#defunct, see 02_PrediXcan_dosages_to_GEMMA.py
import gzip
dosage_path = "/home/angela/px_his_chol/Imputation/UMich/UMich_dosages/"
BIMBAM_path = "/home/angela/px_his_chol/GEMMA/BIMBAM/"

#have list of rsids
with open("/home/angela/px_his_chol/SNPs_in_PrediXcan_models.txt") as file:
    rsids = set() #use set b/c order isn't important
    for line in file:
        line = line.strip()
        rsids.add(line)

for i in range(1,23):
    BIMBAM = gzip.open(BIMBAM_path + "chr" + str(i) + ".txt.gz", "wb")
    for line in gzip.open(dosage_path + "chr" + str(i) + ".maf0.01.r20.8.dosage.txt.gz", "rb"):
        arr = line.strip().split()
        (chr, rs, bp, A1, A2, MAF) = arr[0:6]
        if len(A1) < 2 and len(A2) < 2 and rs in rsids: #exclude multi-allelic and SNPs that don't have proper dosages in all and only SNPs w/in PrediXcan models
            dosages = arr[6:]
            dosages_str = '\t'.join(dosages)
            BIMBAM_format = (rs + "\t" + A1 + "\t" + A2 + "\t" + dosages_str + "\n")
            BIMBAM.write(BIMBAM_format)
    BIMBAM.close()
