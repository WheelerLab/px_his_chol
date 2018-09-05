#Takes information from PrediXcan-style dosages to make SNP annotation and BIMBAM files for GEMMA
import gzip
import numpy

mount = "/home/angela"
dosage_path = mount + "/px_his_chol/Imputation/UMich/UMich_dosages/"
BIMBAM_path = mount + "/px_his_chol/local_anc_GEMMA/BIMBAM/"
anno_path = mount + "/px_his_chol/local_anc_GEMMA/anno/"
PX_SNPs = set(numpy.loadtxt(mount + "/px_his_chol/local_anc_GEMMA/PrediXcan_SNPs.txt", dtype = 'string'))

for i in range(1,23):
    anno = open(anno_path + "anno" + str(i) + ".txt", "w")
    BIMBAM = gzip.open(BIMBAM_path + "chr" + str(i) + ".txt.gz", "wb")
    for line in gzip.open(dosage_path + "chr" + str(i) + ".maf0.01.r20.8.dosage.txt.gz", "rb"):
        arr = line.strip().split()
        (chr, rs, bp, A1, A2, MAF) = arr[0:6]
        if len(A1) < 2 and len(A2) < 2 and rs in PX_SNPs:
            dosages = arr[7:] #skip problematic person #1
            dosages_str = '\t'.join(dosages)
            BIMBAM_format = (rs + "\t" + A1 + "\t" + A2 + "\t" + dosages_str + "\n")
            BIMBAM.write(BIMBAM_format)
            anno.write(rs + "\t" + bp + "\t" + chr + "\n")
    BIMBAM.close()
    anno.close()
    print("Completed with chr " + str(i) + ".")
