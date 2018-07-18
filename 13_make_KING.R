#before running this, make king.kin0 and king.kin using KING
library(data.table)
library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(OmicKriging)
system("/home/angela/px_his_chol/KING/king -b /home/angela/px_his_chol/imputed_PLINK/HCHS_pedigree_chr22.bed --kinship --prefix /home/angela/px_his_chol/imputed_PLINK/relatedness")
snpgdsBED2GDS(bed.fn = "/home/angela/px_his_chol/imputed_PLINK/HCHS_pedigree_chr22.bed", bim.fn = "/home/angela/px_his_chol/imputed_PLINK/HCHS_pedigree_chr22.bim", fam.fn = "/home/angela/px_his_chol/imputed_PLINK/HCHS_pedigree_chr22.fam", out.gdsfn = "/home/angela/px_his_chol/imputed_PLINK/HCHS_pedigree_chr22.gds", family = T)
geno <- GdsGenotypeReader(filename = "/home/angela/px_his_chol/imputed_PLINK/HCHS_pedigree_chr22.gds")
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)
KINGmat <- king2mat(file.kin0 = "/home/angela/px_his_chol/imputed_PLINK/relatedness.kin0", file.kin = "/home/angela/px_his_chol/imputed_PLINK/relatedness.kin", iids = iids)
KINGmat2 <- as.data.frame(KINGmat)
colnames(KINGmat2) <- iids
rownames(KINGmat2) <- iids
fwrite(KINGmat2, "/home/angela/px_his_chol/imputed_PLINK/relatedness_matrix_wIID.txt", sep = "\t", row.names = T, col.names = T, nThread = 40)
fwrite(KINGmat2, "/home/angela/px_his_chol/imputed_PLINK/relatedness_matrix_woIID.txt", sep = "\t", row.names = F, col.names = F, nThread = 40)

#for GCTA b/c they don't like neg. values
KINGmat3 <- KINGmat2
KINGmat3[KINGmat3 < 0] <- 0
fwrite(KINGmat3, "/home/angela/px_his_chol/imputed_PLINK/relatedness_matrix_wIID_noNeg.txt", sep = "\t", row.names = T, col.names = T, nThread = 40)
fwrite(KINGmat3, "/home/angela/px_his_chol/imputed_PLINK/relatedness_matrix_woIID_noNeg.txt", sep = "\t", row.names = F, col.names = F, nThread = 40)

#then GCTA style
KINGmat3 <- as.matrix(KINGmat3)
rownames(KINGmat3) <- iids
colnames(KINGmat3) <- iids
write_GRMBin(KINGmat3, prefix = "/home/angela/px_his_chol/imputed_PLINK/relatedness")
