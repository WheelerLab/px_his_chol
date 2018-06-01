#makes GEMMA-style GRM using KING output
#before running this, make king.kin0 and king.kn using KING
library(GENESIS)
library(GWASTools)
library(SNPRelate)

snpgdsBED2GDS(bed.fn = "/home/angela/px_his_chol/MESAreplication/HIS_pruned_for_PC.bed", bim.fn = "/home/angela/px_his_chol/MESAreplication/HIS_pruned_for_PC.bim", fam.fn = "/home/angela/px_his_chol/MESAreplication/HIS_pruned_for_PC.fam", out.gdsfn = "/home/angela/px_his_chol/KING/MESAreplication/HIS_2.gds", family = T)
geno <- GdsGenotypeReader(filename = "/home/angela/px_his_chol/KING/MESAreplication/HIS_2.gds")
genoData <- GenotypeData(geno)
iids <- getScanID(genoData)
KINGmat <- king2mat(file.kin0 = "/home/angela/px_his_chol/KING/MESAreplication/king.kin0", file.kin = "/home/angela/px_his_chol/KING/MESAreplication/king.kin", iids = iids)
KINGmat2 <- as.data.frame(KINGmat)
colnames(KINGmat2) <- iids
rownames(KINGmat2) <- iids
fwrite(KINGmat2, "/home/angela/px_his_chol/KING/MESAreplication/relatedness_matrix_wIID.txt", sep = "\t", row.names = T, col.names = T, nThread = 40)
fwrite(KINGmat2, "/home/angela/px_his_chol/KING/MESAreplication/relatedness_matrix_woIID.txt", sep = "\t", row.names = F, col.names = F, nThread = 40)
