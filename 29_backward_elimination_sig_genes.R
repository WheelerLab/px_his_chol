#backward elimination of significant genes to determine which ones are independent
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
sig_gene_HCHS <- fread("/home/angela/px_his_chol/MESA_compare/GTEx_WB/sig_gene_HCHS.csv")
sig_tiss_gene <- paste(sig_gene_HCHS$tissue, sig_gene_HCHS$gene, sep = "_")
pheno <- fread('/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt', header = T)
pheno <- pheno %>% select(FID, IID, CHOL_rank, HDL_rank, TRIG_rank, LDL_rank)
phenos <- c("CHOL", "HDL", "TRIG", "LDL")
tissues <- fread("/home/angela/px_yri_chol/PrediXcan/database_tissues.txt", header = F)
tissues <- tissues$V1

pheno_name <- "CHOL"
tiss_name <- "TW_Artery_Coronary_0.5.db"

back_elim <- pheno[,2] #extract only FID and IID

test_pheno <- subset(sig_gene_HCHS, pheno == pheno_name) #subset to just the pheno in question
for(tiss_name in tissues){
  tiss <- fread('/home/angela/px_his_chol/PrediXcan/' %&% tiss_name %&% '_predicted_expression.txt', header = F) #load in predicted expression file
  tiss$V1 <- NULL #remove FID
  tiss_header <- paste(tiss_name, "_", gsub("\\..*", "", c(tiss[1,])), sep = "") #remove all after . and add tiss name
  tiss <- tiss[2:nrow(tiss),] #remove header line
  colnames(tiss) <- tiss_header #assign header line
  colnames(tiss)[1] <- "IID" #rename beginning
  tiss_gene_in_tiss <- intersect(colnames(tiss), sig_tiss_gene) #intersection of tiss-gene combos
  tiss_gene_in_tiss <- c("IID", tiss_gene_in_tiss) #add IID to list
  for_back_elim <- tiss %>% select(tiss_gene_in_tiss) #select IIDs and intersection of tiss-gene combos
  back_elim <- left_join(back_elim, for_back_elim, by = "IID") #add to list to back-elim
}
