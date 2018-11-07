#backward elimination of significant genes to determine which ones are independent
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
sig_gene_HCHS <- fread("/home/angela/px_his_chol/MESA_compare/GTEx_WB/sig_gene_HCHS.csv")
pheno <- fread('/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt', header = T)
pheno <- pheno %>% select(FID, IID, CHOL_rank, HDL_rank, TRIG_rank, LDL_rank)
phenos <- c("CHOL", "HDL", "TRIG", "LDL")
tissues <- fread("/home/angela/px_yri_chol/PrediXcan/database_tissues.txt", header = F)
tissues <- tissues$V1 #get list of tissues
back_elim <- pheno[,2] #extract only FID and IID

#pheno_name <- "CHOL"
#tiss_name <- "TW_Artery_Coronary_0.5.db"

for(pheno_name in phenos){
  test_pheno <- subset(sig_gene_HCHS, pheno == pheno_name) #subset to just the pheno in question
  sig_tiss_gene <- paste(test_pheno$tissue, test_pheno$gene, sep = "_")
  pheno_name_rank <- paste(pheno_name, "_rank", sep = "")
  back_elim <- pheno %>% select(IID, pheno_name_rank)
  print("Running analyses on " %&% pheno_name_rank %&% ".")
  
  for(tiss_name in tissues){
    print("Running analyses on " %&% tiss_name %&% ".")
    tiss <- fread('/home/angela/px_his_chol/PrediXcan/' %&% tiss_name %&% '_predicted_expression.txt', header = F, nThread = 30) #load in predicted expression file
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
  
  print("Started making models for " %&% pheno_name_rank %&% ".")
  back_elim$IID <- NULL #the IID column has served its purpose
  back_elim <- back_elim[complete.cases(back_elim),]
  predictor_genes <- colnames(back_elim)[2:length(colnames(back_elim))] #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
  fmla <- as.formula(paste(pheno_name_rank,  " ~ ", paste(predictor_genes, collapse= "+")))
  fwrite(back_elim, "_all_tiss_gene_before_back_elim.csv", row.names = F, col.names = T, sep = ",", na = NA, quote = F)
  all_tiss_gene <- lm(fmla, data = back_elim) 
  saveRDS(all_tiss_gene, file = pheno_name %&% "_all_tiss_gene.rds")
  print("Finished making full model for " %&% pheno_name_rank %&% ".")
  back_elim_complete <- step(all_tiss_gene, direction = "backward", trace = FALSE) #perform backward analysis on full model
  saveRDS(back_elim_complete, file = pheno_name %&% "_back_elim.rds")
  print("Finished making backward-eliminated model for " %&% pheno_name_rank %&% ".")
}


