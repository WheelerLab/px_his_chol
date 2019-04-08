#process COLOC output and join to novel significant gene associations
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

phenos <- c("CHOL", "HDL", "TRIG", "LDL")
pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")
tiss_abbv <- fread("/home/angela/tiss_abbreviations_MESA_no_db.txt")
colnames(tiss_abbv) <- c("model", "model_2")
tiss_abbv_noTW <- fread("/home/angela/tiss_abbreviations_no_TW.txt")
colnames(tiss_abbv_noTW) <- c("model_3", "model_2")
tiss_abbv <- left_join(tiss_abbv, tiss_abbv_noTW)
HCHS_COLOC_output <- data.frame(gene_id = character(), p3 = character(), p4 = character(), model = character(), pheno = character(), stringsAsFactors = F) 

#MESA pops
for(pheno in phenos){
  for(pop in pops){
    COLOC_output <- fread("/home/angela/px_his_chol/COLOC/results/" %&% pop %&% "_" %&% pheno %&% "_rank.txt.gz", header = T)
    COLOC_output <- COLOC_output %>% dplyr::select(gene_id, p3, p4)
    COLOC_output$model <- pop
    COLOC_output$pheno <- pheno
    HCHS_COLOC_output <- rbind(HCHS_COLOC_output, COLOC_output)
  }
} 

#GTEx tiss
tissues <- c("TW_Adipose_Subcutaneous_0.5.db", "TW_Adipose_Visceral_Omentum_0.5.db", "TW_Adrenal_Gland_0.5.db", "TW_Artery_Aorta_0.5.db", "TW_Artery_Coronary_0.5.db", "TW_Artery_Tibial_0.5.db", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "TW_Brain_Caudate_basal_ganglia_0.5.db", "TW_Brain_Cerebellar_Hemisphere_0.5.db", "TW_Brain_Cerebellum_0.5.db", "TW_Brain_Cortex_0.5.db", "TW_Brain_Frontal_Cortex_BA9_0.5.db", "TW_Brain_Hippocampus_0.5.db", "TW_Brain_Hypothalamus_0.5.db", "TW_Brain_Nucleus_accumbens_basal_ganglia_0.5.db", "TW_Brain_Putamen_basal_ganglia_0.5.db", "TW_Breast_Mammary_Tissue_0.5.db", "TW_Cells_EBV-transformed_lymphocytes_0.5.db", "TW_Cells_Transformed_fibroblasts_0.5.db", "TW_Colon_Sigmoid_0.5.db", "TW_Colon_Transverse_0.5.db", "TW_Esophagus_Gastroesophageal_Junction_0.5.db", "TW_Esophagus_Mucosa_0.5.db", "TW_Esophagus_Muscularis_0.5.db", "TW_Heart_Atrial_Appendage_0.5.db", "TW_Heart_Left_Ventricle_0.5.db", "TW_Liver_0.5.db", "TW_Lung_0.5.db", "TW_Muscle_Skeletal_0.5.db", "TW_Nerve_Tibial_0.5.db", "TW_Ovary_0.5.db", "TW_Pancreas_0.5.db", "TW_Pituitary_0.5.db", "TW_Prostate_0.5.db", "TW_Skin_Not_Sun_Exposed_Suprapubic_0.5.db", "TW_Skin_Sun_Exposed_Lower_leg_0.5.db", "TW_Small_Intestine_Terminal_Ileum_0.5.db", "TW_Spleen_0.5.db", "TW_Stomach_0.5.db", "TW_Testis_0.5.db", "TW_Thyroid_0.5.db", "TW_Uterus_0.5.db", "TW_Vagina_0.5.db", "TW_Whole_Blood_0.5.db")
for(pheno in phenos){
  for(tiss in tissues){
    COLOC_output <- fread("/home/angela/px_his_chol/COLOC/results/" %&% tiss %&% "_" %&% pheno %&% "_rank.txt.gz")
    COLOC_output <- COLOC_output %>% dplyr::select(gene_id, p3, p4)
    COLOC_output$model <- tiss
    COLOC_output$pheno <- pheno
    HCHS_COLOC_output <- rbind(HCHS_COLOC_output, COLOC_output)
  }
}

HCHS_COLOC_output$gene_id <- as.character(HCHS_COLOC_output$gene_id)
HCHS_COLOC_output$p3 <- as.numeric(HCHS_COLOC_output$p3)
HCHS_COLOC_output$p4 <- as.numeric(HCHS_COLOC_output$p4)
HCHS_COLOC_output$model <- as.character(HCHS_COLOC_output$model)
HCHS_COLOC_output$pheno <- as.character(HCHS_COLOC_output$pheno) #sometimes R makes them weird types
colnames(HCHS_COLOC_output) <- c("gene", "P3", "P4", "model", "pheno")
HCHS_COLOC_output <- left_join(HCHS_COLOC_output, tiss_abbv, by = "model")
HCHS_COLOC_output$model <- HCHS_COLOC_output$model_3
HCHS_COLOC_output$model_2 <- NULL 
HCHS_COLOC_output$model_3 <- NULL
fwrite(HCHS_COLOC_output, "/home/angela/px_his_chol/COLOC/results/HCHS_COLOC_all.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)

#restrict to just sig
#HCHS_COLOC_output <- fread("/home/angela/px_his_chol/COLOC/results/HCHS_COLOC_all.csv")
primary_sig_HCHS <- fread("/home/angela/px_his_chol/MESA_compare/primary_sig.csv")
primary_sig_HCHS$pheno <- gsub("_rank", "", primary_sig_HCHS$pheno)
primary_sig_HCHS_COLOC <- left_join(primary_sig_HCHS, HCHS_COLOC_output, by = c("gene", "model", "pheno"))

start_bp <- fread("/home/angela/px_yri_chol/PrediXcan/BP_Chrome.txt")
primary_sig_HCHS_COLOC <- left_join(primary_sig_HCHS_COLOC, start_bp, by = "gene")
primary_sig_HCHS_COLOC <- primary_sig_HCHS_COLOC %>% dplyr::select(CHR, BP, genename, model, pheno, HCHS_beta, HCHS_P, P3, P4)
primary_sig_HCHS_COLOC <- primary_sig_HCHS_COLOC[order(as.numeric(primary_sig_HCHS_COLOC$CHR), as.numeric(primary_sig_HCHS_COLOC$BP)),]
fwrite(primary_sig_HCHS_COLOC, "/home/angela/px_his_chol/COLOC/results/primary_sig_HCHS_COLOC_P4_all.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)
nrow(subset(primary_sig_HCHS_COLOC, P3 <= 0.5))
nrow(subset(primary_sig_HCHS_COLOC, P4 >= 0.5))
primary_sig_HCHS_COLOC_P4 <- subset(primary_sig_HCHS_COLOC, P4 >= 0.5)
primary_sig_HCHS_COLOC_P4 <- primary_sig_HCHS_COLOC_P4[order(primary_sig_HCHS_COLOC_P4$CHR, primary_sig_HCHS_COLOC_P4$BP),]
fwrite(primary_sig_HCHS_COLOC_P4, "/home/angela/px_his_chol/COLOC/results/primary_sig_HCHS_COLOC_P4_0.5.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)
primary_sig_HCHS_COLOC_P3 <- subset(primary_sig_HCHS_COLOC, P3 <= 0.5)
primary_sig_HCHS_COLOC_P3 <- primary_sig_HCHS_COLOC_P3[order(primary_sig_HCHS_COLOC_P3$CHR, primary_sig_HCHS_COLOC_P3$BP),]
fwrite(primary_sig_HCHS_COLOC_P3, "/home/angela/px_his_chol/COLOC/results/primary_sig_HCHS_COLOC_P3_0.5.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)

#sec
secondary_sig_HCHS <- fread("/home/angela/px_his_chol/MESA_compare/secondary_sig.csv")
secondary_sig_HCHS$pheno <- gsub("_rank", "", secondary_sig_HCHS$pheno)
secondary_sig_HCHS_COLOC <- left_join(secondary_sig_HCHS, HCHS_COLOC_output, by = c("gene", "model", "pheno"))

start_bp <- fread("/home/angela/px_yri_chol/PrediXcan/BP_Chrome.txt")
secondary_sig_HCHS_COLOC <- left_join(secondary_sig_HCHS_COLOC, start_bp, by = "gene")
secondary_sig_HCHS_COLOC <- secondary_sig_HCHS_COLOC %>% dplyr::select(CHR, BP, genename, model, pheno, HCHS_beta, HCHS_P, P3, P4)
secondary_sig_HCHS_COLOC <- secondary_sig_HCHS_COLOC[order(as.numeric(secondary_sig_HCHS_COLOC$CHR), as.numeric(secondary_sig_HCHS_COLOC$BP)),]
fwrite(secondary_sig_HCHS_COLOC, "/home/angela/px_his_chol/COLOC/results/secondary_sig_HCHS_COLOC_P4_all.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)
nrow(subset(secondary_sig_HCHS_COLOC, P3 <= 0.5))
nrow(subset(secondary_sig_HCHS_COLOC, P4 >= 0.5))
secondary_sig_HCHS_COLOC_P4 <- subset(secondary_sig_HCHS_COLOC, P4 >= 0.5)
secondary_sig_HCHS_COLOC_P4 <- secondary_sig_HCHS_COLOC_P4[order(secondary_sig_HCHS_COLOC_P4$CHR, secondary_sig_HCHS_COLOC_P4$BP),]
fwrite(secondary_sig_HCHS_COLOC_P4, "/home/angela/px_his_chol/COLOC/results/secondary_sig_HCHS_COLOC_P4_0.5.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)
secondary_sig_HCHS_COLOC_P3 <- subset(secondary_sig_HCHS_COLOC, P3 <= 0.5)
secondary_sig_HCHS_COLOC_P3 <- secondary_sig_HCHS_COLOC_P3[order(secondary_sig_HCHS_COLOC_P3$CHR, secondary_sig_HCHS_COLOC_P3$BP),]
fwrite(secondary_sig_HCHS_COLOC_P3, "/home/angela/px_his_chol/COLOC/results/secondary_sig_HCHS_COLOC_P3_0.5.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)

