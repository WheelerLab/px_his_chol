#process COLOC output and join to novel significant gene associations
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

sig_gene_HCHS <- fread("/home/angela/px_his_chol/MESA_compare/sig_PX_results_MESA.csv")
sig_gene_HCHS$gene <- gsub("\\..", "", sig_gene_HCHS$gene)
phenos <- c("CHOL", "HDL", "TRIG", "LDL")
pops <- c("AFA", "AFHI", "ALL", "CAU", "HIS")
tiss_abbv <- fread("/home/angela/tiss_abbreviations.txt")
HCHS_COLOC_output <- data.frame(gene_id = character(), p3 = character(), p4 = character(), tissue = character(), pheno = character(), stringsAsFactors = F) 

#MESA pops
for(pheno in phenos){
  for(pop in pops){
    COLOC_output <- fread("/home/angela/px_his_chol/COLOC/results/" %&% pop %&% "_" %&% pheno %&% "_rank.txt.gz", header = T)
    COLOC_output <- COLOC_output %>% dplyr::select(gene_id, p3, p4)
    COLOC_output$tissue <- paste(pop, ".db", sep = "")
    COLOC_output$pheno <- pheno
    COLOC_output$gene_id <- gsub("\\..", "", COLOC_output$gene_id)
    HCHS_COLOC_output <- rbind(HCHS_COLOC_output, COLOC_output)
  }
} 

#GTEx tiss
tissues <- c("TW_Adipose_Subcutaneous_0.5.db", "TW_Adipose_Visceral_Omentum_0.5.db", "TW_Adrenal_Gland_0.5.db", "TW_Artery_Aorta_0.5.db", "TW_Artery_Coronary_0.5.db", "TW_Artery_Tibial_0.5.db", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "TW_Brain_Caudate_basal_ganglia_0.5.db", "TW_Brain_Cerebellar_Hemisphere_0.5.db", "TW_Brain_Cerebellum_0.5.db", "TW_Brain_Cortex_0.5.db", "TW_Brain_Frontal_Cortex_BA9_0.5.db", "TW_Brain_Hippocampus_0.5.db", "TW_Brain_Hypothalamus_0.5.db", "TW_Brain_Nucleus_accumbens_basal_ganglia_0.5.db", "TW_Brain_Putamen_basal_ganglia_0.5.db", "TW_Breast_Mammary_Tissue_0.5.db", "TW_Cells_EBV-transformed_lymphocytes_0.5.db", "TW_Cells_Transformed_fibroblasts_0.5.db", "TW_Colon_Sigmoid_0.5.db", "TW_Colon_Transverse_0.5.db", "TW_Esophagus_Gastroesophageal_Junction_0.5.db", "TW_Esophagus_Mucosa_0.5.db", "TW_Esophagus_Muscularis_0.5.db", "TW_Heart_Atrial_Appendage_0.5.db", "TW_Heart_Left_Ventricle_0.5.db", "TW_Liver_0.5.db", "TW_Lung_0.5.db", "TW_Muscle_Skeletal_0.5.db", "TW_Nerve_Tibial_0.5.db", "TW_Ovary_0.5.db", "TW_Pancreas_0.5.db", "TW_Pituitary_0.5.db", "TW_Prostate_0.5.db", "TW_Skin_Not_Sun_Exposed_Suprapubic_0.5.db", "TW_Skin_Sun_Exposed_Lower_leg_0.5.db", "TW_Small_Intestine_Terminal_Ileum_0.5.db", "TW_Spleen_0.5.db", "TW_Stomach_0.5.db", "TW_Testis_0.5.db", "TW_Thyroid_0.5.db", "TW_Uterus_0.5.db", "TW_Vagina_0.5.db", "TW_Whole_Blood_0.5.db")
for(pheno in phenos){
  for(tiss in tissues){
    COLOC_output <- fread("/home/angela/px_his_chol/COLOC/results/" %&% tiss %&% "_" %&% pheno %&% "_rank.txt.gz")
    COLOC_output <- COLOC_output %>% dplyr::select(gene_id, p3, p4)
    COLOC_output$tissue <- tiss
    COLOC_output$pheno <- pheno
    COLOC_output$gene_id <- gsub("\\..", "", COLOC_output$gene_id)
    HCHS_COLOC_output <- rbind(HCHS_COLOC_output, COLOC_output)
  }
}

HCHS_COLOC_output$gene_id <- as.character(HCHS_COLOC_output$gene_id)
HCHS_COLOC_output$p3 <- as.numeric(HCHS_COLOC_output$p3)
HCHS_COLOC_output$p4 <- as.numeric(HCHS_COLOC_output$p4)
HCHS_COLOC_output$tissue <- as.character(HCHS_COLOC_output$tissue)
HCHS_COLOC_output$pheno <- as.character(HCHS_COLOC_output$pheno) #sometimes R makes them weird types
colnames(HCHS_COLOC_output) <- c("gene", "P3", "P4", "tissue", "pheno")
HCHS_COLOC_output <- left_join(HCHS_COLOC_output, tiss_abbv, by = "tissue")
HCHS_COLOC_output$tissue <- HCHS_COLOC_output$tiss
HCHS_COLOC_output$tiss <- NULL 

sig_gene_HCHS_COLOC <- left_join(sig_gene_HCHS, HCHS_COLOC_output, by = c("gene", "tissue", "pheno"))
fwrite(sig_gene_HCHS_COLOC, "/home/angela/px_his_chol/COLOC/results/sig_gene_HCHS_COLOC_P4_all.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)
sig_gene_HCHS_COLOC <- subset(sig_gene_HCHS_COLOC, P4 >= 0.5)
sig_gene_HCHS_COLOC <- sig_gene_HCHS_COLOC[order(sig_gene_HCHS_COLOC$pheno, sig_gene_HCHS_COLOC$genename),]
fwrite(sig_gene_HCHS_COLOC, "/home/angela/px_his_chol/COLOC/results/sig_gene_HCHS_COLOC_P4_0.5.csv", sep = ",", quote = F, na = "NA", row.names = F, col.names = T)
