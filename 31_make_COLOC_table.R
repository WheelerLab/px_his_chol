#process COLOC output and join to novel significant gene associations
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

sig_gene_HCHS_novel <- fread("/home/angela/px_his_chol/COLOC/sig_gene_HCHS_novel.csv")
phenos <- c("CHOL", "HDL", "TRIG", "LDL")
pops <- c("AFA.db", "AFHI.db", "ALL.db", "CAU.db", "HIS.db")
HCHS_novel <- data.frame(gene_id = character(), p4 = character(), tissue = character(), pheno = character()) 

for(pheno in phenos){
  for(pop in pops){
    COLOC_output <- fread("/home/angela/px_his_chol/COLOC/results/" %&% pop %&% "_" %&% pheno %&% "_rank.txt.gz", header = T)
    COLOC_output <- COLOC_output %>% select(gene_id, p4)
    COLOC_output$tissue <- pop
    COLOC_output$pheno <- pheno
    COLOC_output$gene_id <- gsub("\\..", "", COLOC_output$gene_id)
    HCHS_novel <- rbind(HCHS_novel, COLOC_output)
  }
}

tissues <- c("TW_Adipose_Subcutaneous_0.5.db", "TW_Adipose_Visceral_Omentum_0.5.db", "TW_Adrenal_Gland_0.5.db", "TW_Artery_Aorta_0.5.db", "TW_Artery_Coronary_0.5.db", "TW_Artery_Tibial_0.5.db", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "TW_Brain_Caudate_basal_ganglia_0.5.db", "TW_Brain_Cerebellar_Hemisphere_0.5.db", "TW_Brain_Cerebellum_0.5.db", "TW_Brain_Cortex_0.5.db", "TW_Brain_Frontal_Cortex_BA9_0.5.db", "TW_Brain_Hippocampus_0.5.db", "TW_Brain_Hypothalamus_0.5.db", "TW_Brain_Nucleus_accumbens_basal_ganglia_0.5.db", "TW_Brain_Putamen_basal_ganglia_0.5.db", "TW_Breast_Mammary_tissuesue_0.5.db", "TW_Cells_EBV-transformed_lymphocytes_0.5.db", "TW_Cells_Transformed_fibroblasts_0.5.db", "TW_Colon_Sigmoid_0.5.db", "TW_Colon_Transverse_0.5.db", "TW_Esophagus_Gastroesophageal_Junction_0.5.db", "TW_Esophagus_Mucosa_0.5.db", "TW_Esophagus_Muscularis_0.5.db", "TW_Heart_Atrial_Appendage_0.5.db", "TW_Heart_Left_Ventricle_0.5.db", "TW_Liver_0.5.db", "TW_Lung_0.5.db", "TW_Muscle_Skeletal_0.5.db", "TW_Nerve_Tibial_0.5.db", "TW_Ovary_0.5.db", "TW_Pancreas_0.5.db", "TW_Pituitary_0.5.db", "TW_Prostate_0.5.db", "TW_Skin_Not_Sun_Exposed_Suprapubic_0.5.db", "TW_Skin_Sun_Exposed_Lower_leg_0.5.db", "TW_Small_Intestine_Terminal_Ileum_0.5.db", "TW_Spleen_0.5.db", "TW_Stomach_0.5.db", "TW_Testis_0.5.db", "TW_Thyroid_0.5.db", "TW_Uterus_0.5.db", "TW_Vagina_0.5.db", "TW_Whole_Blood_0.5.db")
for(pheno in phenos){
  for(tiss in tissues){
    COLOC_output <- fread("/home/angela/px_his_chol/COLOC/results/" %&% tiss %&% "_" %&% pheno %&% "_rank.txt.gz")
    COLOC_output <- COLOC_output %>% select(gene_id, p4)
    COLOC_output$tissue <- tiss
    COLOC_output$pheno <- pheno
    COLOC_output$gene_id <- gsub("\\..", "", COLOC_output$gene_id)
    HCHS_novel <- rbind(HCHS_novel, COLOC_output)
  }
}

HCHS_novel$gene_id <- as.character(HCHS_novel$gene_id)
HCHS_novel$p4 <- as.numeric(HCHS_novel$p4)
HCHS_novel$tissue <- as.character(HCHS_novel$tissue)
HCHS_novel$pheno <- as.character(HCHS_novel$pheno)
HCHS_novel[HCHS_novel == "p4"] <- NA
HCHS_novel <- HCHS_novel[complete.cases(HCHS_novel),]
HCHS_novel_max <- HCHS_novel[order(-p4),] 
HCHS_novel_max$first_occur <- paste(HCHS_novel_max$gene_id, HCHS_novel_max$tissue, HCHS_novel_max$pheno, sep = "_") #instead of checking across three columns for unique, just check one
HCHS_novel_max_first_occur <- HCHS_novel_max[match(unique(HCHS_novel_max$first_occur), HCHS_novel_max$first_occur),]
HCHS_novel_max_first_occur$first_occur <- NULL
colnames(HCHS_novel_max_first_occur) <- c("gene", "P4", "tissue", "pheno")

sig_gene_HCHS_novel_COLOC <- left_join(sig_gene_HCHS_novel, HCHS_novel_max_first_occur, by = c("gene", "tissue", "pheno"))
