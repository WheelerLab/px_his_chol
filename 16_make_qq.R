#makes qq and aggregate gene lists from GEMMA/PrediXcan output
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
pheno_list <- c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")
tiss_list <- c("TW_Adipose_Subcutaneous_0.5.db", "TW_Adipose_Visceral_Omentum_0.5.db", "TW_Adrenal_Gland_0.5.db", "TW_Artery_Aorta_0.5.db", "TW_Artery_Coronary_0.5.db", "TW_Artery_Tibial_0.5.db", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "TW_Brain_Caudate_basal_ganglia_0.5.db", "TW_Brain_Cerebellar_Hemisphere_0.5.db", "TW_Brain_Cerebellum_0.5.db", "TW_Brain_Cortex_0.5.db", "TW_Brain_Frontal_Cortex_BA9_0.5.db", "TW_Brain_Hippocampus_0.5.db", "TW_Brain_Hypothalamus_0.5.db", "TW_Brain_Nucleus_accumbens_basal_ganglia_0.5.db", "TW_Brain_Putamen_basal_ganglia_0.5.db", "TW_Breast_Mammary_Tissue_0.5.db", "TW_Cells_EBV-transformed_lymphocytes_0.5.db", "TW_Cells_Transformed_fibroblasts_0.5.db", "TW_Colon_Sigmoid_0.5.db", "TW_Colon_Transverse_0.5.db", "TW_Esophagus_Gastroesophageal_Junction_0.5.db", "TW_Esophagus_Mucosa_0.5.db", "TW_Esophagus_Muscularis_0.5.db", "TW_Heart_Atrial_Appendage_0.5.db", "TW_Heart_Left_Ventricle_0.5.db", "TW_Liver_0.5.db", "TW_Lung_0.5.db", "TW_Muscle_Skeletal_0.5.db", "TW_Nerve_Tibial_0.5.db", "TW_Ovary_0.5.db", "TW_Pancreas_0.5.db", "TW_Pituitary_0.5.db", "TW_Prostate_0.5.db", "TW_Skin_Not_Sun_Exposed_Suprapubic_0.5.db", "TW_Skin_Sun_Exposed_Lower_leg_0.5.db", "TW_Small_Intestine_Terminal_Ileum_0.5.db", "TW_Spleen_0.5.db", "TW_Stomach_0.5.db", "TW_Testis_0.5.db", "TW_Thyroid_0.5.db", "TW_Uterus_0.5.db", "TW_Vagina_0.5.db", "TW_Whole_Blood_0.5.db", "AFA.db", "CAU.db", "HIS.db", "AFHI.db", "ALL.db")
ChrENGene <- read.table('/home/angela/px_yri_chol/PrediXcan/ChrENGene_forRenaming.txt', header = T)
ChrENGene$gene <- as.character(ChrENGene$gene)
ChrENGene$gene <- gsub("\\..*", "", ChrENGene$gene)
source("/home/angela/px_yri_chol/GWAS/qqman.r")

for(pheno in pheno_list){
  genes <- t(as.data.frame(c("rs", "n_miss", "af", "beta", "se", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score", "tissue")))
  colnames(genes) <- genes[1,]
  genes <- genes[-1,]
  for(tiss in tiss_list){
    GEMMA_output <- fread("/home/angela/px_his_chol/MESAreplication/new_covars/output/" %&% pheno %&% "_" %&% tiss %&% ".assoc.txt")
    GEMMA_output$tissue <- tiss
    GEMMA_output$chr <- NULL
    GEMMA_output$ps <- NULL
    GEMMA_output$allele1 <- NULL
    GEMMA_output$allele0 <- NULL
    GEMMA_output$p_wald <- as.numeric(GEMMA_output$p_wald)
    GEMMA_output <- GEMMA_output[complete.cases(GEMMA_output),]
    png("/home/angela/px_his_chol/MESAreplication/new_covars/output/" %&% pheno %&% "_" %&% tiss %&% ".png")
    qq(GEMMA_output$p_wald, main = "MESA, " %&% tiss %&% "\n" %&% pheno)
    dev.off()
    genes <- rbind(genes, GEMMA_output)
  }
  genes$rs <- gsub("\\..*", "", genes$rs)
  colnames(genes) <- c("gene", "n_miss", "af", "beta", "se", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score", "tissue")
  genes <- left_join(genes, ChrENGene)
  fwrite(genes, "/home/angela/px_his_chol/MESAreplication/new_covars/output/" %&% pheno %&% "_genes.csv", col.names = T, row.names = F, na = "NA", quote = F)
  png("/home/angela/px_his_chol/MESAreplication/new_covars/output/" %&% pheno %&% "_genes.png")
  genes$p_wald <- as.numeric(genes$p_wald)
  qq(genes$p_wald, main = "MESA, all tissues (GEMMA), " %&% pheno)
  dev.off()
}


