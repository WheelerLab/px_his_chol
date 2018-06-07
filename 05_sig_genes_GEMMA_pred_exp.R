#Finding most significant genes from GEMMA predicted expression
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")
pheno_list <- c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")
tiss_list <- c("TW_Adrenal_Gland_0.5.db", "TW_Artery_Aorta_0.5.db", "TW_Artery_Coronary_0.5.db", "TW_Artery_Tibial_0.5.db", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "TW_Brain_Caudate_basal_ganglia_0.5.db", "TW_Brain_Cerebellar_Hemisphere_0.5.db", "TW_Brain_Cerebellum_0.5.db", "TW_Brain_Cortex_0.5.db", "TW_Brain_Frontal_Cortex_BA9_0.5.db", "TW_Brain_Hippocampus_0.5.db", "TW_Brain_Hypothalamus_0.5.db", "TW_Brain_Nucleus_accumbens_basal_ganglia_0.5.db", "TW_Brain_Putamen_basal_ganglia_0.5.db", "TW_Breast_Mammary_Tissue_0.5.db", "TW_Cells_EBV-transformed_lymphocytes_0.5.db", "TW_Cells_Transformed_fibroblasts_0.5.db", "TW_Colon_Sigmoid_0.5.db", "TW_Colon_Transverse_0.5.db", "TW_Esophagus_Gastroesophageal_Junction_0.5.db", "TW_Esophagus_Mucosa_0.5.db", "TW_Esophagus_Muscularis_0.5.db", "TW_Heart_Atrial_Appendage_0.5.db", "TW_Heart_Left_Ventricle_0.5.db", "TW_Liver_0.5.db", "TW_Lung_0.5.db", "TW_Muscle_Skeletal_0.5.db", "TW_Nerve_Tibial_0.5.db", "TW_Ovary_0.5.db", "TW_Pancreas_0.5.db", "TW_Pituitary_0.5.db", "TW_Prostate_0.5.db", "TW_Skin_Not_Sun_Exposed_Suprapubic_0.5.db", "TW_Skin_Sun_Exposed_Lower_leg_0.5.db", "TW_Small_Intestine_Terminal_Ileum_0.5.db", "TW_Spleen_0.5.db", "TW_Stomach_0.5.db", "TW_Testis_0.5.db", "TW_Thyroid_0.5.db", "TW_Uterus_0.5.db", "TW_Vagina_0.5.db", "TW_Whole_Blood_0.5.db", "AFA.db", "CAU.db", "HIS.db", "AFHI.db", "ALL.db")
GEMMA_output <- "/home/angela/px_his_chol/GEMMA/pred_exp/output/"
ChrENGene <- read.table('/home/angela/px_yri_chol/PrediXcan/ChrENGene_forRenaming.txt', header = T)
ChrENGene$gene <- as.character(ChrENGene$gene)
ChrENGene$gene <- gsub("\\..*", "", ChrENGene$gene)

for(pheno in pheno_list){
  first_tiss <- fread(GEMMA_output %&% pheno %&% "_TW_Adipose_Subcutaneous_0.5.db.assoc.txt")
  first_tiss$FDR <- p.adjust(first_tiss$p_wald)
  first_tiss$tissue <- "TW_Adipose_Subcutaneous_0.5.db"
  first_tiss$bon_sig <- 0.05/dim(first_tiss)[1]
  for(tiss in tiss_list){
    next_tiss <- fread(GEMMA_output %&% pheno %&% "_" %&% tiss %&% ".assoc.txt") 
    next_tiss$FDR <- p.adjust(next_tiss$p_wald)
    next_tiss$tissue <- tiss
    next_tiss$bon_sig <- 0.05/dim(next_tiss)[1]
    first_tiss <- rbind(first_tiss, next_tiss)
  }
  first_tiss$rs <- gsub("\\..*", "", first_tiss$rs)
  colnames(first_tiss) <- c("chr", "gene", "ps", "n_miss", "allele1", "allele0", "af", "beta", "se", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score", "FDR", "tissue", "bon_sig")
  genes <- left_join(first_tiss, ChrENGene)
  genes <- genes[c("CHR", "gene", "genename", "tissue", "beta", "se", "l_remle", "l_mle", "p_wald", "FDR", "bon_sig")]
  fwrite(genes, GEMMA_output %&% pheno %&% "_all.csv", sep = ",")
  sig_0.05 <- subset(genes, FDR < 0.05) #genes that have false discovery rate < 0.05
  fwrite(sig_0.05, GEMMA_output %&% pheno %&% "_FDR_0.05.csv", sep = ",")
  sig_bon <- subset(genes, p_wald < bon_sig) #genes that are bonferroni significant
  fwrite(sig_bon, GEMMA_output %&% pheno %&% "_bon_sig.csv", sep = ",")
}


