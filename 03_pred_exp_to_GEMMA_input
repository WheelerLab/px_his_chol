#converting predicted expression to GEMMA input
#this is very unnecessarily clunky and can probably done much faster with half the code
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")
#tiss_list <- c("TW_Adrenal_Gland_0.5.db", "TW_Artery_Aorta_0.5.db", "TW_Artery_Coronary_0.5.db", "TW_Artery_Tibial_0.5.db", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "TW_Brain_Caudate_basal_ganglia_0.5.db", "TW_Brain_Cerebellar_Hemisphere_0.5.db", "TW_Brain_Cerebellum_0.5.db", "TW_Brain_Cortex_0.5.db", "TW_Brain_Frontal_Cortex_BA9_0.5.db", "TW_Brain_Hippocampus_0.5.db", "TW_Brain_Hypothalamus_0.5.db", "TW_Brain_Nucleus_accumbens_basal_ganglia_0.5.db", "TW_Brain_Putamen_basal_ganglia_0.5.db", "TW_Breast_Mammary_Tissue_0.5.db", "TW_Cells_EBV-transformed_lymphocytes_0.5.db", "TW_Cells_Transformed_fibroblasts_0.5.db", "TW_Colon_Sigmoid_0.5.db", "TW_Colon_Transverse_0.5.db", "TW_Esophagus_Gastroesophageal_Junction_0.5.db", "TW_Esophagus_Mucosa_0.5.db", "TW_Esophagus_Muscularis_0.5.db", "TW_Heart_Atrial_Appendage_0.5.db", "TW_Heart_Left_Ventricle_0.5.db", "TW_Liver_0.5.db", "TW_Lung_0.5.db", "TW_Muscle_Skeletal_0.5.db", "TW_Nerve_Tibial_0.5.db", "TW_Ovary_0.5.db", "TW_Pancreas_0.5.db", "TW_Pituitary_0.5.db", "TW_Prostate_0.5.db", "TW_Skin_Not_Sun_Exposed_Suprapubic_0.5.db", "TW_Skin_Sun_Exposed_Lower_leg_0.5.db", "TW_Small_Intestine_Terminal_Ileum_0.5.db", "TW_Spleen_0.5.db", "TW_Stomach_0.5.db", "TW_Testis_0.5.db", "TW_Thyroid_0.5.db", "TW_Uterus_0.5.db", "TW_Vagina_0.5.db", "TW_Whole_Blood_0.5.db", "AFA.db", "CAU.db", "HIS.db", "AFHI.db", "ALL.db")
tiss_list <- c("AFA.db", "CAU.db", "HIS.db", "AFHI.db", "ALL.db")

for(tiss in tiss_list){
  pred_exp <- fread("/home/angela/px_his_chol/MESAreplication/PrediXcan/" %&% tiss %&% "_predicted_expression.txt", header = F)
  pred_exp$V1 <- NULL
  pred_exp$V2 <- NULL
  pred_exp <- transpose(pred_exp)
  gene_list <- as.data.frame(pred_exp$V1)
  colnames(gene_list) <- "V0"
  pred_exp <- cbind(gene_list, pred_exp)
  pred_exp$V0 <- NA
  pred_exp$V1 <- NA
  colnames(gene_list) <- "V"
  pred_exp <- cbind(gene_list, pred_exp)
  fwrite(pred_exp, "/home/angela/px_his_chol/MESAreplication/GEMMA/pred_exp/" %&% tiss, na = "NA", sep = "\t", col.names = F)
}
