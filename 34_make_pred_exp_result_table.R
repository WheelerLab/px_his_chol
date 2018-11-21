#extracts both all and significant (P < 9.654e-6) PrediXcan results; adds gene name and starting bp
library(data.table)
library(dplyr)
"%&%" = function(a,b) paste(a,b,sep="")

phenos <- c("CHOL", "HDL", "TRIG", "LDL")
tissues <- c("AFA.db", "AFHI.db", "ALL.db", "CAU.db", "HIS.db", "TW_Adipose_Subcutaneous_0.5.db", "TW_Adipose_Visceral_Omentum_0.5.db", "TW_Adrenal_Gland_0.5.db", "TW_Artery_Aorta_0.5.db", "TW_Artery_Coronary_0.5.db", "TW_Artery_Tibial_0.5.db", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "TW_Brain_Caudate_basal_ganglia_0.5.db", "TW_Brain_Cerebellar_Hemisphere_0.5.db", "TW_Brain_Cerebellum_0.5.db", "TW_Brain_Cortex_0.5.db", "TW_Brain_Frontal_Cortex_BA9_0.5.db", "TW_Brain_Hippocampus_0.5.db", "TW_Brain_Hypothalamus_0.5.db", "TW_Brain_Nucleus_accumbens_basal_ganglia_0.5.db", "TW_Brain_Putamen_basal_ganglia_0.5.db", "TW_Breast_Mammary_Tissue_0.5.db", "TW_Cells_EBV-transformed_lymphocytes_0.5.db", "TW_Cells_Transformed_fibroblasts_0.5.db", "TW_Colon_Sigmoid_0.5.db", "TW_Colon_Transverse_0.5.db", "TW_Esophagus_Gastroesophageal_Junction_0.5.db", "TW_Esophagus_Mucosa_0.5.db", "TW_Esophagus_Muscularis_0.5.db", "TW_Heart_Atrial_Appendage_0.5.db", "TW_Heart_Left_Ventricle_0.5.db", "TW_Liver_0.5.db", "TW_Lung_0.5.db", "TW_Muscle_Skeletal_0.5.db", "TW_Nerve_Tibial_0.5.db", "TW_Ovary_0.5.db", "TW_Pancreas_0.5.db", "TW_Pituitary_0.5.db", "TW_Prostate_0.5.db", "TW_Skin_Not_Sun_Exposed_Suprapubic_0.5.db", "TW_Skin_Sun_Exposed_Lower_leg_0.5.db", "TW_Small_Intestine_Terminal_Ileum_0.5.db", "TW_Spleen_0.5.db", "TW_Stomach_0.5.db", "TW_Testis_0.5.db", "TW_Thyroid_0.5.db", "TW_Uterus_0.5.db", "TW_Vagina_0.5.db", "TW_Whole_Blood_0.5.db")

#for gene names, start positions, and tissue abbreviations
ChrENGene <- read.table('/home/angela/px_yri_chol/PrediXcan/ChrENGene_forRenaming.txt', header = T) #for adding gene name
ChrENGene$gene <- as.character(ChrENGene$gene)
ChrENGene$gene <- gsub("\\..*", "", ChrENGene$gene) #remove decimal and after
BP_Chrome <- read.table("/home/angela/px_yri_chol/PrediXcan/BP_Chrome.txt", header = T, stringsAsFactors = F) #for adding starting bp
BP_Chrome <- BP_Chrome %>% dplyr::select(gene_name, BP)
colnames(BP_Chrome) <- c("genename", "start_bp")
BP_Chrome$genename <- as.character(BP_Chrome$genename)
tiss_abb <- fread("/home/angela/tiss_abbreviations.txt") $for adding tissue abbreviations
PX_results <- data.frame(chr = numeric(), gene = character(), genename = character(), start_bp = numeric(), pheno = character(), tissue = character(), p_wald = numeric(), beta = numeric(), se = numeric(), stringsAsFactors = F)  #df for collecting all the results

for(pheno in phenos){
  for(tissue in tissues){
    pred_exp_output <- fread("/home/angela/px_his_chol/GEMMA/pred_exp/output/" %&% pheno %&% "_rank_" %&% tissue %&% ".assoc.txt") 
    pred_exp_output <- pred_exp_output %>% dplyr::select(rs, beta, se, p_wald)
    pred_exp_output$dec <- pred_exp_output$rs
    pred_exp_output$rs <- gsub("\\..*", "", pred_exp_output$rs)
    pred_exp_output <- left_join(pred_exp_output, ChrENGene, by = c("rs" = "gene")) #add gene name
    pred_exp_output$rs <- NULL
    pred_exp_output$pheno <- pheno
    pred_exp_output$tissue <- tissue
    pred_exp_output <- left_join(pred_exp_output, tiss_abb, by = "tissue") #add tissue abbreviations
    pred_exp_output$genename <- as.character(pred_exp_output$genename)
    pred_exp_output <- left_join(pred_exp_output, BP_Chrome, by = "genename") #add starting base pair
    pred_exp_output <- pred_exp_output %>% dplyr::select(CHR, start_bp, dec, genename, pheno, tiss, p_wald, beta, se) #extract necessary data
    colnames(pred_exp_output)[1:6] <- c("chr", "start_bp", "gene", "genename","pheno", "tissue")
    PX_results <- rbind(PX_results, pred_exp_output)
  }
}

#write various iterations of signficance
fwrite(PX_results, "/home/angela/px_his_chol/GEMMA/pred_exp/output/all_PX_results.csv", row.names = F, col.names = T, sep = ",", quote = F, na = "NA")
sig_PX_results <- subset(PX_results, p_wald <= 0.05)
fwrite(sig_PX_results, "/home/angela/px_his_chol/GEMMA/pred_exp/output/0.05_PX_results.csv", row.names = F, col.names = T, sep = ",", quote = F, na = "NA")
sig_PX_results <- subset(PX_results, p_wald <= 9.654e-6)
sig_PX_results <- sig_PX_results[order(sig_PX_results$chr, sig_PX_results$start_bp, sig_PX_results$genename),]
fwrite(sig_PX_results, "/home/angela/px_his_chol/GEMMA/pred_exp/output/sig_PX_results.csv", row.names = F, col.names = T, sep = ",", quote = F, na = "NA")
