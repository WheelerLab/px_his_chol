#makes Manhattan plot with ancestry mixture overlays by ancestry
library(data.table)
library(dplyr)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")

plot_admixture_mapping <- function(chr, pheno, tissue, gene){
  GWAS_file <- "/home/angela/px_his_chol/GEMMA/output/" %&% pheno %&% "_sig_snps_p-0.txt"
  NAT_file <- "/home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/output/chr" %&% chr %&% "_" %&% pheno %&% "_NAT.assoc.txt"
  IBS_file <- "/home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/output/chr" %&% chr %&% "_" %&% pheno %&% "_IBS.assoc.txt"
  YRI_file <- "/home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/output/chr" %&% chr %&% "_" %&% pheno %&% "_YRI.assoc.txt"
  SNPs_file <- "/home/angela/px_his_chol/MESA_compare/GTEx_WB/SNPs_n_samples/" %&% tissue %&% "_SNPs.txt"

  GWAS <- fread(GWAS_file, stringsAsFactors = F)
  NAT <- fread(NAT_file, stringsAsFactors = F)
  IBS <- fread(IBS_file, stringsAsFactors = F)
  YRI <- fread(YRI_file, stringsAsFactors = F)
  gene_SNPs <- fread(SNPs_file, stringsAsFactors = F)

  gene_SNPs <- subset(gene_SNPs, V3 == gene)
  SNPs_in_models <- as.data.frame(gene_SNPs$V1)
  colnames(SNPs_in_models) <- "rs"
  SNPs_in_models$predictor <- "In predictor"
  GWAS <- left_join(GWAS, SNPs_in_models)
  GWAS$predictor[is.na(GWAS$predictor)] <- "Not in predictor" #fill not predictor SNPs
  
  pdf("/home/angela/px_his_chol/Manuscript_figures/" %&% gene %&% "_chr" %&% chr %&% "_" %&% pheno %&% "_" %&% tissue %&% ".pdf", width = 6, height = 4)
  #png("/home/angela/px_his_chol/Manuscript_figures/" %&% gene %&% "_chr" %&% chr %&% "_" %&% pheno %&% "_" %&% tissue %&% ".png")#, width = 6, height = 4)
  admixture_plot <- ggplot() + 
    geom_point(data = GWAS, aes(x = ps, y = -log10(p_wald), fill = predictor), shape = 21, stroke = 0, size = 3, color = "transparent") +
    geom_line(data = NAT, aes(x = ps, y = -log10(p_wald), color = "NAT")) +
    geom_line(data = IBS, aes(x = ps, y = -log10(p_wald), color = "IBS")) +
    geom_line(data = YRI, aes(x = ps, y = -log10(p_wald), color = "YRI")) + 
    coord_cartesian(xlim = c(min(subset(GWAS, predictor == "In predictor")$ps), max(subset(GWAS, predictor == "In predictor")$ps)), ylim = c(0, 10)) + 
    xlab("bp") + 
    ylab("-log10(p)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + 
    guides(color = guide_legend(title = "Ancestry"), fill = guide_legend(title = "Location")) + 
    scale_color_brewer(type = "div", palette = "Dark2") + 
    scale_fill_brewer(palette = "Accent", direction = 1)
  print(admixture_plot)
  dev.off()
}

#CHOL-16-CCL22-ESPMSL
plot_admixture_mapping("16", "CHOL_rank", "TW_Esophagus_Muscularis_0.5.db", "CCL22")

#HDL-19-CCDC159-FIBRBLS
plot_admixture_mapping("19", "HDL_rank", "TW_Cells_Transformed_fibroblasts_0.5.db", "CCDC159")

#LDL-16-HP-BRNACC
plot_admixture_mapping("16", "LDL_rank", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "HP")

#HDL-16-ENKD1-BRNACC
plot_admixture_mapping("16", "HDL_rank", "TW_Brain_Anterior_cingulate_cortex_BA24_0.5.db", "ENKD1")

#HDL-16-HSD11B2-BRNACC
plot_admixture_mapping("16", "HDL_rank", "TW_Testis_0.5.db", "HSD11B2")


#TRIG-14-RPGRIP1-LIVER *INSIGNIFICANT IN GLGC
#LDL-11-TMEM258-SKINS
