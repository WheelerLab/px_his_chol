#make admixture mapping plot for gene in input
library(data.table)
library(dplyr)
library(ecoflux)
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
  GWAS <- GWAS[order(GWAS$predictor, decreasing = T),] 
  
  pdf("/home/angela/px_his_chol/Manuscript_figures/" %&% gene %&% "_chr" %&% chr %&% "_" %&% pheno %&% "_" %&% tissue %&% ".pdf", width = 6, height = 4)
  #png("/home/angela/px_his_chol/Manuscript_figures/" %&% gene %&% "_chr" %&% chr %&% "_" %&% pheno %&% "_" %&% tissue %&% ".png")#, width = 6, height = 4)
  admixture_plot <- ggplot() + 
    geom_point(data = GWAS, aes(x = ps, y = -log10(p_wald), fill = predictor), shape = 21, stroke = 0, size = 2, color = "transparent") +
    geom_line(data = NAT, aes(x = ps, y = -log10(p_wald), color = "NAT")) +
    geom_line(data = IBS, aes(x = ps, y = -log10(p_wald), color = "IBS")) +
    geom_line(data = YRI, aes(x = ps, y = -log10(p_wald), color = "YRI")) + 
    coord_cartesian(xlim = c(min(subset(GWAS, predictor == "In predictor")$ps), max(subset(GWAS, predictor == "In predictor")$ps)), ylim = c(0, 10)) + 
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) + #https://stackoverflow.com/questions/42323247/how-to-force-axis-values-to-scientific-notation-in-ggplot?rq=1
    xlab("bp") + 
    ylab("-log10(p)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 15)) + 
    guides(color = guide_legend(title = "Ancestry"), fill = guide_legend(title = "Location")) + 
    scale_color_brewer(type = "div", palette = "Set1", direction = -1) + 
    scale_fill_grey()
  print(admixture_plot)
  dev.off()
}

#HDL-16-COQ9-OVARY
plot_admixture_mapping("16", "HDL_rank", "TW_Ovary_0.5.db", "COQ9")
