library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
"%&%" = function(a,b) paste(a,b,sep="")
setwd("/home/angela/px_his_chol/MESA_compare/")

#4/12 - Separate by model instead of by pheno

three_pops <- fread("full_join_MESA_HIS_CAU.csv")
MESA_models <- c("AFA", "CAU", "HIS", "AFHI", "ALL")
three_pops <- subset(three_pops, model %in% MESA_models)
thres <- c(1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
phenos <- c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")

#beta per pheno
beta_results_overall <- data.frame(MESA_HIS_beta = numeric(), MESA_CAU_beta = numeric(), thres = numeric(), pheno = character(), model = character())
for(k in MESA_models){
  #for(j in 1:length(phenos)){
  #j <- 1
  #three_pops_filtered <- subset(three_pops, pheno == phenos[j] & model == k)
  three_pops_filtered <- subset(three_pops, model == k)
  beta_results <- data.frame(MESA_HIS_beta = numeric(), MESA_CAU_beta = numeric(), thres = numeric())
  beta_df <- three_pops_filtered %>% dplyr::select("HCHS_P", "HCHS_beta", "MESA_HIS_beta", "MESA_CAU_beta")
  for(i in thres){
    thres_df <- subset(beta_df, HCHS_P < i)
    thres_df$HCHS_P <- NULL
    if(nrow(thres_df[complete.cases(thres_df)]) > 1){
      cor_df <- as.data.frame(cor(thres_df, use = "complete.obs", method = "spearman"))
      cor_df$thres <- i
      cor_df$HCHS_beta <- NULL
      beta_results <- rbind(beta_results, cor_df[1,])
    }
  }
  #fwrite(beta_results, phenos[j] %&% "_cor_beta.csv", row.names = F, quote = F)
  melted <- melt(beta_results, id.vars = "thres")
  #melted$pheno <- gsub("_rank", "", phenos[j])
  melted$model <- k
  beta_results_overall <- rbind(beta_results_overall, melted)
  #} 
  
  
}
beta_results_overall$thres <- as.factor(beta_results_overall$thres)
#beta_results_overall$pheno <- as.factor(beta_results_overall$pheno)
#beta_results_overall$pheno <- factor(beta_results_overall$pheno , levels = c("CHOL", "HDL", "TRIG", "LDL"))
beta_results_overall$model <- factor(beta_results_overall$model , levels = c("AFA", "CAU", "HIS", "AFHI", "ALL"))
beta_results_overall$variable <- as.character(beta_results_overall$variable)

beta <- ggplot(data = beta_results_overall, aes(x = thres, y = value, group = variable)) + 
  geom_line(aes(color = variable)) +
  scale_x_discrete(limits = rev(levels(beta_results_overall$thres))) + 
  #scale_color_brewer(palette = "Set1", labels = c("CAU", "HIS")) + 
  scale_color_viridis(labels = c("CAU", "HIS"), discrete = T, option = "D", end = 0.8) + 
  facet_wrap(~ model, nrow = 1) + 
  labs(color = "Pop.", y = "Correlation of betas", x = "HCHS P < #") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1))

pdf("/home/angela/px_his_chol/Manuscript_figures/Fig4.pdf", width = 8, height = 3)
print(beta)
dev.off()

tiff("/home/angela/px_his_chol/Manuscript_figures/Fig4.tiff", width = 15.24, height = 7.62, units = 'cm', res = 300, compression = 'lzw')
print(beta)
dev.off()
