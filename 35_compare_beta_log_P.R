library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
"%&%" = function(a,b) paste(a,b,sep="")
setwd("/home/angela/px_his_chol/MESA_compare/")

three_pops <- fread("/home/angela/px_his_chol/MESA_compare/three_pops_filtered.csv")
thres <- c(1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)

#beta per pop
for(j in c("CHOL_rank", "HDL_rank", "TRIG_rank", "LDL_rank")){
  three_pops_filtered <- subset(three_pops, pheno == j)
  beta_results <- data.frame(MESA_HIS_beta = numeric(), MESA_CAU_beta = numeric(), thres = numeric())
  beta_df <- three_pops_filtered %>% dplyr::select("HCHS_P", "HCHS_beta", "MESA_HIS_beta", "MESA_CAU_beta")
  for(i in thres){
    thres_df <- subset(beta_df, HCHS_P < i)
    thres_df$HCHS_P <- NULL
    if(nrow(thres_df[complete.cases(thres_df)])){
      cor_df <- as.data.frame(cor(thres_df, use = "complete.obs"))
      cor_df$thres <- i
      cor_df$HCHS_beta <- NULL
      beta_results <- rbind(beta_results, cor_df[1,])
    }
  }
  fwrite(beta_results, j %&% "_cor_beta.csv", row.names = F, quote = F)

  #logP
  log_P_results <- data.frame(MESA_HIS_log_P = numeric(), MESA_CAU_log_P = numeric(), thres = numeric())
  log_P_df <- three_pops_filtered %>% dplyr::select("HCHS_P", "HCHS_log_P", "MESA_HIS_log_P", "MESA_CAU_log_P")
  for(i in thres){
    thres_df <- subset(log_P_df, HCHS_P < i)
    thres_df$HCHS_P <- NULL
    if(nrow(thres_df[complete.cases(thres_df)])){
      cor_df <- as.data.frame(cor(thres_df, use = "complete.obs"))
      cor_df$thres <- i
      cor_df$HCHS_log_P <- NULL
      log_P_results <- rbind(log_P_results, cor_df[1,])
    }
  }
  fwrite(log_P_results, j %&% "_cor_log_P.csv", row.names = F, quote = F)
} 

#beta overall
beta_results <- data.frame(MESA_HIS_beta = numeric(), MESA_CAU_beta = numeric(), thres = numeric())
beta_df <- three_pops %>% dplyr::select("HCHS_P", "HCHS_beta", "MESA_HIS_beta", "MESA_CAU_beta")
for(i in thres){
  thres_df <- subset(beta_df, HCHS_P < i)
  thres_df$HCHS_P <- NULL
  if(nrow(thres_df[complete.cases(thres_df)])){
    cor_df <- as.data.frame(cor(thres_df, use = "complete.obs"))
    cor_df$thres <- i
    cor_df$HCHS_beta <- NULL
    beta_results <- rbind(beta_results, cor_df[1,])
  }
}
fwrite(beta_results, "cor_beta.csv", row.names = F, quote = F)

#logP overall
log_P_results <- data.frame(MESA_HIS_log_P = numeric(), MESA_CAU_log_P = numeric(), thres = numeric())
log_P_df <- three_pops %>% dplyr::select("HCHS_P", "HCHS_log_P", "MESA_HIS_log_P", "MESA_CAU_log_P")
for(i in thres){
  thres_df <- subset(log_P_df, HCHS_P < i)
  thres_df$HCHS_P <- NULL
  if(nrow(thres_df[complete.cases(thres_df)])){
    cor_df <- as.data.frame(cor(thres_df, use = "complete.obs"))
    cor_df$thres <- i
    cor_df$HCHS_log_P <- NULL
    log_P_results <- rbind(log_P_results, cor_df[1,])
  }
}
fwrite(log_P_results, "cor_log_P.csv", row.names = F, quote = F)

#PLOTTING
log_P_results$thres <- as.factor(log_P_results$thres)
melted <- melt(log_P_results, id.vars = "thres")

log_P <- ggplot(data = melted, aes(x = thres, y = value, group = variable)) + 
  geom_line(aes(color = variable)) +
  scale_x_discrete(limits = rev(levels(melted$thres))) + 
  scale_color_brewer(palette = "Set1", labels = c("HIS", "CAU")) + 
  labs(title = "HCHS and MESA -log10(P) values (PX)", color = "Pop.", y = "Correlation of -log10(P)", x = "HCHS P < #") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))
png("log_P.png")
print(log_P)
dev.off()

beta_results$thres <- as.factor(beta_results$thres)
melted <- melt(beta_results, id.vars = "thres")

beta <- ggplot(data = melted, aes(x = thres, y = value, group = variable)) + 
  geom_line(aes(color = variable)) +
  scale_x_discrete(limits = rev(levels(melted$thres))) + 
  scale_color_brewer(palette = "Set1", labels = c("HIS", "CAU")) + 
  labs(title = "HCHS and MESA beta values (PX)", color = "Pop.", y = "Correlation of betas", x = "HCHS P < #") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))
png("beta.png")
print(beta)
dev.off()
