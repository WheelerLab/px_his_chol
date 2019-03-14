#make plots comparing replication of -log(P) and beta b/w HCHS and MESA HIS and MESA CAU at various P thresholds
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd("/home/angela/px_his_chol/MESA_compare/")

data <- fread("log_P.csv")
data$HCHS_P_sig_thres <- as.factor(data$HCHS_P_sig_thres)
data_log_P <- data %>% dplyr::select(HCHS_P_sig_thres, MESA_HIS_log_P, MESA_CAU_log_P) 
melted <- melt(data_log_P, id.vars = "HCHS_P_sig_thres")

log_P <- ggplot(data = melted, aes(x = HCHS_P_sig_thres, y = value, group = variable)) + 
  geom_line(aes(color = variable)) +
  scale_x_discrete(limits = rev(levels(melted$HCHS_P_sig_thres))) + 
  scale_color_brewer(palette = "Set1", labels = c("HIS", "CAU")) + 
  labs(title = "Correlation of HCHS and MESA -log10(P) values", color = "Pop.", y = "Correlation of -log10(P)", x = "HCHS P < #") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))
log_P

data_beta <- data %>% dplyr::select(HCHS_P_sig_thres, MESA_HIS_beta, MESA_CAU_beta) 
melted <- melt(data_beta, id.vars = "HCHS_P_sig_thres")

beta <- ggplot(data = melted, aes(x = HCHS_P_sig_thres, y = value, group = variable)) + 
  geom_line(aes(color = variable)) +
  scale_x_discrete(limits = rev(levels(melted$HCHS_P_sig_thres))) + 
  scale_color_brewer(palette = "Set1", labels = c("HIS", "CAU")) + 
  labs(title = "Correlation of HCHS and MESA beta values", color = "Pop.", y = "Correlation of betas", x = "HCHS P < #") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5))
beta
