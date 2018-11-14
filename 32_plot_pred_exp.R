#make scatterplot of predicted expression vs. phenotype
library(data.table)
library(ggplot2)
library(dplyr)

#HDL CETP in ARTCRN, HCHS
x <- fread('/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt', header = T)
yGTEx <- fread('/home/angela/px_his_chol/PrediXcan/TW_Artery_Coronary_0.5.db_predicted_expression.txt', header = T)
xyGTEx <- inner_join(x, yGTEx, by = "IID")
pdf("/home/angela/px_his_chol/Manuscript_figures/HDL_CETP_ARTCRN_HCHS.pdf", width = 4.25, height = 4)
ggplot(xyGTEx, aes(x = xyGTEx$ENSG00000087237.6, y = HDL_rank)) + 
  geom_jitter(size = 0.75, color = "#06cdd2ff") + 
  stat_smooth(method="lm", se = TRUE, color = "#999999ff", fullrange = TRUE) + 
  scale_x_continuous(name = "Predicted gene expression") + 
  scale_y_continuous(name = "HDL (rank normalized)") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) #+ 
#ggtitle("CETP in ARTCRN, HCHS")
dev.off()

#HDL CETP in ARTCRN, MESA HIS
x <- fread('/home/angela/px_his_chol/MESAreplication/MESA_HIS_pheno.txt')
yGTEx <- fread('/home/angela/px_his_chol/MESAreplication/PrediXcan/TW_Artery_Coronary_0.5.db_predicted_expression.txt')
xyGTEx <- inner_join(x, yGTEx, by = "IID")

pdf("/home/angela/px_his_chol/Manuscript_figures/HDL_CETP_ARTCRN_MESA_HIS.pdf", width = 4.25, height = 4)
ggplot(xyGTEx, aes(x = xyGTEx$ENSG00000087237.6, y = HDL_rank)) + 
  geom_jitter(size = 0.75, color = "#06cdd2ff") + 
  stat_smooth(method="lm", se = TRUE, color = "#999999ff", fullrange = TRUE) + 
  scale_x_continuous(name = "Predicted gene expression") + 
  scale_y_continuous(name = "HDL (rank normalized)") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) #+ 
#ggtitle("CETP in ARTCRN, MESA HIS")
dev.off()

#CHOL CCL22 in ESPMSL, HCHS
x <- fread('/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt', header = T)
yGTEx <- fread('/home/angela/px_his_chol/PrediXcan/TW_Esophagus_Muscularis_0.5.db_predicted_expression.txt', header = T)
xyGTEx <- inner_join(x, yGTEx, by = "IID")
pdf("/home/angela/px_his_chol/Manuscript_figures/CHOL_CCL22_ESPMSL_HCHS.pdf", width = 4.25, height = 4)
ggplot(xyGTEx, aes(x = xyGTEx$ENSG00000102962.4, y = CHOL_rank)) + 
  geom_jitter(size = 0.75, color = "#06cdd2ff") + 
  stat_smooth(method="lm", se = TRUE, color = "#999999ff", fullrange = TRUE) + 
  scale_x_continuous(name = "Predicted gene expression") + 
  scale_y_continuous(name = "CHOL (rank normalized)") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) #+ 
dev.off()

#CHOL CCL22 in ESPMSL, MESA HIS
x <- fread('/home/angela/px_his_chol/MESAreplication/MESA_HIS_pheno.txt')
yGTEx <- fread('/home/angela/px_his_chol/MESAreplication/PrediXcan/TW_Esophagus_Muscularis_0.5.db_predicted_expression.txt')
xyGTEx <- inner_join(x, yGTEx, by = "IID")

pdf("/home/angela/px_his_chol/Manuscript_figures/CHOL_CCL22_ESPMSL_MESA_HIS.pdf", width = 4.25, height = 4)
ggplot(xyGTEx, aes(x = xyGTEx$ENSG00000102962.4, y = CHOL_rank)) + 
  geom_jitter(size = 0.75, color = "#06cdd2ff") + 
  stat_smooth(method="lm", se = TRUE, color = "#999999ff", fullrange = TRUE) + 
  scale_x_continuous(name = "Predicted gene expression") + 
  scale_y_continuous(name = "CHOL (rank normalized)") + 
  theme_bw() + 
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) #+ 
dev.off()
