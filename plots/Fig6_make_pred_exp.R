library(cowplot)
library(data.table)
library(ggplot2)
library(dplyr)
#HDL ICAM1 in ESPMSL, HCHS
HCHS_x <- fread('/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt', header = T)
HCHS_yGTEx <- fread('/home/angela/px_his_chol/PrediXcan/HIS.db_predicted_expression.txt', header = T)
HCHS_xyGTEx <- inner_join(HCHS_x, HCHS_yGTEx, by = "IID")
#pdf("/home/angela/px_his_chol/Manuscript_figures/HDL_ICAM1_ESPMSL_HCHS.pdf", width = 4.25, height = 4)
HCHS <- ggplot(HCHS_xyGTEx, aes(x = HCHS_xyGTEx$ENSG00000087237.6, y = HDL_rank)) + 
  geom_jitter(size = 0.75, color = "#b3cde3") + 
  geom_density_2d(color = "#8c96c6") + 
  stat_smooth(method="lm", se = TRUE, fullrange = TRUE, color = "#8856a7") + 
  scale_x_continuous(name = "Predicted gene expression") + 
  scale_y_continuous(name = "HDL (rank normalized)") + 
  theme_bw() + 
  theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5)) +
  ggtitle("HCHS/SoL")
#dev.off()

#HDL ICAM1 in ESPMSL, MESA HIS
MESA_HIS_x <- fread('/home/angela/px_his_chol/MESAreplication/MESA_HIS_pheno.txt')
MESA_HIS_yGTEx <- fread('/home/angela/px_his_chol/MESAreplication/PrediXcan/HIS.db_predicted_expression.txt')
MESA_HIS_xyGTEx <- inner_join(MESA_HIS_x, MESA_HIS_yGTEx, by = "IID")

#pdf("/home/angela/px_his_chol/Manuscript_figures/HDL_ICAM1_ESPMSL_MESA_HIS.pdf", width = 4.25, height = 4)
MESA_HIS <- ggplot(MESA_HIS_xyGTEx, aes(x = MESA_HIS_xyGTEx$ENSG00000087237.6, y = HDL_rank)) + 
  geom_jitter(size = 0.75, color = "#b3cde3") + 
  geom_density_2d(color = "#8c96c6") + 
  stat_smooth(method="lm", se = TRUE, fullrange = TRUE, color = "#8856a7") + 
  scale_x_continuous(name = "Predicted gene expression") + 
  scale_y_continuous(name = "HDL (rank normalized)") + 
  theme_bw() + 
  theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5))+
  ggtitle("MESA HIS") 
#dev.off()

pdf("/home/angela/px_his_chol/Manuscript_figures/Fig6.pdf", width = 5.5, height = 3)
plot_grid(HCHS, MESA_HIS, labels = c('A', 'B'), scale = c(.9, .9, .9, .9))
dev.off()

tiff("/home/angela/px_his_chol/Manuscript_figures/Fig6.tiff", width = 13.97, height = 7.62, units = 'cm', res = 300, compression = 'lzw')
plot_grid(HCHS, MESA_HIS, labels = c('A', 'B'), scale = c(.9, .9, .9, .9))
dev.off()