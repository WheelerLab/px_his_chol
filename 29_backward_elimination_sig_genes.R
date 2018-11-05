#backward elimination of significant genes to determine which ones are independent
library(dplyr)

pheno <- fread('/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt', header = T)
pheno <- pheno %>% select(FID, IID, CHOL_rank, HDL_rank, TRIG_rank, LDL_rank)
yGTEx <- fread('/home/angela/px_his_chol/PrediXcan/TW_Artery_Coronary_0.5.db_predicted_expression.txt', header = T)
