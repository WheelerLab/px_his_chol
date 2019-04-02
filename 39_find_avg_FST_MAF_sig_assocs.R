#find the average FST (MESA_HIS, MESA_CAU) for each gene in the similar, not similar, and NA categories
library(data.table)
library(dplyr)
setwd("/home/angela/px_his_chol/MESA_compare/compare_NA_genes/FST/")

primary_sig_assocs_FST_MAF <- fread("primary_sig_assocs_FST_MAF.csv")
primary_sigs <- fread("/home/angela/px_his_chol/MESA_compare/HCHS_primary_sig.csv")
primary_sigs <- primary_sigs %>% dplyr::select(genename, model, pheno)
primary_sigs$avg_FST_MESA_HIS <- -1000
primary_sigs$avg_FST_MESA_CAU <- -1000
primary_sigs$avg_MAF_MESA_HIS <- -1000
primary_sigs$avg_MAF_MESA_CAU <- -1000

#add category to primary_sigs
categories <- primary_sig_assocs_FST_MAF %>% dplyr::select(genename, model, pheno, category)
primary_sigs <- left_join(primary_sigs, categories, by = c("genename", "model", "pheno"))
primary_sigs <- primary_sigs[!duplicated(primary_sigs), ]

for(i in 1:nrow(primary_sigs)){ #find the average FST in MESA HIS and MESA CAU for each model
  gn <- as.character(primary_sigs[i, 1]) #genename
  tiss <- as.character(primary_sigs[i, 2]) #tissue
  phe <- as.character(primary_sigs[i, 3]) #pheno
  
  assoc <- subset(primary_sig_assocs_FST_MAF, genename == gn & phe == pheno & model == tiss)
  avg_FST_MESA_HIS <- sum(abs(assoc$HCHS_MESA_HIS_FST), na.rm = T) / sum(!is.na(assoc$HCHS_MESA_HIS_FST))
  avg_FST_MESA_CAU <- sum(abs(assoc$HCHS_MESA_CAU_FST), na.rm = T) / sum(!is.na(assoc$HCHS_MESA_CAU_FST))
  avg_MAF_MESA_HIS <- sum(abs(assoc$HCHS_MAF - assoc$MESA_HIS_MAF), na.rm = T) / sum(!is.na(assoc$HCHS_MAF - assoc$MESA_HIS_MAF))
  avg_MAF_MESA_CAU <- sum(abs(assoc$HCHS_MAF - assoc$MESA_CAU_MAF), na.rm = T) / sum(!is.na(assoc$HCHS_MAF - assoc$MESA_CAU_MAF))
  
  primary_sigs[i, 4] <- avg_FST_MESA_HIS
  primary_sigs[i, 5] <- avg_FST_MESA_CAU
  primary_sigs[i, 6] <- avg_MAF_MESA_HIS
  primary_sigs[i, 7] <- avg_MAF_MESA_CAU
}

#split up models into categories of sig in all reps, not in EUR, and NA in EUR
similar <- subset(primary_sigs, category == "similar")
dissonant <- subset(primary_sigs, category == "dissonant")
EUR_both_NA <- subset(primary_sigs, category == "EUR_both_NA")
dissonant_EUR_both_NA <- subset(primary_sigs, category == "EUR_both_NA" | category == "dissonant")

#mean(dissonant$avg_FST_MESA_HIS)
#mean(EUR_both_NA$avg_FST_MESA_HIS)
mean(similar$avg_FST_MESA_HIS)
mean(dissonant_EUR_both_NA$avg_FST_MESA_HIS)
mean(similar$avg_MAF_MESA_HIS)
mean(dissonant_EUR_both_NA$avg_MAF_MESA_HIS)

#mean(dissonant$avg_FST_MESA_CAU)
#mean(EUR_both_NA$avg_FST_MESA_CAU)
mean(similar$avg_FST_MESA_CAU)
mean(dissonant_EUR_both_NA$avg_FST_MESA_CAU)
mean(similar$avg_MAF_MESA_CAU)
mean(dissonant_EUR_both_NA$avg_MAF_MESA_CAU)


