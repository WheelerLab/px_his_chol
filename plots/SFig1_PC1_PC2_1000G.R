#make PC1 vs. PC2 merged w/ 1000G, colored by region
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)

#regions
HCHS_region <- fread("/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt", header = T)
HCHS_region <- HCHS_region %>% dplyr::select("IID", "REGION")
HCHS_region$REGION <- as.character(HCHS_region$REGION)
HCHS_region$Cohort <- "HCHS/SoL"
for(ind in 1:nrow(HCHS_region)){ #Translate region numbers to actual region names
  if(is.na(HCHS_region[ind, 2])){
    next
  }else if(as.character(HCHS_region[ind, 2]) == "1"){
    HCHS_region[ind, 2] = "Central American"
  }else if(as.character(HCHS_region[ind, 2]) == "2"){
    HCHS_region[ind, 2] = "Cuban"
  }else if(as.character(HCHS_region[ind, 2]) == "3"){
    HCHS_region[ind, 2] = "Dominican"
  }else if(as.character(HCHS_region[ind, 2]) == "4"){
    HCHS_region[ind, 2] = "Mexican"
  }else if(as.character(HCHS_region[ind, 2]) == "5"){
    HCHS_region[ind, 2] = "Puerto Rican"
  }else if(as.character(HCHS_region[ind, 2]) == "6"){
    HCHS_region[ind, 2] = "South American"
  }
}
X1000G_region <- fread("/home/angela/1000G/fams/1000G_pops_fam.csv")
X1000G_region <- X1000G_region %>% dplyr::select("Sample", "Population")
colnames(X1000G_region) <- c("IID", "REGION")
X1000G_region$Cohort <- "1000G"
region <- rbind(HCHS_region, X1000G_region)

#pcs
kingpc <- fread("/home/angela/px_his_chol/1000G_PCs/kingpc.ped")
kingpc <- kingpc %>% dplyr::select(V1, V7, V8)
colnames(kingpc) <- c("IID", "PC1", "PC2")

#pcs & regions
pcs_region <- left_join(kingpc, region, by = "IID")
pcs_region$REGION <- gsub("IBS", "EUR", pcs_region$REGION)
pcs_region$REGION <- gsub("CEU", "EUR", pcs_region$REGION)
pcs_region$REGION <- gsub("MXL", "NAT", pcs_region$REGION)
pcs_region$REGION <- gsub("PEL", "NAT", pcs_region$REGION)
#pcs_region$REGION <- gsub("YRI", "1000G YRI", pcs_region$REGION)
colnames(pcs_region)[4] <- "Region"
pcs_region <- pcs_region[complete.cases(pcs_region),]

X1000G_pcs_region <- subset(pcs_region, Cohort == "1000G")
HCHS_pcs_region <- subset(pcs_region, Cohort == "HCHS/SoL")
colnames(HCHS_pcs_region)[4] <- "region"

'
#plot
PC1_PC2_RdYlBu <- ggplot() + 
  geom_point(data = HCHS_pcs_region, aes(x = PC1, y = PC2, color = Region)) +
  facet_wrap(~ Region) + 
  geom_point(data = X1000G_pcs_region, aes(x = PC1, y = PC2, color = region)) +
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  scale_color_brewer(palette = "RdYlBu")
PC1_PC2_RdYlBu

pdf("/home/angela/px_his_chol/Manuscript_figures/SFig1_PC1_PC2_RdYl_Bu.pdf", width = 6, height = 4)
print(PC1_PC2_RdYl_Bu)
dev.off()

tiff("/home/angela/px_his_chol/Manuscript_figures/SFig1_PC1_PC2_RdYl_Bu.tiff", width = 15.24, height = 7.62, units = "cm", res = 300, compression = "lzw")
print(PC1_PC2_RdYl_Bu)
dev.off()
'
PC1_PC2_Inferno <- ggplot() + 
  geom_point(data = HCHS_pcs_region, aes(x = PC1, y = PC2), color = "grey") +
  facet_wrap(~ region) + 
  geom_point(data = X1000G_pcs_region, aes(x = PC1, y = PC2, color = Region)) +
  theme_bw() + 
  theme(text = element_text(size = 15)) +
  scale_color_viridis(discrete = T, option = "inferno", begin = 0.25, end = 0.9)
PC1_PC2_Inferno

pdf("/home/angela/px_his_chol/Manuscript_figures/SFig1.pdf", width = 9, height = 6)
print(PC1_PC2_Inferno)
dev.off()

tiff("/home/angela/px_his_chol/Manuscript_figures/SFig1.tiff", width = 15.24, height = 7.62, units = 'cm', res = 300, compression = 'lzw')
print(PC1_PC2_Inferno)
dev.off()
