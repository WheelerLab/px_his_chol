#make principal component eigenvalue plot as seen in PAGE fig. 4B
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
pcs <- fread("/home/angela/px_his_chol/KING/kingpc.ped")
region <- fread("/home/angela/px_his_chol/editedPheno/11_all_lipid_rank_12236_FID_IID.txt", header = T)
region <- region %>% dplyr::select("IID", "REGION")
pcs <- left_join(region, pcs, c("IID" = "V2"))
pcs$V1 <- NULL
pcs$V3 <- NULL
pcs$V4 <- NULL
pcs$V5 <- NULL
pcs$V6 <- NULL
pcs$FID <- NULL #Columns leftover from PLINK FID, irrelevant rn
colnames(pcs) <- c("IID", "REGION", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
pcs <- pcs[complete.cases(pcs),]
pcs <- pcs %>% dplyr::select(IID, REGION, PC1, PC2, PC3, PC4, PC5)
for(ind in 1:nrow(pcs)){ #Translate region numbers to actual region names
  if(as.character(pcs[ind,2]) == "1"){
    pcs[ind,2] = "Central American"
  }else if(as.character(pcs[ind,2]) == "2"){
    pcs[ind,2] = "Cuban"
  }else if(as.character(pcs[ind,2]) == "3"){
    pcs[ind,2] = "Dominican"
  }else if(as.character(pcs[ind,2]) == "4"){
    pcs[ind,2] = "Mexican"
  }else if(as.character(pcs[ind,2]) == "5"){
    pcs[ind,2] = "Puerto Rican"
  }else if(as.character(pcs[ind,2]) == "6"){
    pcs[ind,2] = "South American"
  }
}
pcs_t <- data.table::transpose(pcs) #why did this suddenly break
pcs_t <- pcs_t[-c(1:2),]
IDs <- pcs$IID
#PCs <- as.data.frame(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20"))
PCs <- as.data.frame(c("PC1", "PC2", "PC3", "PC4", "PC5"))
colnames(PCs) <- "PCs"
pcs_t <- cbind(PCs, pcs_t)
colnames(pcs_t) <- c("PCs", IDs)
melted <- melt(pcs_t, id.vars = "PCs")
#melted$PCs <- factor(melted$PCs, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20"))
melted$PCs <- factor(melted$PCs, levels = c("PC1", "PC2", "PC3", "PC4", "PC5"))
region <- pcs %>% dplyr::select(IID, REGION)
melted <- left_join(melted, region, by = c("variable" = "IID"))

#rotation plot
fig1 <- ggplot(data = melted, aes(x = PCs, y = as.numeric(value), group = variable)) + 
  geom_line(aes(color = REGION)) + 
  #scale_color_brewer(palette = "Set1") + 
  scale_color_viridis(discrete = T, option = "magma") + 
  labs(color = "Region", fill = "Region", y = "eigenvalue") + 
  theme_bw() + 
  guides(color = guide_legend(override.aes = list(size = 6))) + #make legend lines bigger
  theme(text = element_text(size = 15), legend.key = element_rect(colour = 'black', size = 0.2))
fig1

pdf("/home/angela/px_his_chol/Manuscript_figures/Fig1.pdf", width = 6, height = 3)
print(fig1)
dev.off()

tiff("/home/angela/px_his_chol/Manuscript_figures/Fig1.tiff", width = 15.24, height = 7.62, units = 'cm', res = 300, compression = 'lzw')
print(fig1)
dev.off()

#PC1 vs PC2
PC1_PC2 <- pcs %>% dplyr::select(REGION, PC1, PC2)
colnames(PC1_PC2) <- c("Region", "PC1", "PC2")
SFig1 <- ggplot() + 
  geom_point(data = PC1_PC2, aes(x = PC1, y = PC2, col = Region))+
  theme_bw() + 
  theme(text = element_text(size=15)) + 
  scale_colour_brewer(palette="Set1")

pdf("/home/angela/px_his_chol/Manuscript_figures/SFig1.pdf", width = 6, height = 3)
print(SFig1)
dev.off()

tiff("/home/angela/px_his_chol/Manuscript_figures/SFig1.tiff", width = 15.24, height = 7.62, units = 'cm', res = 300, compression = 'lzw')
print(SFig1)
dev.off()

