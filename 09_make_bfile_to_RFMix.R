#converts from PLINK --recode beagle output to RFMix input
library(data.table)
library(dplyr)
wd <- "/home/angela/px_his_chol/RFMix/RFMix_v1.5.4/"
setwd(wd)
remove_SNP_threshold <- 0

for(i in 1:22){
  cleanPLINK = function(fromPLINK){ #cleans up output from PLINK to be suitable for RFMix
    fromPLINK$V1 <- NULL
    fromPLINK = fromPLINK[-c(2:3),] #remove pheno row
    #for(hap in 1:ncol(fromPLINK)){ #add _A and _B to haplotypes
    #  if(hap %% 2 == 0){ #if column is even
    #    fromPLINK[1, hap] = paste0((fromPLINK[[1, hap]]), "_A")
    #  }else{ #if column is odd
    #    fromPLINK[1, hap] = paste0((fromPLINK[[1, hap]]), "_B")
    #  }
    #}
    fromPLINK <- left_join(SNPs_in_adm, fromPLINK, by = "V2")
    fromPLINK[fromPLINK == 0] <- NA #replace all 0's w/ NA's
    #remove individuals w/ low call rates
    return(fromPLINK)
  }
  
  makeHCHS <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/10m/HIS_PX_SNPs --chr ", i ," --mind 0.99 --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/HCHS/HCHS", sep = "")
  system(makeHCHS)
  HCHS <- fread(paste("HCHS/HCHS.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  SNPs_in_adm <- HCHS$V2
  SNPs_in_adm <- SNPs_in_adm[-c(2:3)]
  SNPs_in_adm <- as.data.frame(SNPs_in_adm)
  colnames(SNPs_in_adm) <- "V2"
  HCHS <- cleanPLINK(HCHS)
  print(paste("HCHS, chromosome ", i, " has finished processing.", sep = ""))
  #numHCHS <- ncol(HCHS) - 1 #set after removing low call rate persons
  
  makeIBS <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/IBS --chr ", i ," --mind 0.99 --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/refPanel/IBS", sep = "")
  system(makeIBS)
  IBS <- fread(paste("refPanel/IBS.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  IBS <- cleanPLINK(IBS)
  missing <- colSums(is.na(IBS))
  ref_missing <- as.numeric(names(table(missing)[2])) #num of. missing haplotypes for 1000G
  numIBS <- ncol(IBS) - 1
  print(paste("IBS, chromosome ", i, " has finished processing.", sep = ""))
  
  #remove individuals in test pop w/ low call rate
  HCHS <- HCHS[, names(which(colSums(is.na(HCHS)) <= remove_SNP_threshold))] #set this threshold as necessary
  numHCHS <- ncol(HCHS) - 1
  IIDs <- unique(transpose(as.data.frame(HCHS[1,]))) #individuals remaining
  fwrite(IIDs, paste("/home/angela/px_his_chol/RFMix/RFMix_v1.5.4/IIDs_chr_", i, "_", remove_SNP_threshold, ".txt", sep = ""), col.names = F, row.names = F, quote = F, na = "NA")
  
  makePEL <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/PEL --chr ", i ," --mind 0.99 --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/refPanel/PEL", sep = "")
  system(makePEL)
  PEL <- fread(paste("refPanel/PEL.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  PEL <- cleanPLINK(PEL)
  numPEL <- ncol(PEL) - 1
  print(paste("PEL, chromosome ", i, " has finished processing.", sep = ""))
  
  makeYRI <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/YRI --chr ", i ," --mind 0.99 --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/refPanel/YRI", sep = "")
  system(makeYRI)
  YRI <- fread(paste("refPanel/YRI.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  YRI <- cleanPLINK(YRI)
  numYRI <- ncol(YRI) - 1
  print(paste("YRI, chromosome ", i, " has finished processing. Now joining all pops. by matching SNPs.", sep = ""))
  
  all <- left_join(HCHS, IBS, by = "V2")
  all <- left_join(all, PEL, by = "V2")
  all <- left_join(all, YRI, by = "V2")
  all <- all[-1,] #remove FID row
  
  #position data
  print("All populations have been added together. Now making marker positions file.")
  RU_map <- fread(paste("genome_map/rutgers_map_v3/RUMapv3_B137_chr", i, ".txt", sep = ""), header = T)
  RU_map <- RU_map %>% select(Marker_name, Sex_averaged_start_map_position)
  colnames(RU_map) <- c("V2", "cM")
  all <- left_join(all, RU_map, by = "V2")
  
  print("Removing all NA SNPs.")
  all_noNA <- all[complete.cases(all),] #only take SNPs in all people and that have cM
  print(paste(dim(all_noNA)[1], "SNPs remain in the dataset."))
  cMs <- all_noNA %>% select(V2, cM)
  fwrite(as.data.frame(cMs$cM), paste("HCHS/marker_positions_chr", i, "_", remove_SNP_threshold ,".txt", sep = ""), col.names = F, quote = F, na = "NA")
  all_noNA$V2 <- NULL #remove SNP col
  all_noNA$cM <- NULL
  allT <- as.matrix(t(all_noNA)) #can't use transpose() for some unknown reason
  
  binarize = function(snp){
    asFactor = factor(snp)
    allele1 = levels(asFactor)[1]
    allele2 = levels(asFactor)[2]
    result = rep(NA, length(snp))
    result[snp == allele1] = '0'
    result[snp == allele2] = '1'
    return(result)
  }
  
  all01 = apply(allT, 2, binarize)
  
  # Transpose back
  alleles = as.data.frame(t(all01))
  write.table(alleles, paste(wd, "HCHS/alleles", i, "_", remove_SNP_threshold ,".txt", sep=""), quote=F, sep="", col.names=F, row.names=F)
  print(paste("Chromosome ", i, " has fully processed.", sep = ""))
}

# Classes
classes = c(rep(0,numHCHS), rep(1, numIBS), rep(2, numPEL), rep(3,numYRI))
classesM = matrix(classes, nrow=1)
write.table(classesM, paste(wd, "HCHS/classes", i, "_", remove_SNP_threshold ,".txt", sep = ""), quote=F, row.names=F, col.names=F)
print("Classes file is complete. Ending program.")
