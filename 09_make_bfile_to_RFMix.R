#converts from PLINK --recode beagle output to RFMix input
library(data.table)
library(dplyr)
wd <- "/home/angela/px_his_chol/RFMix/RFMix_v1.5.4/"
setwd(wd)
  
for(i in 1:22){
  cleanPLINK = function(fromPLINK){ #cleans up output from PLINK to be suitable for RFMix
    fromPLINK$V1 <- NULL
    fromPLINK = fromPLINK[-c(2:3),] #remove IID and pheno row
    for(hap in 1:ncol(fromPLINK)){ #add _A and _B to haplotypes
      if(hap %% 2 == 0){ #if column is even
        fromPLINK[1, hap] = paste0((fromPLINK[[1, hap]]), "_A")
      }else{ #if column is odd
        fromPLINK[1, hap] = paste0((fromPLINK[[1, hap]]), "_B")
      }
    }
    return(fromPLINK)
  }
  
  makeHCHS <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/HCHS/HIS3 --chr ", i ," --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/HCHS/HCHS", sep = "")
  system(makeHCHS)
  HCHS <- fread(paste("HCHS/HCHS.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  HCHS <- cleanPLINK(HCHS)
  numHCHS <- ncol(HCHS) - 1
  print(paste("HCHS chromosome ", i, " is done processing.", sep = ""))

  makeIBS <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/IBS --chr ", i ," --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/refPanel/IBS", sep = "")
  system(makeIBS)
  IBS <- fread(paste("refPanel/IBS.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  IBS <- cleanPLINK(IBS)
  numIBS <- ncol(IBS) - 1
  print(paste("IBS chromosome ", i, " is done processing.", sep = ""))

  makePEL <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/PEL --chr ", i ," --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/refPanel/PEL", sep = "")
  system(makePEL)
  PEL <- fread(paste("refPanel/PEL.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  PEL <- cleanPLINK(PEL)
  numPEL <- ncol(PEL) - 1
  print(paste("PEL chromosome ", i, " is done processing.", sep = ""))
  
  makeYRI <- paste("plink --bfile /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/YRI --chr ", i ," --recode beagle --out /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/refPanel/YRI", sep = "")
  system(makeYRI)
  YRI <- fread(paste("refPanel/YRI.chr-", i, ".dat", sep = ""), header=F, stringsAsFactors=F, colClasses = "character")
  YRI <- cleanPLINK(YRI)
  numYRI <- ncol(YRI) - 1
  print(paste("YRI chromosome ", i, "is done processing. Starting processing of the allele matrix.", sep = ""))

  all <- left_join(HCHS, IBS, by = "V2") #join by SNPs in the admixed pop
  all <- left_join(all, PEL, by = "V2")
  all <- left_join(all, YRI, by = "V2")
  SNPs <- matrix(all$V2[2:ncols(all),], nrow=1)
  all$V2 <- NULL #remove SNP col
  IIDs <- as.character(all[1,])
  allT <- as.matrix(t(all)) #can't use transpose() for some unknown reason
  rownames(allT) <- IIDs
  allT <- allT[,-1]
  allT[allT == 0] <- NA #if 0 is in the matrix, causes error in asFactor

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
  write.table(alleles, paste(wd, "HCHS/alleles", i, ".txt", sep=""), quote=F, sep="", col.names=F, row.names=F)
  print(paste("Allele matrix for chromosome ", i, " is complete.", sep = ""))

  #make centimorgan map because this is apparently 1980
  HCHS_map <- fread(paste("HCHS/HCHS.chr-", i, ".map", sep = ""), header=F, stringsAsFactors=F)
  colnames(HCHS_map) <- c("rsid", "Position(bp)", "A1", "A0")
  UCSC_map <- fread(paste("calculatecM/genetic_map_GRCh37_chr", i, ".txt", sep = ""), header = T)
  map <- left_join(HCHS_map, UCSC_map)
  fwrite(as.data.frame(map[,7]), paste("HCHS/marker_positions_chr", i, ".txt", sep = ""), col.names = F, quote = F, na = "NA")
  print(paste("Marker position map for chromosome ", i, " is complete.", sep = ""))
}

# Classes
classes = c(rep(0,numHCHS), rep(1, numIBS), rep(2, numPEL), rep(3,numYRI))
classesM = matrix(classes, nrow=1)
write.table(classesM, paste(wd, "HCHS/classes.txt", sep = ""), quote=F, row.names=F, col.names=F)
print("Classes file is complete. Ending program.")
