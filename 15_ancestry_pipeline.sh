#Personal version of https://github.com/armartin/ancestry_pipeline
#Using MOSAIC instead of RFMix
#RUN ALL CHRS
for i in {1..22};
do

##NAT
  #1. Extract chr with PLINK
  cd /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/
  /usr/local/bin/plink --bfile 1000G_HCHS_geno_0.01_maf_0.05_ordered --keep NAT.txt --chr ${i} --make-bed --out NAT_chr${i}
  #2. Phase w/ HAPI-UR
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p NAT_chr${i} -w 64 -o NAT_chr${i}
  #3. Rename HAPI-UR output for MOSAIC input
  mv NAT_chr${i}.phgeno NATgenofile.${i}
  printf ":sites:" > rates.${i}
  wc -l < NATgenofile.${i} >> rates.${i} #add number of sites
  awk '{print $4}' NAT_chr${i}.phsnp | paste -sd ' ' >> rates.${i} #add row of positions
  awk '{print $3}' NAT_chr${i}.phsnp | paste -sd ' ' >> rates.${i} #add row of cM
  mv NAT_chr${i}.phsnp snpfile.${i}
  
##IBS
  #4. Extract chr with PLINK
  /usr/local/bin/plink --bfile 1000G_HCHS_geno_0.01_maf_0.05_ordered --keep IBS.txt --chr ${i} --make-bed --out IBS_chr${i}
  #5. Phase w/ HAPI-UR
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p IBS_chr${i} -w 64 -o IBS_chr${i}
  #6. Rename HAPI-UR output for MOSAIC input
  mv IBS_chr${i}.phgeno IBSgenofile.${i}
  
##YRI
  #7. Extract chr with PLINK
  /usr/local/bin/plink --bfile 1000G_HCHS_geno_0.01_maf_0.05_ordered --keep YRI.txt --chr ${i} --make-bed --out YRI_chr${i}
  #8. Phase w/ HAPI-UR
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p YRI_chr${i} -w 64 -o YRI_chr${i}
  #9. Rename HAPI-UR output for MOSAIC input
  mv YRI_chr${i}.phgeno YRIgenofile.${i}
  
##HCHS
  #10. Extract chr with PLINK
  /usr/local/bin/plink --bfile 1000G_HCHS_geno_0.01_maf_0.05_ordered --keep HCHS.txt --chr ${i} --make-bed --out HCHS_chr${i}
  
  #FOR TESTING PURPOSES KEEP SMALL SUBSET
  #/usr/local/bin/plink --bfile HCHS_chr${i} --thin-indiv-count 100 --make-bed --out HCHS_chr${i}
  
  #11. Phase w/ HAPI-UR
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p HCHS_chr${i} -w 64 -o HCHS_chr${i}
  #12. Rename HAPI-UR output for MOSAIC input
  mv HCHS_chr${i}.phgeno HCHSgenofile.${i}
  
##Run MOSAIC
  mkdir -p MOSAIC_RESULTS #why doesn't MOSAIC automatically make these
  mkdir -p MOSAIC_PLOTS
  /usr/bin/Rscript /home/angela/px_his_chol/MOSAIC/mosaic.R HCHS /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/sep_pops/ -a 3 -n 24000 -c ${i}:${i}
done

#how in the world do you read MOSAIC output
