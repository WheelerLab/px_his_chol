#Personal version of https://github.com/armartin/ancestry_pipeline
#Phasing w/ HAPI-UR and local ancestry inference w/ RFMix
for i in {1..22};
do
  #1. Extract chr with PLINK
  cd /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/
  /usr/local/bin/plink --bfile 1000G_HCHS_geno_0.01_maf_0.05_ordered --chr ${i} --make-bed --out chr${i}

  #2. Phase w/ HAPI-UR
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p chr${i} -w 111 -o phase_chr${i}

  #3. Make additional files for RFMix input
  awk '{print $3}' phase_chr${i}.phsnp > phase_chr${i}.snp_locations
  /usr/bin/Rscript /home/angela/px_his_chol/ancestry_pipeline/make_classes_from_HAPI-UR.R /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/phase_chr${i}.phind /home/angela/px_his_chol/ancestry_pipeline/HCHS/ordered_pops.txt HCHS

  #4. Local ancestry inference w/ RFMix
  cd /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/ #requires being run within RFMix folder
  /home/angela/anaconda2/bin/python RunRFMix.py PopPhased /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/phase_chr${i}.phgeno /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/HCHS.classes /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/phase_chr${i}.snp_locations --output-name /home/angela/px_his_chol/local_anc_GEMMA/RFMix_output/chr${i} -n 5 --num-threads 15 --window-size 0.025
done
