#Personal version of https://github.com/armartin/ancestry_pipeline

#RUN ALL CHRS
for i in {1..22};
do
  #1. Extract chr with PLINK
  cd /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/
  /usr/local/bin/plink --bfile 1000G_HCHS_geno_0.01_maf_0.05_ordered --chr ${i} --make-bed --out chr${i}
  #limit # of SNPs here?

  #2. Phase w/ HAPI-UR
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p chr${i} -w 256 -o phase_chr${i}

  #3. Make additional files for RFMix input
  awk '{print $3}' phase_chr${i}.phsnp > phase_chr${i}.snp_locations
  /usr/bin/Rscript /home/angela/px_his_chol/ancestry_pipeline/make_classes_from_HAPI-UR.R /home/angela/p

  #4. Run RFMix (takes very long time)
  cd /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/
  python RunRFMix.py -e 2 -w 0.2 --num-threads 40 --use-reference-panels-in-EM --forward-backward PopPha
done

#will take many days, check up again when classes start
