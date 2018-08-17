#!/bin/bash
#PBS -N test_anc_pipeline_1000G_ref
#PBS -S /bin/bash
#PBS -l walltime=500:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=200gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR
#Personal version of https://github.com/armartin/ancestry_pipeline

#RUN ALL CHRS
for i in {1..22};
do
  #1. Extract chr with PLINK
  cd /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/
  /usr/local/bin/plink --bfile 1000G_HCHS_geno_0.01_maf_0.05_ordered --chr ${i} --make-bed --out chr${i}

  #2. Phase w/ HAPI-UR
  #rate-limiting step
  /home/angela/px_his_chol/HAPI-UR/hapi-ur-1.01/hapi-ur -p chr${i} -w 256 -o phase_chr${i}
  #Chromosome 22 (23754 variants and 12477 people) took 34552.859 seconds in HAPI-UR

  #3. Make additional files for RFMix input
  awk '{print $3}' phase_chr${i}.phsnp > phase_chr${i}.snp_locations
  /usr/bin/Rscript /home/angela/px_his_chol/ancestry_pipeline/make_classes_from_HAPI-UR.R /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/phase_chr${i}.phind /home/angela/px_his_chol/ancestry_pipeline/HCHS/ordered_pops.txt HCHS
done

for i in {1..22};
do
  #4. Run RFMix
  #rate-limiting step
  cd /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/
  python RunRFMix.py \
  -e 2 \
  -w 0.2 \
  --num-threads 50 \
  --use-reference-panels-in-EM \
  --forward-backward \
  PopPhased \
  /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/phase_chr${i}.phgeno \
  /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/HCHS.classes \
  /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/phase_chr${i}.snp_locations \
  -o /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/phase_chr${i}.rfmix
done

#will take many days, check up again when classes start
