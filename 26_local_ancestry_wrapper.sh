#Goes from the RFMix local ancestry estimates, PrediXcan-style dosages, and HAPI-UR phasings to full admixture mapping in GEMMA
#Uses 25_RFMix_loc_anc.py -> 21_make_GEMMA_input.R -> 02_PrediXcan_dosages_to_GEMMA.py -> 22_GEMMA_wrapper_admixture_mapping.py

echo "Converting RFMix output from chr. ${1}."
python 25_RFMix_loc_anc.py --admixed_hap haps.txt --phsnp /home/angela/px_his_chol/ancestry_pipeline/HCHS/no_NativeAmerican-h/PrediXcan_SNPs/phase_chr${1}.phsnp --Viterbi chr${1}.0.Viterbi.txt --output_prefix RFMix_for_GEMMA_${1}
echo "Completed converting RFMix output. Now making GEMMA input."
Rscript 21_make_GEMMA_input.R --ind_file_name RFMix_for_GEMMA_${1}_ind.txt --output ${1}
echo "Completed making GEMMA input. Now subsetting dosages to BIMBAM and anno."
python 02_PrediXcan_dosages_to_GEMMA.py --dosage_path /home/angela/px_his_chol/Imputation/UMich/UMich_dosages/ --local_anc_samples RFMix_for_GEMMA_${1}_ind.txt --local_anc_SNPs RFMix_for_GEMMA_${1}_SNPs.txt --chr ${1}
echo "Completed subsetting dosages. Now running admixture mapping in GEMMA."
python 22_GEMMA_wrapper_admixture_mapping.py --snptable RFMix_for_GEMMA_${1}.csv --relatedness relatedness_${1}.txt --BIMBAM BIMBAM/chr${1}.txt.gz --anno anno/anno${1}.txt --pheno pheno_${1}.txt --covariates covariates_${1}.txt --output chr${1}
echo "Completed admixture mapping. Have a nice day!"

