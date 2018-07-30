#Personal version of https://github.com/armartin/ancestry_pipeline

#RUNNING THE SAMPLE TOY DATA
cd /home/angela/px_his_chol/ancestry_pipeline/toy_data/
for i in {1..22};
#for i in {22..22};
do 
  #Note: the original .bed/.bim/.fam already has combined reference/test
  #1. Split into chromosomes
  plink --bfile ACB_example --chr ${i} --make-bed --geno 0.01 --out ACB_example_chr${i};
  #2. Check data is aligned to reference
  shapeit -check --input-ref /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.hap.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.legend.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3.sample -B ACB_example_chr${i} --input-map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --output-log ACB_example_chr${i}.mendel
  #3. Phasing the dataset using the reference panel of haplotypes (long step)
  shapeit --input-ref /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.hap.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.legend.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3.sample -B ACB_example_chr${i} --duohmm --input-map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --exclude-snp ACB_example_chr${i}.mendel.snp.strand.exclude --output-max ACB_example_chr${i}.haps.gz ACB_example_chr${i}.sample
  #4. Convert to RFMix input
  python ../shapeit2rfmix.py --shapeit_hap_ref ACB_example_chr${i}.haps.gz --shapeit_hap_admixed ACB_example_chr${i}.haps.gz --shapeit_sample_ref ACB_example_chr${i}.sample --shapeit_sample_admixed ACB_example_chr${i}.sample --ref_keep ACB_example.ref --admixed_keep ACB_example.notref --chr ${i} --genetic_map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --out ACB_example
done

#5. Fix classes file for RFMix
python ../classes.py --ref CEU_example.keep,YRI_example.keep --sample CEU_YRI_ACB.sample --out CEU_YRI_ACB.classes
  #MUST KEEP CONSISTENT WITH INDVIDUAL ORDER IN THE ALLELES FILE

#6. Run RFMix (long step)
cd /home/angela/px_his_chol/RFMix/RFMix_v1.5.4/
for i in {1..22};
do
  python RunRFMix.py -e 2 -w 0.2 --num-threads 40 --use-reference-panels-in-EM --forward-backward PopPhased /home/angela/px_his_chol/ancestry_pipeline/toy_data/ACB_example_chr${i}.alleles /home/angela/px_his_chol/ancestry_pipeline/toy_data/CEU_YRI_ACB.classes /home/angela/px_his_chol/ancestry_pipeline/toy_data/ACB_example_chr${i}.snp_locations -o /home/angela/px_his_chol/ancestry_pipeline/toy_data/CEU_YRI_ACB_chr${i}.rfmix
done








#RUNNING ACTUAL DATA



