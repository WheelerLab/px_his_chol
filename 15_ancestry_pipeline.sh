#Personal version of https://github.com/armartin/ancestry_pipeline

#RUNNING THE SAMPLE TOY DATA
cd /home/angela/px_his_chol/ancestry_pipeline/toy_data/
#1. Split into chromosomes
for i in {1..22};
do 
  plink --bfile ACB_example --chr ${i} --make-bed --out ACB_example_chr${i};
done

#2. SHAPEIT2 check
for i in {1..22}; 
do 
  shapeit -check --input-ref /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.hap.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.legend.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3.sample -B ACB_example_chr${i} --input-map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --output-log ACB_example_chr${i}.mendel
done

#3. SHAPEIT2 phasing
for i in {1..22}
do
  shapeit --input-ref /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.hap.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3_chr${i}.legend.gz /home/angela/1000GP_Phase3_combined/1000GP_Phase3.sample -B ACB_example_chr${i} --duohmm --input-map /home/angela/1000GP_Phase3_combined/genetic_map_chr${i}_combined_b37.txt --exclude-snp ACB_example_chr${i}.mendel.snp.strand.exclude --output-max ACB_example_chr${i}.haps.gz ACB_example_chr${i}.sample
done










#RUNNING ACTUAL DATA



