#Pipeline to merge and convert all downloaded genotypes from SGDP and INMEGEN

#SGDP
cd /home/angela/px_his_chol/SGDP_filtered/ 
declare -a arr=("LP6005441-DNA_A04" "LP6005441-DNA_A12" "LP6005441-DNA_B04" "LP6005441-DNA_B12" "LP6005441-DNA_E10" "LP6005441-DNA_F10" "LP6005441-DNA_G06" "LP6005441-DNA_G07" "LP6005441-DNA_H06" "LP6005441-DNA_H07" "LP6005443-DNA_A12" "LP6005443-DNA_E11" "LP6005443-DNA_F11" "LP6005443-DNA_G11" "LP6005443-DNA_H11" "LP6005519-DNA_D01" "LP6005519-DNA_G02" "LP6005677-DNA_D01" "LP6005677-DNA_E01" "LP6005677-DNA_F01" "SS6004476" "SS6004479") 

for ind in "${arr[@]}"
do
  echo "Starting filtering on" $ind
  #1. Filter by FL = 9
    #"filter level in range 0-9 or no value (non-integer: N,?) with zero being least reliable; to threshold at FL=n, use all levels n-9"
  bcftools filter --threads 40 -i 'FORMAT/FL=="9"' --output $ind.FL9.vcf /home/wheelerlab3/Data/SGDP/$ind.annotated.nh.vcf.gz
  bgzip $ind.FL9.vcf
  #2. Remove non-genotype information
  bcftools annotate --threads 40 -x INFO,^FORMAT/GT $ind.FL9.vcf.gz > $ind.FL9.noINFO.vcf
  bgzip $ind.FL9.noINFO.vcf
  echo "Done filtering on" $ind
  bcftools index --threads 40 -f --tbi -o $ind.FL9.noINFO.vcf.tbi $ind.FL9.noINFO.vcf.gz
  echo "Done indexing on" $ind
done

#3. Merge samples
bcftools merge --threads 40 --merge all LP*.FL9.noINFO.vcf.gz SS*.FL9.noINFO.vcf.gz > SGDP.FL9.noINFO.vcf
bgzip SGDP.FL9.noINFO.vcf
bcftools index --threads 40 -f --tbi SGDP.FL9.noINFO.vcf.gz > SGDP.FL9.noINFO.vcf.tbi
echo "SGDP samples have been merged"

#NATIVE MEXICAN
cd /home/angela/px_his_chol/INMEGEN/
declare -a arr=("MAY1" "MAY2" "NAH1" "NAH2" "TAR1" "TAR2" "TEP1" "TEP2" "TOT1" "TOT2" "ZAP1" "ZAP2")

for ind in "${arr[@]}"
do
  #4. Filter by passing SNPs
  echo "Starting filtering on" $ind
  bcftools filter -i  --threads 40'FORMAT/FT=="PASS"' /home/wheelerlab3/Data/INMEGEN/12G/$ind.vcf.gz > $ind.PASS.vcf
  bgzip $ind.PASS.vcf
  #5. Remove non-genotype information
  bcftools annotate --threads 40 -x INFO,^FORMAT/GT $ind.PASS.vcf.gz > $ind.PASS.noINFO.vcf
  bgzip $ind.PASS.noINFO.vcf
  echo "Done filtering on" $ind
  bcftools index --threads 40 -f --tbi $ind.PASS.noINFO.vcf.gz > $ind.PASS.noINFO.vcf.tbi
  echo "Done indexing on" $ind
done

#6. Merge samples
bcftools merge --threads 40 --merge all *.PASS.noINFO.vcf.gz > NativeMexican.PASS.noINFO.vcf
bgzip NativeMexican.PASS.noINFO.vcf
bcftools index --threads 40 -f --tbi NativeMexican.PASS.noINFO.vcf.gz > NativeMexican.PASS.noINFO.vcf.tbi
echo "Native Mexican samples have been merged"

#7. Merge SGDP and Native Mexican
cd /home/angela/px_his_chol/SGDP_filtered/
bcftools merge --threads 40 SGDP.FL9.noINFO.vcf.gz /home/angela/px_his_chol/INMEGEN/NativeMexican.PASS.noINFO.vcf.gz > NativeAmerican.PASS.noINFO.vcf
bgzip NativeAmerican.PASS.noINFO.vcf

#8. Add rsid, REF, and ALT
bcftools annotate --threads 40 -a /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423.vcf.gz -c CHROM,POS,ID NativeAmerican.PASS.noINFO.vcf.gz > NativeAmerican.PASS.noINFO.rsid.vcf
bgzip NativeAmerican.PASS.noINFO.rsid.vcf
bcftools index --threads 40 -f --tbi NativeAmerican.PASS.noINFO.rsid.vcf.gz > NativeAmerican.PASS.noINFO.rsid.vcf.tbi
echo "rsids have been added"

#9. Filter to SNPs only
bcftools view --threads 40 --known -m2 -M2 -v snps NativeAmerican.PASS.noINFO.rsid.vcf.gz > NativeAmerican.PASS.noINFO.rsid.SNPs.vcf
bgzip NativeAmerican.PASS.noINFO.rsid.SNPs.vcf
bcftools index --threads 40 -f --tbi NativeAmerican.PASS.noINFO.rsid.SNPs.vcf.gz > NativeAmerican.PASS.noINFO.rsid.SNPs.vcf.tbi
echo "Samples have been filtered to SNPs only"

#10. Convert to PLINK
plink --vcf NativeAmerican.PASS.noINFO.rsid.SNPs.vcf.gz --vcf-half-call h --make-bed --out NativeAmerican-h
  #half-called genotypes treated as haploid/homozygous in Mexican (Complete Genomics data)




