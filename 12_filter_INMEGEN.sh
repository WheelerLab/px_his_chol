declare -a arr=("MAY1" "MAY2" "NAH1" "NAH2" "TAR1" "TAR2" "TOT1" "TOT2" "ZAP1" "ZAP2")
cd /home/angela/px_his_chol/INMEGEN/

for ind in "${arr[@]}"
do
  #1. Filter by passing SNPs
  echo "Starting filtering on" $ind
  bcftools filter -i 'FORMAT/FT=="PASS"' /home/wheelerlab3/Data/INMEGEN/12G/$ind.vcf.gz > $ind.PASS.vcf
  bgzip $ind.PASS.vcf
  echo "Done filtering on" $ind
  bcftools index --threads 40 -f --tbi $ind.PASS.vcf.gz > $ind.PASS.vcf.tbi
  echo "Done indexing on" $ind
  #2. Add rsids
  bcftools annotate --threads 40 -a /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423.vcf.gz -c CHROM,POS,ID,ALT $ind.PASS.vcf.gz > $ind.PASS.anno.vcf
  echo "Done annotating on" $ind
  #3. SNPs only
  bcftools view --threads 40 --known -m2 -M2 -v snps MAY1.PASS.anno.vcf > $ind.PASS.anno.known.vcf
  bgzip $ind.PASS.anno.vcf
  bgzip $ind.PASS.anno.known.vcf
  echo "Done viewing on" $ind
done

#then merge
