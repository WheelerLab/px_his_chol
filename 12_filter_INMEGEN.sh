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
  bgzip $ind.PASS.anno.vcf
  echo "Done annotating on" $ind
  #3. SNPs only
  bcftools view --threads 40 --known -m2 -M2 -v snps $ind.PASS.anno.vcf.gz > $ind.PASS.anno.known.vcf
  bgzip $ind.PASS.anno.known.vcf
  bcftools index --threads 40 -f --tbi $ind.PASS.anno.known.vcf.gz > $ind.PASS.anno.known.vcf.tbi
  echo "Done viewing on" $ind
done

#then merge
bcftools merge --merge all MAY1.PASS.anno.known.vcf.gz MAY2.PASS.anno.known.vcf.gz NAH1.PASS.anno.known.vcf.gz NAH2.PASS.anno.known.vcf.gz TAR1.PASS.anno.known.vcf.gz TAR2.PASS.anno.known.vcf.gz TOT1.PASS.anno.known.vcf.gz TOT2.PASS.anno.known.vcf.gz ZAP1.PASS.anno.known.vcf.gz ZAP2.PASS.anno.known.vcf.gz > NativeMexican.vcf
bgzip NativeMexican.vcf
bcftools index --threads 40 -f --tbi NativeMexican.vcf.gz > NativeMexican.vcf.tbi
plink --vcf NativeMexican.vcf.gz --vcf-half-call h --make-bed --out NativeMexican-h
  #half-called genotypes treated as haploid/homozygous (Complete Genomics data)
