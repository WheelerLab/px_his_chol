#Convert and filter individual vcf files to PLINK format
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
done

#2. Merge samples
bcftools merge --merge all MAY1.PASS.vcf.gz MAY2.PASS.vcf.gz NAH1.PASS.vcf.gz NAH2.PASS.vcf.gz TAR1.PASS.vcf.gz TAR2.PASS.vcf.gz TOT1.PASS.vcf.gz TOT2.PASS.vcf.gz ZAP1.PASS.vcf.gz ZAP2.PASS.vcf.gz > NativeMexican.vcf
bgzip NativeMexican.vcf
bcftools index --threads 40 -f --tbi NativeMexican.PASS.vcf.gz > NativeMexican.PASS.vcf.tbi
echo "Samples have been merged"

#3. Add rsid, REF, and ALT
bcftools annotate --threads 40 -a /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423.vcf.gz -c CHROM,POS,ID NativeMexican.PASS.vcf.gz > NativeMexican.PASS.rsid.vcf
  #add rsid
bgzip NativeMexican.PASS.rsid.vcf
bcftools index --threads 40 -f --tbi NativeMexican.PASS.rsid.vcf.gz > NativeMexican.PASS.rsid.vcf.tbi
echo "rsids have been added"
bcftools +fixref NativeMexican.PASS.rsid.vcf.gz -Ob -o NativeMexican.PASS.match.vcf --threads 40 -- -d -f /home/angela/human_g1k_v37.fasta -i /home/angela/px_his_chol/SGDP_filtered/anno/All_20180423.vcf.gz
  #add REF/ALT
bcftools sort NativeMexican.PASS.match.vcf -Ob -o NativeMexican.PASS.match.sorted.vcf
echo "Samples have been sorted"

#4. Filter to SNPs only
bcftools view --threads 40 --known -m2 -M2 -v snps NativeMexican.PASS.match.sorted.vcf > NativeMexican.PASS.match.sorted.known.vcf
bgzip NativeMexican.PASS.match.sorted.vcf
bgzip NativeMexican.PASS.match.sorted.known.vcf
bcftools index --threads 40 -f --tbi NativeMexican.PASS.match.sorted.known.vcf.gz > NativeMexican.PASS.match.sorted.known.vcf.tbi
echo "Samples have been filtered to SNPs only"

#5. Convert to PLINK
plink --vcf NativeMexican.PASS.match.sorted.known.vcf.gz --vcf-half-call h --make-bed --out NativeMexican-h
  #half-called genotypes treated as haploid/homozygous (Complete Genomics data)

