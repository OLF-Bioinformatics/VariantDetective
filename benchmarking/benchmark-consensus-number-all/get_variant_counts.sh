NUMBER=${1}
SAMPLE=sim_${NUMBER}_${2}X
cp vcf_header sim_snp_${NUMBER}_fixed.vcf
grep -v '#' sim_snp_${NUMBER}.vcf | awk -v OFS='\t' '{$2+=1}1' >> sim_snp_${NUMBER}_fixed.vcf
bgzip -c sim_snp_${NUMBER}_fixed.vcf > sim_snp_${NUMBER}_fixed.vcf.gz
tabix -p vcf sim_snp_${NUMBER}_fixed.vcf.gz
bgzip -c $SAMPLE/snp_set1.vcf > $SAMPLE/snp_set1.vcf.gz
tabix -p vcf $SAMPLE/snp_set1.vcf.gz
bgzip -c $SAMPLE/snp_set2.vcf > $SAMPLE/snp_set2.vcf.gz
tabix -p vcf $SAMPLE/snp_set2.vcf.gz
bgzip -c $SAMPLE/snp_set3.vcf > $SAMPLE/snp_set3.vcf.gz
tabix -p vcf $SAMPLE/snp_set3.vcf.gz

vcf-isec -p $SAMPLE/snp_set1_ sim_snp_${NUMBER}_fixed.vcf.gz $SAMPLE/snp_set1.vcf.gz
vcf-isec -p $SAMPLE/snp_set2_ sim_snp_${NUMBER}_fixed.vcf.gz $SAMPLE/snp_set2.vcf.gz
vcf-isec -p $SAMPLE/snp_set3_ sim_snp_${NUMBER}_fixed.vcf.gz $SAMPLE/snp_set3.vcf.gz

zcat sim_snp_${NUMBER}_fixed.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set2.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set3.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set1_0_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set2_0_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set3_0_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set1_0.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set2_0.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set3_0.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set1_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set2_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_set3_1.vcf.gz | grep -v '#' | wc -l
