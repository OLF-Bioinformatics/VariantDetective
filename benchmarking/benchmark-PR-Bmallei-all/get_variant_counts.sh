NUMBER=${1}
SAMPLE=sim_${NUMBER}_${2}X
cp vcf_header sim_snp_${NUMBER}_fixed.vcf
grep -v '#' sim_snp_${NUMBER}.vcf | awk -v OFS='\t' '{$2+=1}1' >> sim_snp_${NUMBER}_fixed.vcf
bgzip -c sim_snp_${NUMBER}_fixed.vcf > sim_snp_${NUMBER}_fixed.vcf.gz
tabix -p vcf sim_snp_${NUMBER}_fixed.vcf.gz
bgzip -c $SAMPLE/snp_indel/snp_final.vcf > $SAMPLE/snp_indel/snp_final.vcf.gz
tabix -p vcf $SAMPLE/snp_indel/snp_final.vcf.gz
vcf-isec -p $SAMPLE/variantdetective_ sim_snp_${NUMBER}_fixed.vcf.gz $SAMPLE/snp_indel/snp_final.vcf.gz
vcf-isec -p $SAMPLE/freebayes_ sim_snp_${NUMBER}_fixed.vcf.gz $SAMPLE/snp_indel/freebayes/freebayes.filt.vcf.gz 
vcf-isec -p $SAMPLE/hc_ sim_snp_${NUMBER}_fixed.vcf.gz $SAMPLE/snp_indel/haplotypecaller/haplotypecaller.filt.vcf.gz
vcf-isec -p $SAMPLE/clair3_ sim_snp_${NUMBER}_fixed.vcf.gz $SAMPLE/snp_indel/clair3/clair3.filt.vcf.gz

zcat sim_snp_${NUMBER}_fixed.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_indel/snp_final.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_indel/freebayes/freebayes.filt.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_indel/haplotypecaller/haplotypecaller.filt.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/snp_indel/clair3/clair3.filt.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/variantdetective_0_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/freebayes_0_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/hc_0_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/clair3_0_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/variantdetective_0.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/freebayes_0.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/hc_0.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/clair3_0.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/variantdetective_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/freebayes_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/hc_1.vcf.gz | grep -v '#' | wc -l
zcat $SAMPLE/clair3_1.vcf.gz | grep -v '#' | wc -l
