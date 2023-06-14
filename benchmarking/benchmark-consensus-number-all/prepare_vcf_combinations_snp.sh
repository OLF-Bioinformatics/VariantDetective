cd $1 

vcf-isec -p snp_combination_ snp_indel/freebayes/freebayes.filt.vcf.gz snp_indel/haplotypecaller/haplotypecaller.filt.vcf.gz snp_indel/clair3/clair3.filt.vcf.gz

tabix -p vcf snp_combination_0_1_2.vcf.gz
tabix -p vcf snp_combination_0_1.vcf.gz
tabix -p vcf snp_combination_0_2.vcf.gz
tabix -p vcf snp_combination_1_2.vcf.gz
tabix -p vcf snp_combination_0.vcf.gz
tabix -p vcf snp_combination_1.vcf.gz
tabix -p vcf snp_combination_2.vcf.gz

bcftools concat -a snp_combination_0_1_2.vcf.gz \
snp_combination_0_1.vcf.gz \
snp_combination_0_2.vcf.gz \
snp_combination_1_2.vcf.gz \
snp_combination_0.vcf.gz \
snp_combination_1.vcf.gz \
snp_combination_2.vcf.gz \
-o snp_set1.vcf

bcftools concat -a snp_combination_0_1_2.vcf.gz \
snp_combination_0_1.vcf.gz \
snp_combination_0_2.vcf.gz \
snp_combination_1_2.vcf.gz \
-o snp_set2.vcf

bcftools concat -a snp_combination_0_1_2.vcf.gz \
-o snp_set3.vcf


