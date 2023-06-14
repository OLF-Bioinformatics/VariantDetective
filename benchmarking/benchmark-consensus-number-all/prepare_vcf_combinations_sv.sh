cd $1 

vcf-isec -p snp_combination_ snp_indel/freebayes/freebayes.filt.vcf.gz snp_indel/haplotypecaller/haplotypecaller.filt.vcf.gz snp_indel/clair3/clair3.filt.vcf.gz
ls structural_variant/nanovar/*pass.vcf \
   structural_variant/cutesv/variants.vcf \
   structural_variant/nanosv/variants.vcf \
   structural_variant/svim/variants.filt.vcf \
   > sv_vcf_list

for i in {1,2,3,4}
do
   SURVIVOR merge sv_vcf_list 1000 $i 1 1 0 50 combined_sv_$i.vcf
done
