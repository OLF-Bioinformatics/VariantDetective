nucmername=$(ls *nucmer.vcf)
parsnpname=$(ls *parsnp.vcf)
#cut --complement -f 10 *.varcall_relative_to_representative.parsnp.vcf > tmp.vcf
#mv tmp.vcf $parsnpname
#echo '##fileformat=VCFv4.2' | cat - $parsnpname > temp && mv temp $parsnpname
#sed -i 's/v4.1/v4.2/' $nucmername
#sed -i 's/sample1/SAMPLE/' $nucmername
#sed -i -E 's/(#C.+\s)(.+)$/\1SAMPLE/' $parsnpname

bgzip -c -@ 40 $parsnpname > $parsnpname.gz
tabix -p vcf $parsnpname.gz
bgzip -c -@ 40 $nucmername > $nucmername.gz
tabix -p vcf $nucmername.gz
vcf-isec -p truecalls_ $parsnpname.gz $nucmername.gz
tabix -p vcf truecalls_0_1.vcf.gz
bgzip -c -@ 40 snp_indel/snp_final.vcf > snp_indel/snp_final.vcf.gz
tabix -p vcf snp_indel/snp_final.vcf.gz
vcf-isec -p final_comparison_ truecalls_0_1.vcf.gz snp_indel/snp_final.vcf.gz
vcf-isec -p freebayes_ truecalls_0_1.vcf.gz snp_indel/freebayes/freebayes.filt.vcf.gz
vcf-isec -p hc_ truecalls_0_1.vcf.gz snp_indel/haplotypecaller/haplotypecaller.filt.vcf.gz
vcf-isec -p clair3_ truecalls_0_1.vcf.gz snp_indel/clair3/clair3.filt.vcf.gz

zcat final_comparison_0_1.vcf.gz | grep -v '#' | wc -l
zcat final_comparison_1.vcf.gz | grep -v '#' | wc -l
zcat final_comparison_0.vcf.gz | grep -v '#' | wc -l
zcat freebayes_0_1.vcf.gz | grep -v '#' | wc -l
zcat freebayes_1.vcf.gz | grep -v '#' | wc -l
zcat freebayes_0.vcf.gz | grep -v '#' | wc -l
zcat hc_0_1.vcf.gz | grep -v '#' | wc -l
zcat hc_1.vcf.gz | grep -v '#' | wc -l
zcat hc_0.vcf.gz | grep -v '#' | wc -l
zcat clair3_0_1.vcf.gz | grep -v '#' | wc -l
zcat clair3_1.vcf.gz | grep -v '#' | wc -l
zcat clair3_0.vcf.gz | grep -v '#' | wc -l
