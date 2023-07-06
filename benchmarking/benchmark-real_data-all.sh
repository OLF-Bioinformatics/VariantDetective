dir='benchmark-real_data-fasta'
mkdir -p $dir
#for i in {rbhstw00053,rbhstw00059,rbhstw00127,rbhstw00131,rbhstw00167,rbhstw00189,rbhstw00277,rbhstw00309,rbhstw00340,rbhstw00350,rhb10c07,rhb11c04}
for i in rbhstw00059
do
    ./variantdetective.py snp_indel \
    	-g $dir/${i}/${i}.fasta \
    	-r $dir/${i}/*.fa \
    	-t 40 \
    	-o $dir/${i}/ \
	--readcov 100x \
	--readlen 100000,2000 \
#     rm $dir/${i}/*/*bam*
done
