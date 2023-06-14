dir='benchmark-real_data-bwa-snp'
mkdir -p $dir
for i in {cft073,mgh78578,rbhstw00029,rbhstw00053,rbhstw00059,rbhstw00122,rbhstw00127,rbhstw00128,rbhstw00131,rbhstw00167,rbhstw00189,rbhstw00277,rbhstw00309,rbhstw00340,rbhstw00350,rhb10c07,rhb11c04,rhb14c01}
do
    ./variantdetective.py snp_indel \
    	-1 $dir/${i}/${i}.1.fq.gz \
        -2 $dir/${i}/${i}.2.fq.gz \
    	-r $dir/${i}/*.fa \
    	-t 40 \
    	-o $dir/${i}/
     rm $dir/${i}/*/*bam*
done
