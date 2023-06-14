dir='benchmark-real_data-sv'
mkdir -p $dir
#for i in rbhstw00189
for i in {cft073,mgh78578,rbhstw00029,rbhstw00053,rbhstw00059,rbhstw00122,rbhstw00127,rbhstw00128,rbhstw00131,rbhstw00167,rbhstw00189,rbhstw00277,rbhstw00309,rbhstw00340,rbhstw00350,rhb10c07,rhb11c04,rhb14c01}
do
    ./variantdetective.py structural_variant \
    	-l tmp2/${i}.fastq \
    	-r benchmark-real_data-snp/${i}/*.fa \
    	-t 40 \
    	-o $dir/${i}/
     rm $dir/${i}/*/*bam*
done
