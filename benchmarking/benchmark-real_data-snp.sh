dir='benchmark-real_data-snp'
mkdir -p $dir
for i in rbhstw00189
do
    ./variantdetective.py snp_indel \
    	-1 $dir/${i}/${i}.1.fq.gz \
        -2 $dir/${i}/${i}.2.fq.gz \
    	-r $dir/${i}/*.fa \
    	-t 40 \
    	-o $dir/${i}/
     rm $dir/${i}/*/*bam*
done
