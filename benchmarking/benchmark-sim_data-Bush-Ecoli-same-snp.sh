dir='benchmark-sim_data-Bush-Ecoli-same-snp'
for i in {09-00049,1223,14EC017,14EC020,2012C-4227,2013C-4465,2016C-3878,210221272,2149,28RC1,317,350,51008369SK1,5CRE51,746,789,9000,94-3024,95JB1,ABWA45,AR_0006,AR_0017,AR_0055,AR_0069,AR_0151,AR_0369,AR436,ATCC25922,BH100Lsubstr.MG2017}
do
    ./variantdetective.py snp_indel \
    	-g $dir/${i}/${i}_mutated_simulation.fasta \
        -r $dir/${i}/${i}_simulated.fasta \
      	-t 40 \
    	-o $dir/${i}/ \
		--readcov 100x \
		--readlen 100000,2000 \
    rm $dir/${i}/*/*bam*
done


