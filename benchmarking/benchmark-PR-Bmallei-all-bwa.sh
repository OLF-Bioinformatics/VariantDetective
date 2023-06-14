dir='benchmark-PR-Bmallei-bwa-all'
mkdir -p $dir
cp parameter_sv.txt $dir/parameter_sv.txt
cp parameter_snp.txt $dir/parameter_snp.txt
j=50X
for i in {1,2,3,4,5}
do
  #SURVIVOR simSV GCA_000011705.1.fa parameter_snp.txt 0.00017 0 $dir/sim_snp_${i}
  #SURVIVOR simSV $dir/sim_snp_${i}.fasta parameter_sv.txt 0 0 $dir/sim_snp_sv_${i}
    ./variantdetective.py snp_indel \
    	-g $dir/sim_snp_sv_${i}.fasta \
    	-r GCA_000011705.1.fa \
    	-t 24 \
        --readcov $j \
    	-o $dir/sim_${i}_${j}
     rm $dir/sim_${i}_${j}/GCA_000011705.1*
     rm $dir/sim_${i}_${j}/*fast*
     rm $dir/sim_${i}_${j}/*/*bam*


done
