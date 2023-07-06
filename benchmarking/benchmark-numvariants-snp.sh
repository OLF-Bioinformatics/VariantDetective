dir='benchmark-numvariants-snp'
mkdir -p $dir
cp parameter_snp.txt $dir/parameter_snp.txt
for j in {0.000017,0.000034,0.00017,0.00034,0.0017,0.0034,0.017}
do
  for i in {1..5}
  do 
    SURVIVOR simSV GCA_000011705.1.fa $dir/parameter_snp.txt $j 0 $dir/sim_$i
    start=`date +%s%N`
    ./variantdetective.py snp_indel \
    	-g $dir/sim_$i.fasta \
    	-r GCA_000011705.1.fa \
    	-t 24 \
    	-o $dir/sim_$i
    end=`date +%s%N`
    rm -r $dir/sim_$i*
    echo `expr $end - $start` >> $dir/results-$dir-$j.txt
  done
done
