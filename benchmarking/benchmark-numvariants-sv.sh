dir='benchmark-numvariants-sv'
mkdir -p $dir
for j in {21,46,66,128,210}
do
  for i in {1..5}
  do 
    SURVIVOR simSV GCA_000011705.1.fa $dir/parameter_sv_$j.txt 0 0 $dir/sim_$i
    start=`date +%s%N`
    ./variantdetective.py structural_variant \
    	-g $dir/sim_$i.fasta \
    	-r GCA_000011705.1.fa \
    	-t 24 \
        --readcov 50X \
    	-o $dir/sim_$i
    end=`date +%s%N`
    rm -r $dir/sim_$i*
    echo `expr $end - $start` >> $dir/results-$dir-$j.txt
  done
done
