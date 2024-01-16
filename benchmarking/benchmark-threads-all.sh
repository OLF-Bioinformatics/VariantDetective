dir='benchmark-threads-all'
mkdir -p $dir
cp parameter_sv.txt $dir/parameter_sv.txt
for j in {1,2,4,8,12,24,48}
do
  for i in {1..5}
  do 
    SURVIVOR simSV GCA_000011705.1.fa $dir/parameter_sv.txt 0.00017 0 $dir/sim_$i
    start=`date +%s%N`
    ./variantdetective.py all_variants \
    	-g $dir/sim_$i.fasta \
    	-r GCA_000011705.1.fa \
    	-t $j \
    	-o $dir/sim_$i
    end=`date +%s%N`
    rm -r $dir/sim_$i*
    echo `expr $end - $start` >> $dir/results-$dir-$j.txt
  done
done
