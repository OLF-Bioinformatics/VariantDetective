dir='benchmark-bwa-vs-mm-all'
mkdir -p $dir
cp parameter_sv.txt $dir/parameter_sv.txt
cp parameter_snp.txt $dir/parameter_snp.txt
j=50X
for i in {1..5}
do
  SURVIVOR simSV GCA_000011705.1.fa parameter_snp.txt 0.00017 0 $dir/sim_snp_${i}
  SURVIVOR simSV $dir/sim_snp_${i}.fasta parameter_sv.txt 0 0 $dir/sim_snp_sv_${i}
  # Run using BWA
  ./variantdetective.py all_variants \
    	-g $dir/sim_snp_sv_${i}.fasta \
    	-r GCA_000011705.1.fa \
    	-t 24 \
        --readcov $j \
        --assembler bwa \
    	-o $dir/sim_${i}_${j}_bwa
     rm $dir/sim_${i}_${j}_bwa/GCA_000011705.1*
     rm $dir/sim_${i}_${j}_bwa/*fast*
     rm $dir/sim_${i}_${j}_bwa/*/*bam*
  # Run using minimap2
  ./variantdetective.py all_variants \
    	-g $dir/sim_snp_sv_${i}.fasta \
    	-r GCA_000011705.1.fa \
    	-t 24 \
        --readcov $j \
        --assembler minimap2 \
    	-o $dir/sim_${i}_${j}_mm
     rm $dir/sim_${i}_${j}_mm/GCA_000011705.1*
     rm $dir/sim_${i}_${j}_mm/*fast*
     rm $dir/sim_${i}_${j}_mm/*/*bam*
done