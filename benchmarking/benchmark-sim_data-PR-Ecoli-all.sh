dir='benchmark-PR-Ecoli-all'
mkdir -p $dir
cp parameter_sv_ecoli.txt $dir/parameter_sv_ecoli.txt
cp parameter_snp.txt $dir/parameter_snp.txt
for i in {1..5}
do
  SURVIVOR simSV NC_000913.3.fasta parameter_snp.txt 0.00017 0 $dir/sim_snp_${i}
  SURVIVOR simSV $dir/sim_snp_${i}.fasta parameter_sv_ecoli.txt 0 0 $dir/sim_snp_sv_${i}
  variantdetective all_variants \
    	-g $dir/sim_snp_sv_${i}.fasta \
    	-r NC_000913.3.fasta \
    	-t 24 \
      --readcov 50X \
    	-o $dir/sim_${i}_${j}
  rm $dir/sim_${i}_${j}/NC_000913*
  rm $dir/sim_${i}_${j}/*fast*
  rm $dir/sim_${i}_${j}/*/*bam*
  rm $dir/sim_${i}_${j}/structural_variant/nanovar/*bam*
done
