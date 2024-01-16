for file in benchmarking/benchmark-Bmallei-fasta-all/*.fasta
do
     name=${file::-6}
     echo ${file}
     ./variantdetective.py all_variants \
        -g ${file} \
        -r benchmarking/benchmark-Bmallei-fasta-all/GCA_000011705.1.fa \
        -t 24 \
        --readcov 50X \
        -o $name
     rm $name/GCA_000011705.1*
     rm $name/*fast*
     rm $name/*/*bam*
done
cat 