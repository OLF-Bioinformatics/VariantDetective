dir='benchmark-real_data_HG002-sv-v2'
mkdir -p $dir
./variantdetective.py structural_variant \
    -l $dir/NA12878-12.5mil.fq.gz \
    -r $dir/hg19_ucsc_main.fa \
    -t 24 \
    -o $dir/NA12878_results

