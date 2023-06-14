dir='benchmark-real_data_HG002-sv'
mkdir -p $dir
./variantdetective.py structural_variant \
    -l $dir/NA12878.pacbio.sanger.fq.gz \
    -r $dir/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    -t 24 \
    -o $dir/NA12878_results
