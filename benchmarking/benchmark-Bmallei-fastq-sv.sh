dir='benchmark-Bmallei-fastq-sv'
mkdir -p $dir

for sra in {SRR1618494,SRR1618669,SRR1618686,SRR1618500,SRR1618350,SRR1617359,SRR2146906,SRR2146901,SRR2146903,SRR2147669,SRR8283092,SRR8072934,SRR8072936,SRR8072939,ERR9616715};
do
echo $sra;
./variantdetective.py structural_variant \
        -l ${sra}.fastq \
        -r GCA_000011705.1.fa \
        -t 24 \
        -o $dir/$sra
     rm $dir/$sra/GCA_000011705.1*
     rm $dir/$sra/*fast*
     rm $dir/$sra/*/*bam*
done
