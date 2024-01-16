dir='benchmark-Bmallei-fastq-snp'
mkdir -p $dir

for sra in {SRR1618492,SRR1618671,SRR1618688,SRR1618499,SRR1618349,SRR1616952,SRR2146904,SRR2146899,SRR2146902,SRR2147667,SRR8283094,SRR8072932,SRR8072935,SRR8072938,ERR9616711};
do
echo $sra;
./variantdetective.py snp_indel \
        -1 ${sra}_1.fastq \
	-2 ${sra}_2.fastq \
        -r GCA_000011705.1.fa \
        -t 24 \
        -o $dir/$sra
     rm $dir/$sra/GCA_000011705.1*
     rm $dir/$sra/*fast*
     rm $dir/$sra/*/*bam*
done
