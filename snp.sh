freebayes -f GCA_000011705.1.fa alignment.mm.sorted.bam > freebayes.vcf



picard CreateSequenceDictionary R=GCA_000011705.1.fa O=GCA_000011705.1.dict
picard AddOrReplaceReadGroups I=alignment.mm.sorted.bam O=alignment.mm.rg.sorted.bam RGID=1 RGLB=sim RGPL=ont RGSM=sim RGPU=sim
samtools index alignment.mm.rg.sorted.bam
gatk HaplotypeCaller -R GCA_000011705.1.fa -I alignment.mm.rg.sorted.bam -O haplotypecaller.vcf


run_clair3.sh -b alignment.mm.sorted.bam -f GCA_000011705.1.fa -t 20 -p 'ont' -o clair3.vcf -m ~/miniconda3/envs/snp/bin/models/ont --include_all_ctgs --no_phasing_for_fa

mamba install -c anaconda -c conda-forge -c bioconda \
picard \
gatk \
freebayes \
clair3 