
mamba install -c anaconda -c conda-forge -c bioconda \
picard \
gatk4 \
freebayes \
clair3 \
vcftools \
svim \
cutesv \
nanovar \
nanosv \
minimap2 \
bedtools \
survivor

picard CreateSequenceDictionary R=GCA_000011705.1.fa O=GCA_000011705.1.dict
picard AddOrReplaceReadGroups I=alignment.mm.sorted.bam O=alignment.mm.rg.sorted.bam RGID=1 RGLB=SAMPLE RGPL=ont RGSM=SAMPLE RGPU=SAMPLE
samtools index alignment.mm.rg.sorted.bam

freebayes -f GCA_000011705.1.fa alignment.mm.rg.sorted.bam > freebayes.vcf

gatk HaplotypeCaller -R GCA_000011705.1.fa -I alignment.mm.rg.sorted.bam -O haplotypecaller.vcf


run_clair3.sh -b alignment.mm.rg.sorted.bam -f GCA_000011705.1.fa -t 20 -p 'ont' -o clair3 -m ~/miniconda3/envs/snp-testing/bin/models/ont --include_all_ctgs --no_phasing_for_fa
cp clair3/merge_output.vcf.gz ./clair3.vcf.gz
gunzip clair3.vcf.gz 
bgzip clair3.vcf
tabix -p vcf clair3.vcf.gz

bgzip freebayes.vcf
tabix -p vcf freebayes.vcf.gz

bgzip haplotypecaller.vcf
tabix -p vcf haplotypecaller.vcf.gz

vcf-isec -p venn freebayes.vcf.gz haplotypecaller.vcf.gz clair3.vcf.gz






minimap2 -t 40 -ax sr ../GCA_000011705.1.fa BmalleiZagreb_S10_L001_R1_001.fastq.gz BmalleiZagreb_S10_L001_R2_001.fastq.gz | samtools view -Sb - -@ 40 -o shortreads_alignment.bam
samtools sort -n -@ 40 shortreads_alignment.bam -o shortreads_alignment.sorted.bam
samtools fixmate -m -@ 40 shortreads_alignment.sorted.bam shortreads_alignment.fm.bam
samtools sort -@ 40 shortreads_alignment.fm.bam -o shortreads_alignment.fm.sorted.bam
samtools markdup -@ 40 -r shortreads_alignment.fm.sorted.bam shortreads_alignment.dup.fm.sorted.bam

picard CreateSequenceDictionary R=GCA_000011705.1.fa O=GCA_000011705.1.dict
picard AddOrReplaceReadGroups I=shortreads_alignment.dup.fm.sorted.bam O=shortreads_alignment.rg.dup.fm.sorted.bam RGID=1 RGLB=SAMPLE RGPL=ont RGSM=SAMPLE RGPU=SAMPLE
samtools index shortreads_alignment.rg.dup.fm.sorted.bam

freebayes -f GCA_000011705.1.fa shortreads_alignment.rg.dup.fm.sorted.bam -p 1 > freebayes.vcf
vcffilter -f "QUAL > 20" freebayes.vcf > freebayes.filt.vcf
bgzip freebayes.filt.vcf
tabix -p vcf freebayes.filt.vcf.gz

gatk HaplotypeCaller -R GCA_000011705.1.fa -I shortreads_alignment.rg.dup.fm.sorted.bam  -O haplotypecaller.vcf -ploidy 1
vcffilter -f "QD > 20" haplotypecaller.vcf > haplotypecaller.filt.vcf
bgzip haplotypecaller.filt.vcf
tabix -p vcf haplotypecaller.filt.vcf.gz

run_clair3.sh -b shortreads_alignment.rg.dup.fm.sorted.bam -f GCA_000011705.1.fa -t 20 -p 'ilmn' -o clair3 -m ~/miniconda3/envs/snp-testing/bin/models/ilmn --include_all_ctgs --no_phasing_for_fa --haploid_precise
cp clair3/merge_output.vcf.gz ./clair3.vcf.gz
gunzip clair3.vcf.gz 
vcffilter -f "QUAL > 20 & FILTER = PASS" clair3.vcf > clair3.filt.vcf
bgzip clair3.filt.vcf
tabix -p vcf clair3.filt.vcf.gz



vcf-isec -p venn_ freebayes.filt.vcf.gz haplotypecaller.filt.vcf.gz clair3.filt.vcf.gz
tabix -p vcf venn_0_1_2.vcf.gz
tabix -p vcf venn_0_1.vcf.gz
tabix -p vcf venn_0_2.vcf.gz
tabix -p vcf venn_1_2.vcf.gz
bcftools concat -a venn_0_1_2.vcf.gz venn_0_1.vcf.gz venn_0_2.vcf.gz  venn_1_2.vcf.gz -o venn_atleast2.vcf
bgzip venn_atleast2.vcf
tabix -p vcf venn_atleast2.vcf.gz

vcf-isec -p comb_ venn_atleast2.vcf.gz ../shortreads_test/venn_atleast2.vcf.gz
