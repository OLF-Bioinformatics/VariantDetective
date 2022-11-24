<p align='center'><img src='images/logo.png' alt="VariantDetective" width="75%"></p>

This program is designed to identify short variants and structural variants. Variants can be identified from genomic sequences (FASTA) or from combinations of short and/or long reads (FASTQ). If genomic sequences are provided as input, long reads will be simulated to detect variants.

This tool makes use of other open-source variant callers and creates consensus sets (support from at least 2 callers) in order to validate a variant. Summary files for short variants and structural variants are generated outlining the different types of variants found in the sample.

## Installation

### Conda Installation
`conda install -c bioconda -c charronp variantdetective`

## Quick 

## List of Commands

- <b>[all_variants](#all_variants):</b> Identify structural variants (SV) from long reads (FASTQ) and SNPs/indels from short reads (FASTQ), or both types of variants from genome sequence (FASTA). If genome sequence (FASTA) is provided, reads will be simulated to predict SV, SNPs and indels.
- <b>[structural_variant](structural_variant):</b> Identify structural variants (SV) from long reads (FASTQ) or genome sequence (FASTA).  If genome sequence (FASTA) is provided, reads will be simulated to predict SVs.
- <b>[snp_indel](snp_indel):</b> Identify SNPs/indels from short reads (FASTQ) or genome sequence (FASTA). If genome sequence (FASTA) is provided instead, reads will be simulated to predict SNPs and indels.

## Commands and Options

