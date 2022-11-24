<p align='center'><img src='images/logo.png' alt="VariantDetective" width="75%"></p>

This program is designed to identify short variants and structural variants. Variants can be identified from genomic sequences (FASTA) or from combinations of short and/or long reads (FASTQ). If genomic sequences are provided as input, long reads will be simulated to detect variants.

This tool makes use of other open-source variant callers and creates consensus sets (support from at least 2 callers) in order to validate a variant. Summary files for short variants and structural variants are generated outlining the different types of variants found in the sample.

## Table of Contents
  - [Requirements](#requirements)
  - [Installation](#installation)
     - [Conda Installation](#conda-installation)
  - [Quick Usage](#quick-usage)
  - [List of Commands](#list-of-commands)
  - [Detailed usage](#detailed-usage)
     - [all_variants](#allvariants)
     - [structural_variant](#structuralvariant)
     - [snp_indel](#snpindel)
  - [Reporting Issues](#reporting-issues)
  - [Citing VariantDetective](#citing-variantdetective)

## Requirements
MUST ADD REQUIREMENTS

## Installation
ADD INSTALLATION FROM SOURCE
### Conda Installation
`conda install -c bioconda -c charronp variantdetective`

## Quick Usage
MUST ADD QUICK USAGE

## List of Commands

- <b>[all_variants](#all_variants):</b> Identify structural variants (SV) from long reads (FASTQ) and SNPs/indels from short reads (FASTQ), or both types of variants from genome sequence (FASTA). If genome sequence (FASTA) is provided, reads will be simulated to predict SV, SNPs and indels.
- <b>[structural_variant](#structural_variant):</b> Identify structural variants (SV) from long reads (FASTQ) or genome sequence (FASTA).  If genome sequence (FASTA) is provided, reads will be simulated to predict SVs.
- <b>[snp_indel](#snp_indel):</b> Identify SNPs/indels from short reads (FASTQ) or genome sequence (FASTA). If genome sequence (FASTA) is provided instead, reads will be simulated to predict SNPs and indels.

## Detailed Usage

### all_variants

```
usage: variantdetective-runner.py all_variants -r FASTA [-1 FASTQ] [-2 FASTQ] [-l FASTQ] [-g FASTA]                            
                                               [--readcov READCOV] [--readlen READLEN]
                                               [--mincov_sv MINCOV_SV] [--minlen_sv MINLEN_SV]
                                               [--mincov_snp MINCOV_SNP] [-o OUT] [-t THREADS]
                                               [-h] [-v] 

Identify structural variants (SV) from long reads (FASTQ) and SNPs/indels from
short reads (FASTQ). If genome sequence (FASTA) is provided instead, simulate
reads and predict SV, SNPs and indels.

Input:
  -r FASTA, --reference FASTA
                        Path to reference genome in FASTA. Required
  -1 FASTQ, --short1 FASTQ
                        Path to pair 1 of short reads FASTQ file. Must be
                        combined with -2 and -l
  -2 FASTQ, --short2 FASTQ
                        Path to pair 2 of short reads FASTQ file. Must be
                        combined with -1 and -l
  -l FASTQ, --long FASTQ
                        Path to long reads FASTQ file. Must be combined with
                        -1 and -2
  -g FASTA, --genome FASTA
                        Path to query genomic FASTA file. Can't be combined
                        with -l, -1 or -2
  
Simulate:
  --readcov READCOV     Either an absolute value (e.g. 250M) or a relative
                        depth (e.g. 50x) (default: 50x)
  --readlen READLEN     Fragment length distribution (mean,stdev) (default:
                        15000,13000)

Structural Variant Call:
  --mincov_sv MINCOV_SV
                        Minimum number of reads required to call SV (default:
                        2)
  --minlen_sv MINLEN_SV
                        Minimum length of SV to be detected (default: 25)

SNP/Indel Call:
  --mincov_snp MINCOV_SNP
                        Minimum number of reads required to call SNP/Indel
                        (default: 2)

Other:
  -o OUT, --out OUT     Output directory. Will be created if it does not exist
  -t THREADS, --threads THREADS
                        Number of threads used for job (default: 1)

Help:
  -h, --help            Show this help message and exit
  -v, --version         Show program version number and exit
```

### structural_variant

```
usage: variantdetective-runner.py structural_variant -r FASTA [-l FASTQ] [-g FASTA]                            
                                               [--readcov READCOV] [--readlen READLEN]
                                               [--mincov_sv MINCOV_SV] [--minlen_sv MINLEN_SV]
                                               [-o OUT] [-t THREADS] [-h] [-v] 

Identify structural variants (SV) from long reads (FASTQ) or genome sequence
(FASTA). If input is FASTA, long reads will be simulated to detect SVs.

Input:
  -r FASTA, --reference FASTA
                        Path to reference genome in FASTA. Required
  -l FASTQ, --long FASTQ
                        Path to long reads FASTQ file. Can't be combined with
                        -g
  -g FASTA, --genome FASTA
                        Path to query genomic FASTA file. Can't be combined
                        with -l
  
Simulate:
  --readcov READCOV     Either an absolute value (e.g. 250M) or a relative
                        depth (e.g. 50x) (default: 50x)
  --readlen READLEN     Fragment length distribution (mean,stdev) (default:
                        15000,13000)

Structural Variant Call:
  --mincov_sv MINCOV_SVLicense
                        Minimum length of SV to be detected (default: 25)

Other:
  -o OUT, --out OUT     Output directory. Will be created if it does not exist
  -t THREADS, --threads THREADS
                        Number of threads used for job (default: 1)

Help:
  -h, --help            Show this help message and exit
  -v, --version         Show program version number and exit
```

### snp_indel

```
usage: variantdetective-runner.py structural_variant -r FASTA [-1 FASTQ] [-2 FASTQ] [-g FASTA]                                                    
                                               [--readcov READCOV] [--readlen READLEN]
                                               [--mincov_snp MINCOV_SNP] [-o OUT]
                                               [-t THREADS] [-h] [-v]
                                               
Identify SNPs/indels from short reads (FASTQ). If genome sequence (FASTA) is
provided instead, simulate reads and predict SNPs and indels.

Input:
  -r FASTA, --reference FASTA
                        Path to reference genome in FASTA. Required
  -1 [FASTQ], --short1 [FASTQ]
                        Path to pair 1 of short reads FASTQ file. Must be
                        combined with -2
  -2 [FASTQ], --short2 [FASTQ]
                        Path to pair 2 of short reads FASTQ file. Must be
                        combined with -1
  -g FASTA, --genome FASTA
                        Path to query genomic FASTA file. Can't be combined
                        with -1 or -2
  
Simulate:
  --readcov READCOV     Either an absolute value (e.g. 250M) or a relative
                        depth (e.g. 50x) (default: 50x)
  --readlen READLEN     Fragment length distribution (mean,stdev) (default:
                        15000,13000)

SNP/Indel Call:
  --mincov_snp MINCOV_SNP
                        Minimum number of reads required to call SNP/Indel
                        (default: 2)

Other:
  -o OUT, --out OUT     Output directory. Will be created if it does not exist
  -t THREADS, --threads THREADS
                        Number of threads used for job (default: 1)

Help:
  -h, --help            Show this help message and exit
  -v, --version         Show program version number and exit
```

## Reporting Issues

If you have any issues installing or running VariantDetective, or would like a new feature added to the tool, please open an issue here on GitHub. 

## Citing VariantDetective

Citation will be included later.