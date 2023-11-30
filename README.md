<p align='center'><img src='res/logo.png' alt="VariantDetective" width="75%"></p>

This program is designed to identify short variants and structural variants. Variants can be identified from genomic sequences (FASTA) or from combinations of short and/or long reads (FASTQ). If genomic sequences are provided as input, long reads will be simulated to detect variants.

This tool makes use of other open-source variant callers and creates consensus sets in order to validate a variant. Summary files for short variants and structural variants are generated outlining the different types of variants found in the sample.

## Author

Phil Charron \<<phil.charron@inspection.gc.ca>\>

## Table of Contents
  - [Installation](#installation)
     - [Installation from Source](#installation-from-source)
     - [Conda Installation](#conda-installation)
  - [Quick Usage](#quick-usage)
  - [List of Commands](#list-of-commands)
  - [Variant Callers](#variant-callers)
     - [Short Variant Callers](#short-variant-callers)
     - [Structural Variant Callers](#structural-variant-callers)
  - [Long Read Simulator](#long-read-simulator)
  - [Parameters](#parameters)
  - [Outputs](#outputs)
     - [Output files - snp_indel directory](#output-files---snp_indel-directory)
     - [Output files - structural_variant directory](#output-files---structural_variant-directory)  
  - [Reporting Issues](#reporting-issues)
  - [Citing VariantDetective](#citing-variantdetective)


## Installation

All software and tools used by VariantDetective can be found in the [spec-file.txt](spec-file.txt). VariantDetective can be installed via pip after creating the conda environment to support it or via conda.

### Installation from Source
VariantDetective can be installed from source using the following method.
```
# Download VariantDetective repository
git clone https://github.com/OLF-Bioinformatics/VariantDetective.git
cd VariantDetective
# Create conda variant for tools
conda create -n variantdetective -y && conda activate variantdetective
# Install specific versions of tools
conda install -n variantdetective --file spec-file.txt
# Install VariantDetective
pip install -e .
```

### Conda Installation
```
conda create -n vd -y
conda activate vd
conda install -c bioconda -c conda-forge -c charronp variantdetective
```

## Quick Usage

**Finding snps/indels and structural variants from an assembled genome (FASTA)**

```
variantdetective all_variants -r REFERENCE.fasta -g SAMPLE.fasta
```

**Finding snps/indels and structural variants from raw reads (FASTQ)**

```
variantdetective all_variants -r REFERENCE.fasta -1 SHORT_READ_1.fastq -2 SHORT_READ_2.fastq -l LONG_READ.fastq
```

**Finding snps/indels from an assembled genome (FASTA)**

```
variantdetective snp_indel -r REFERENCE.fasta -g SAMPLE.fasta
```

**Finding snps/indels from raw reads (FASTQ)**

```
variantdetective snp_indel -r REFERENCE.fasta -1 SHORT_READ_1.fastq -2 SHORT_READ_2.fastq 
```

**Finding structural variants from an assembled genome (FASTA)**

```
variantdetective structural_variant -r REFERENCE.fasta -g SAMPLE.fasta
```

**Finding structural variants from raw reads (FASTQ)**

```
variantdetective structural_variant -r REFERENCE.fasta -l LONG_READ.fastq
```

## List of Commands
| Command | Description |
| --- | --- |
| `variantdetective all_variants` | Identify structural variants (SV) from long reads (FASTQ) and SNPs/indels from short reads (FASTQ), or both types of variants from genome sequence (FASTA). If genome sequence (FASTA) is provided, reads will be simulated to predict SV, SNPs and indels. |
| `variantdetective structural_variant` | Identify structural variants (SV) from long reads (FASTQ) or genome sequence (FASTA).  If genome sequence (FASTA) is provided, reads will be simulated to predict SVs. |
| `variantdetective snp_indel` |  Identify SNPs/indels from short reads (FASTQ) or genome sequence (FASTA). If genome sequence (FASTA) is provided instead, reads will be simulated to predict SNPs and indels. |

## Variant Callers
VariantDetective makes use of published open-source variant callers and creates consensus sets in order to validate a variant.

### Short Variant Callers
- [Clair3](https://github.com/HKU-BAL/Clair3)
- [Freebayes](https://github.com/freebayes/freebayes)
- [GATK4 HaplotypeCaller](https://github.com/broadinstitute/gatk)

Intersections of VCF files are created using the [VCFtools](https://github.com/vcftools) `vcf-isec` tool. The final VCF output consensus file containing variants found in at least 2 variant callers (default) is created using the [BCFtools](https://github.com/samtools/bcftools) `concat` tool.

### Structural Variant Callers
- [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [NanoSV](https://github.com/mroosmalen/nanosv)
- [NanoVar](https://github.com/cytham/nanovar)
- [SVIM](https://github.com/eldariont/svim)

The consensus VCF file is created using the [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) `merge` tool. Parameters for merging structural variants are a maximum allowed distance of 1 kbp between breakpoints and calls supported by at least 3 variant callers (default) where they agree on both type and strand.  

## Long Read Simulator
When a genomic FASTA file is provided as query input, long reads are simulated in order to detect variants. The long read simulation tool is adapted from [Badread](https://github.com/rrwick/Badread), a tool that creates simulated reads. It has been modified to create perfectly matching reads to the reference file and to allow multi-threading to speed up the process.

## Parameters
All input files can be uncompressed (.fasta/.fastq) or gzipped (.fastq.gz/.fastq.gz)

| Options | Available Command | Description | Default | 
| --- | :---: | --- | :---: |
| `-r FASTA` | `all_variants`<br>`structural_variant`<br>`snp_indel` | Path to reference genome in FASTA. Required | - |
| `-g FASTA` |  `all_variants`<br>`structural_variant`<br>`snp_indel` | Path to query genomic FASTA file. Can't be combined with `-1`, `-2` or `-l`| - |
| `-1 FASTQ`<br>`--short1 FASTQ` | `all_variants`<br>`snp_indel` | Path to pair 1 of short reads FASTQ file. Must always be combined with `-2`. If running `all_variants`, must be combined with `-l`| - |
| `-2 FASTQ`<br>`--short2 FASTQ` | `all_variants`<br>`snp_indel` | Path to pair 2 of short reads FASTQ file. Must always be combined with `-1`. If running `all_variants`, must be combined with `-l`| - |
| `-l FASTQ`<br>`--long FASTQ` | `all_variants`<br>`structural_variant` | Path to long reads FASTQ file. If running `all_variants`, must be combined with `-1` and `-2`| - |
| `--readcov READCOV` |  `all_variants`<br>`structural_variant`<br>`snp_indel` | If using `-g` as input, define the absolute amount of simulated reads (e.g. 250M) or relative simulated read depth (e.g. 50x) | `50x` | 
| `--readlen MEAN,STDEV` |  `all_variants`<br>`structural_variant`<br>`snp_indel` | If using `-g` as input, define the mean length and standard deviation of simulated reads | `15000,13000` |
| `--mincov_snp MINCOV_SNP` |  `all_variants`<br>`snp_indel` | Minimum number of reads required to call SNP/Indel | `2` |
| `--minqual_snp MINQUAL_SNP` |  `all_variants`<br>`snp_indel` | Minimum quality of SNP/Indel to be filtered out | `20` |
| `--assembler {bwa,minimap2}` |  `all_variants`<br>`snp_indel` | Choose which assembler (bwa or minimap2) to use when using paired-end short reads | `bwa` |
| `--snp_consensus SNP_CONSENSUS` |  `all_variants`<br>`snp_indel` | Specifies the minimum number of tools required to detect an SNP or Indel to include it in the consensus list | `2` |
| `--mincov_sv MINCOV_SV` | `all_variants`<br>`structural_variant` | Minimum number of reads required to call SV  | `2` |
| `--minlen_sv MINLEN_SV` | `all_variants`<br>`structural_variant` | Minimum length of SV to be detected | `25` |
| `--minqual_sv MINQUAL_SV` |  `all_variants`<br>`structural_variant` | Minimum quality of SV to be filtered out from SVIM | `15` |
| `--sv_consensus SV_CONSENSUS` |  `all_variants`<br>`structural_variant` | Specifies the minimum number of tools required to detect an SV to include it in the consensus list | `3` |
| `-o OUT`<br>`--out OUT` | `all_variants`<br>`structural_variant`<br>`snp_indel` |  Output directory. Will be created if it does not exist | `./` |
| `-t THREADS` <br> `--threads THREADS` | `all_variants` <br> `structural_variant` <br> `snp_indel` | Number of threads used for job | `1` |
| `-h` <br> `--help` | `all_variants`<br>`structural_variant`<br>`snp_indel` | Show help message and exit | - |
| `-v` <br> `--version` |  `all_variants`<br>`structural_variant`<br>`snp_indel`| Show program version number and exit | - |

## Outputs

All input files will be copied to the output folder. Within the output folder, directories containing the `structural_variant` and `snp_indel` results will be created.

### Output files - `snp_indel` directory

| Output |  Description |
|---:|---|
| `snp_final.vcf` | Variants that were found in at least 2 variant callers in VCF format |
| `snp_final.csv` | Variants that were found in at least 2 variant callers in CSV format |
| `snp_final.tab` | Variants that were found in at least 2 variant callers in TSV format |
| `snp_final_summary.txt` | Summary of different short variant types found in snp_final files |
| `freebayes.haplotypecaller.clair3.vcf.gz` | Variants in common between all variants callers |
| `freebayes.clair3.vcf.gz` | Variants in common between Freebayes and Clair3 |
| `freebayes.haplotypecaller.vcf.gz` | Variants in common between Freebayes and HaplotypeCaller |
| `haplotypecaller.unique.vcf.gz` | Variants in common between HaplotypeCaller and Clair3 |
| `clair3.unique.vcf.gz` | Variants only found by Clair3 |
| `freebayes.unique.vcf.gz` | Variants only found by Freebayes |
| `haplotypecaller.unique.vcf.gz` | Variants only found by HaplotypeCaller | 
| `alignment.mm.rg.sorted.bam` |  Alignment in BAM format |
| `alignment.mm.rg.sorted.bam.bai` | Index file of alignments |
| `clair3/` | Directory containing files related to Clair3 variant calling |
| `freebayes/` | Directory containing files related to Freebayes variant calling |
| `haplotypecaller/` | Directory containing files related to HaplotypeCaller variant calling |

### Output files - `structural_variant` directory

| Output |  Description |
|---:|---|
| `combined_sv.vcf` | Variants that were found in at least 2 variant callers in VCF format |
| `combined_sv.csv` | Variants that were found in at least 2 variant callers in CSV format |
| `combined_sv.tab` | Variants that were found in at least 2 variant callers in TSV format |
| `combined_sv_summary.txt` | Summary of different structural variant types found in combined_sv files |
| `alignment.mm.sorted.bam` |  Alignment in BAM format |
| `alignment.mm.sorted.bam.bai` | Index file of alignments |
| `cutesv/` | Directory containing files related to cuteSV variant calling |
| `nanosv/` | Directory containing files related to NanoSV variant calling |
| `nanovar/` | Directory containing files related to NanoVar variant calling |
| `svim/` | Directory containing files related to SVIM variant calling |


## Reporting Issues

If you have any issues installing or running VariantDetective, or would like a new feature added to the tool, please open an issue here on GitHub. 

## Citing VariantDetective

Our manuscript describing this tool is currently under review.