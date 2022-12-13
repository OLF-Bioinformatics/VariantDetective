"""
This file contains code for the snp_indel subcommand.

Copyright (C) 2022 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/philcharron-cfia/VariantDetective
"""

import os
import sys
import io
import pandas as pd
from .tools import get_new_filename, run_process

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def generate_tab_csv_summary(vcf, output_dir):
    CHROM = vcf.iloc[:,0]
    POS = vcf.iloc[:,1]
    FORMAT = vcf.iloc[:,9].str.split(':', expand=True)
    INFO = vcf.iloc[:,7].str.split(';', expand=True)
    REF = vcf.iloc[:,3]
    ALT = vcf.iloc[:,4]
    SUPPORT = pd.Series([None] * len(vcf), name='SUPPORT')
    TYPE = INFO.iloc[:,40].str.split('=',expand=True).iloc[:,1]
    TYPE.name = 'TYPE'
    for i, v in TYPE.items():
        if (v == None):
            if (len(REF[i]) < len(ALT[i])):
                TYPE[i] = 'ins'
            elif (len(REF[i]) > len(ALT[i])):
                TYPE[i] = 'del'
            elif (len(REF[i]) == 1):
                TYPE[i] = 'snp'
            else:
                TYPE[i] = 'complex'
            RAW_SUPPORT = FORMAT.iloc[i,1].split(",")
            SUPPORT[i] = "REF=" + RAW_SUPPORT[0] + ";ALT=" + RAW_SUPPORT[1]
        else:
            RAW_SUPPORT = FORMAT.iloc[i,2].split(",")
            SUPPORT[i] = "REF=" + RAW_SUPPORT[0] + ";ALT=" + RAW_SUPPORT[1]

    TAB_DATA = pd.concat([CHROM, POS, TYPE, REF, ALT, SUPPORT], axis=1)

    VARIANT_TYPES = pd.Series(['SNP', 'DEL', 'INS', 'MNP', 'COMPLEX','TOTAL'], name = 'TYPE')
    VARIANT_DATA = pd.Series([(TYPE=='snp').sum(), (TYPE=='del').sum(), (TYPE=='ins').sum(), (TYPE=='mnp').sum(), (TYPE=='complex').sum(), len(TYPE)], name = 'COUNT')
    SUMMARY_DATA = pd.concat([VARIANT_TYPES, VARIANT_DATA], axis=1)

    TAB_DATA.to_csv(output_dir + '/snp_final.csv', index=False)
    TAB_DATA.to_csv(output_dir + '/snp_final.tab', sep='\t', index=False)
    SUMMARY_DATA.to_csv(output_dir + '/snp_final_summary.txt', sep='\t', index=False)


def snp_indel(args, snp_input, output=sys.stderr):
    print('Running snp_indel tool', file=output)
    reference = get_new_filename(args.reference, args.out)   
    snp_indel_outdir = os.path.join(args.out, 'snp_indel')
    haplotypecaller_outdir = os.path.join(snp_indel_outdir, 'haplotypecaller')
    freebayes_outdir = os.path.join(snp_indel_outdir, 'freebayes')
    clair3_outdir = os.path.join(snp_indel_outdir, 'clair3')
    
    if not os.path.isdir(snp_indel_outdir):
        os.makedirs(snp_indel_outdir)
    if not os.path.isdir(haplotypecaller_outdir):
        os.makedirs(haplotypecaller_outdir)
    if not os.path.isdir(freebayes_outdir):
        os.makedirs(freebayes_outdir)
    if not os.path.isdir(clair3_outdir):
        os.makedirs(clair3_outdir)

    # Map reads if using short reads
    if isinstance(snp_input, list):
        bam_file_dir = snp_indel_outdir
        rgpl = 'ILLUMINA'
        print('Running minimap2...', end=' ', file=output)
        command = 'minimap2 -t ' + str(args.threads) + ' -ax sr ' + \
        reference + ' ' + snp_input[0] + ' ' + snp_input[1] + \
        ' | samtools view -Sb - -@ ' + str(args.threads) + \
        ' | samtools sort -n - -@ ' +  str(args.threads) + \
        ' | samtools fixmate -m - - -@ ' + str(args.threads) + \
        ' | samtools sort - -@ ' + str(args.threads) + \
        ' | samtools markdup -r - -@ ' + str(args.threads) + ' ' + \
        bam_file_dir + '/alignment.mm.sorted.bam'
        run_process(command, "Error: Issue with mapping short reads")
        print('Complete', file=output)
    else:
        bam_file_dir = os.path.join(args.out, 'structural_variant')
        rgpl = 'ONT'

    # Run Picard 
    command = 'picard CreateSequenceDictionary R=' + reference
    run_process(command, "Error: Picard CreateSequenceDictionary failed")
    command = 'picard AddOrReplaceReadGroups I=' + bam_file_dir + '/alignment.mm.sorted.bam O=' + \
        snp_indel_outdir + '/alignment.mm.rg.sorted.bam RGID=1 RGLB=SAMPLE RGSM=SAMPLE RGPU=SAMPLE RGPL=' + rgpl
    run_process(command, "Error: Picard AddOrReplaceReadGroups failed")
    command = 'samtools index ' + snp_indel_outdir + '/alignment.mm.rg.sorted.bam'
    run_process(command, "Error: Samtools index failed")
    try:
        command = 'rm ' + snp_indel_outdir + '/alignment.mm.sorted.bam'
        run_process(command, "Error: Samtools index failed")
    except:
        pass
    
    # Run Freebayes
    print('Running Freebayes...', end=' ', file=output)
    command = 'freebayes -f ' + reference + ' ' + \
            snp_indel_outdir + '/alignment.mm.rg.sorted.bam -p 1 > ' + \
            freebayes_outdir + '/freebayes.vcf'
    run_process(command, "Error: Issue with Freebayes")

    command = 'vcffilter -f "QUAL > 20" ' + freebayes_outdir + '/freebayes.vcf > ' + \
        freebayes_outdir + '/freebayes.filt.vcf'
    run_process(command, "Error: Issue with filtering Freebayes VCF")

    command = 'bgzip -c -@ ' + str(args.threads) + ' ' + \
        freebayes_outdir + '/freebayes.filt.vcf > ' + \
        freebayes_outdir + '/freebayes.filt.vcf.gz'
    run_process(command, "Error: Issue with BGZIP Freebayes VCF")

    command = 'tabix -p vcf ' + freebayes_outdir + '/freebayes.filt.vcf.gz'
    run_process(command, "Error: Issue with Tabix Freebayes VCF")
    print('Complete', file=output)

    # Run HaplotypeCaller
    print('Running HaplotypeCaller...', end=' ', file=output)
    command = 'gatk HaplotypeCaller -R ' + reference + \
        ' -I ' + snp_indel_outdir + '/alignment.mm.rg.sorted.bam' + \
        ' -O ' + haplotypecaller_outdir + '/haplotypecaller.vcf' + \
        ' -ploidy 1'
    run_process(command, "Error: Issue with HaplotypeCaller")

    command = 'vcffilter -f "QD > 20" ' + haplotypecaller_outdir + '/haplotypecaller.vcf > ' + \
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf'
    run_process(command, "Error: Issue with filtering HaplotypeCaller VCF")

    command = 'bgzip -c -@ ' + str(args.threads) + ' ' + \
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf > ' + \
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz'
    run_process(command, "Error: Issue with BGZIP HaplotypeCaller VCF")

    command = 'tabix -p vcf ' + haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz'
    run_process(command, "Error: Issue with Tabix HaplotypeCaller VCF")
    print('Complete', file=output)

    # Run Clair3
    print('Running Clair3...', end=' ', file=output)
    command = 'run_clair3.sh -f ' + reference + \
        ' -b ' + snp_indel_outdir + '/alignment.mm.rg.sorted.bam' + \
        ' -o ' + clair3_outdir + \
        ' -p "ilmn" -m variantdetective/clair3_models/ilmn --include_all_ctgs ' + \
        ' --no_phasing_for_fa --haploid_precise -t ' + str(args.threads)
    run_process(command, "Error: Issue with Clair3")
    
    command = 'mv ' + clair3_outdir + '/merge_output.vcf.gz ' + clair3_outdir + '/clair3.vcf.gz'
    run_process(command, "Error: Issue with Clair3")

    command = 'gunzip ' + clair3_outdir + '/clair3.vcf.gz'
    run_process(command, "Error: Issue with Clair3")

    command = 'vcffilter -f "QUAL > 20 & FILTER = PASS" ' + clair3_outdir + '/clair3.vcf > ' + \
        clair3_outdir + '/clair3.filt.vcf'
    run_process(command, "Error: Issue with filtering Clair3 VCF")

    command = 'bgzip -c -@ ' + str(args.threads) + ' ' + \
        clair3_outdir + '/clair3.filt.vcf > ' + \
        clair3_outdir + '/clair3.filt.vcf.gz'
    run_process(command, "Error: Issue with BGZIP Clair3 VCF")

    command = 'tabix -p vcf ' + clair3_outdir + '/clair3.filt.vcf.gz'
    run_process(command, "Error: Issue with Tabix Clair3 VCF")
    print('Complete', file=output)

    # Combine Variants
    print('Combining variants...', end=' ', file=output)
    command = 'vcf-isec -p ' + snp_indel_outdir + '/snp_ ' + \
        freebayes_outdir + '/freebayes.filt.vcf.gz ' + \
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz ' + \
        clair3_outdir + '/clair3.filt.vcf.gz'
    run_process(command, "Error: Issue with vcf-isec")

    command = 'tabix -p vcf ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz && ' + \
        'tabix -p vcf ' + snp_indel_outdir + '/snp_0_1.vcf.gz && ' + \
        'tabix -p vcf ' + snp_indel_outdir + '/snp_0_2.vcf.gz && ' + \
        'tabix -p vcf ' + snp_indel_outdir + '/snp_1_2.vcf.gz'
    run_process(command, "Error: Issue with Tabix")

    command = 'bcftools concat -a ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + \
        snp_indel_outdir + '/snp_0_1.vcf.gz ' + \
        snp_indel_outdir + '/snp_0_2.vcf.gz ' + \
        snp_indel_outdir + '/snp_1_2.vcf.gz ' + \
        '-o ' + snp_indel_outdir + '/snp_final.vcf'
    run_process(command, "Error: Issue with bcftools")
    
    command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + \
        snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
        'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + \
        snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
        'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + \
        snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
        'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + \
        snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
        'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + \
        snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
        'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + \
        snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
        'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + \
        snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
        'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
    run_process(command, "Error: Issue with renaming VCF files")
    
    generate_tab_csv_summary(read_vcf(snp_indel_outdir + '/snp_final.vcf'), snp_indel_outdir)
    print('Complete', file=output)
    