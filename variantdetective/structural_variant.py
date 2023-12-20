"""
This file contains code for the structural_variant subcommand.

Copyright (C) 2024 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/OLF-Bioinformatics/VariantDetective
"""

import os
import sys
import io
import pandas as pd
import datetime
from .tools import get_new_filename, run_process, read_vcf, generate_tab_csv_sv_summary

def structural_variant(args, input_reads, output=sys.stderr):
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tStarting structural_variant tool', file=output)
    reference = get_new_filename(args.reference, args.out)   
    structural_variant_outdir = os.path.join(args.out, 'structural_variant')
    nanovar_outdir = os.path.join(structural_variant_outdir, 'nanovar')
    nanosv_outdir = os.path.join(structural_variant_outdir, 'nanosv')
    svim_outdir = os.path.join(structural_variant_outdir, 'svim')
    cutesv_outdir = os.path.join(structural_variant_outdir, 'cutesv')
    
    if not os.path.isdir(structural_variant_outdir):
        os.makedirs(structural_variant_outdir)
    if not os.path.isdir(nanovar_outdir):
        os.makedirs(nanovar_outdir)
    if not os.path.isdir(nanosv_outdir):
        os.makedirs(nanosv_outdir)
    if not os.path.isdir(svim_outdir):
        os.makedirs(svim_outdir)
    if not os.path.isdir(cutesv_outdir):
        os.makedirs(cutesv_outdir)

    # Run NanoVar 
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning NanoVar...', file=output)
    if args.long is not None and args.data_type_sv == "pacbio":
        command = 'nanovar -t ' + str(args.threads) +  ' ' + \
            input_reads  + ' ' + \
            reference + ' ' + \
            nanovar_outdir + \
            ' -x pacbio-clr' + \
            ' -c ' + str(args.mincov_sv) + \
            ' -l ' + str(args.minlen_sv)
    else:
        command = 'nanovar -t ' + str(args.threads) +  ' ' + \
            input_reads  + ' ' + \
            reference + ' ' + \
            nanovar_outdir + \
            ' -c ' + str(args.mincov_sv) + \
            ' -l ' + str(args.minlen_sv)
    run_process(command)

    # Sort bam file
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tSorting BAM file...', file=output)
    command = 'samtools sort -@ ' + str(args.threads) + ' ' + \
            nanovar_outdir + '/*-mm.bam > ' + \
            structural_variant_outdir + '/alignment.sorted.bam'
    run_process(command)
    
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tIndexing sorted BAM file...', file=output)
    command = 'samtools index ' + structural_variant_outdir + '/alignment.sorted.bam'
    run_process(command)

    # Run NanoSV
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning NanoSV...', file=output)
    command = 'samtools faidx ' + reference
    run_process(command)
    
    command = 'cut -f 1,2 ' + reference + '.fai > ' + nanosv_outdir + '/chrom.sizes'
    run_process(command)

    command = 'bedtools random -l 1 -g ' + nanosv_outdir + '/chrom.sizes > ' + nanosv_outdir + '/reference.bed'
    run_process(command)

    command = 'NanoSV -t ' + str(args.threads) + \
            ' -o ' + nanosv_outdir + '/variants.vcf ' + \
            ' -s samtools' + \
            ' -b ' + nanosv_outdir + '/reference.bed ' + \
            structural_variant_outdir + '/alignment.sorted.bam' 
    run_process(command)
    
    # Run SVIM
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning SVIM...', file=output)
    command = 'svim alignment ' + svim_outdir + ' ' + \
            structural_variant_outdir + '/alignment.sorted.bam ' + \
            reference + \
            ' --min_sv_size ' + str(args.minlen_sv)
    run_process(command)
    
    command = "bcftools view -i 'QUAL >= " + str(args.minqual_sv) + "' " + \
            svim_outdir + '/variants.vcf > ' + \
            svim_outdir + '/variants.filt.vcf'
    run_process(command)

    # Run CuteSV
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning CuteSV...', file=output)
    command = 'cuteSV ' + structural_variant_outdir + '/alignment.sorted.bam ' + \
            reference + ' ' + \
            cutesv_outdir + '/variants.vcf ' + \
            cutesv_outdir + \
            ' -t ' + str(args.threads) + \
            ' -s ' + str(args.mincov_sv) + \
            ' -l ' + str(args.minlen_sv) + \
            ' -L -1'  
    run_process(command)
    
    # Run SURVIVOR
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning SURVIVOR...', file=output)
    command = 'ls ' + nanovar_outdir + '/*pass.vcf ' + \
            cutesv_outdir + '/variants.vcf ' + \
            nanosv_outdir + '/variants.vcf ' + \
            svim_outdir + '/variants.filt.vcf > ' + \
            structural_variant_outdir + '/vcf_list'
    run_process(command)

    command = 'SURVIVOR merge ' + structural_variant_outdir + \
            '/vcf_list 1000 ' + str(args.sv_consensus) + ' 1 1 0 ' + str(args.minlen_sv) + ' ' \
            + structural_variant_outdir + '/combined_sv.vcf'
    run_process(command)

    generate_tab_csv_sv_summary(read_vcf(structural_variant_outdir + '/combined_sv.vcf'), structural_variant_outdir)
    command = 'rm ' + structural_variant_outdir + '/vcf_list'
    run_process(command)
    
    return structural_variant_outdir + '/alignment.sorted.bam'
