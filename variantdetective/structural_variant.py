import os
import sys
from subprocess import Popen, PIPE, STDOUT

from .tools import get_new_filename

def run_process(command, error_message):
    process = Popen([command],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT, shell=True, executable='/bin/bash')
    exitcode = process.wait()
    if exitcode != 0:
        raise Exception(error_message)
    

def structural_variant(args, input_reads, output=sys.stderr):
    print('Running structural_variant tool...', file=output)
    reference = get_new_filename(args.reference, args.out)   
    structural_variant_outdir = os.path.join(args.out, 'structural_variant')
    nanovar_outdir = os.path.join(structural_variant_outdir, 'nanovar')
    nanosv_outdir = os.path.join(structural_variant_outdir, 'nanosv')
    svim_outdir = os.path.join(structural_variant_outdir, 'svim')
    cutesv_outdir= os.path.join(structural_variant_outdir, 'cutesv')
    
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
    print('Running NanoVar...', file=output)
    command = 'nanovar -t ' + str(args.threads) +  ' ' + \
            input_reads  + ' ' + \
            reference + ' ' + \
            nanovar_outdir + \
            ' -c ' + str(args.mincov_sv) + \
            ' -l ' + str(args.minlen_sv)
    print(command, file=output)
    run_process(command, "Error: NanoVar failed")
    
    # Sort bam file
    print('Sorting BAM file...', file=output)
    command = 'samtools sort -@ ' + str(args.threads) + ' ' + \
            nanovar_outdir + '/*-mm.bam > ' + \
            structural_variant_outdir + '/alignment.mm.sorted.bam'
    run_process(command, "Error: Issue with sorting BAM file")

    print('Indexing sorted BAM file...', file=output)
    command = 'samtools index ' + structural_variant_outdir + '/alignment.mm.sorted.bam'
    run_process(command, "Error: Issue with indexing BAM file")

    # Run NanoSV
    print('Running NanoSV...', file=output)
    command = 'samtools faidx ' + reference
    run_process(command, "Error: Issue with generating index for reference genome")
    
    command = 'cut -f 1,2 ' + reference + '.fai > ' + nanosv_outdir + '/chrom.sizes'
    run_process(command, "Error: Issue with getting reference genome chromosome sizes")

    command = 'bedtools random -l 1 -g ' + nanosv_outdir + '/chrom.sizes > ' + nanosv_outdir + '/reference.bed'
    run_process(command, "Error: Issue with generating reference genome bed file")

    command = 'NanoSV -t ' + str(args.threads) + \
            ' -o ' + nanosv_outdir + '/variants.vcf ' + \
            ' -s samtools' + \
            ' -b ' + nanosv_outdir + '/reference.bed ' + \
            structural_variant_outdir + '/alignment.mm.sorted.bam' 
    print(command, file=output)
    run_process(command, "Error: NanoSV failed")
    
    # Run SVIM
    print('Running SVIM...', file=output)
    command = 'svim alignment ' + svim_outdir + ' ' + \
            structural_variant_outdir + '/alignment.mm.sorted.bam ' + \
            reference + \
            ' --min_sv_size ' + str(args.minlen_sv)
    print(command, file=output)
    run_process(command, "Error: SVIM failed")
    
    command = "bcftools view -i 'QUAL >= 15' " + \
            svim_outdir + '/variants.vcf > ' + \
            svim_outdir + '/variants.filt.vcf'
    run_process(command, "Error: filtering SVIM output failed")

    # Run CuteSV
    print('Running CuteSV...', file=output)
    command = 'cuteSV ' + structural_variant_outdir + '/alignment.mm.sorted.bam ' + \
            reference + ' ' + \
            cutesv_outdir + '/variants.vcf ' + \
            cutesv_outdir + \
            ' -t ' + str(args.threads) + \
            ' -s ' + str(args.mincov_sv) + \
            ' -l ' + str(args.minlen_sv) + \
            ' -L -1'  
    print(command, file=output)
    run_process(command, "Error: CuteSV failed")

    # Run SURVIVOR
    print('Running SURVIVOR...', file=output)
    command = 'ls ' + nanovar_outdir + '/*pass.vcf ' + \
            cutesv_outdir + '/variants.vcf ' + \
            nanosv_outdir + '/variants.vcf ' + \
            svim_outdir + '/variants_15.vcf > ' + \
            structural_variant_outdir + '/vcf_list'
    run_process(command, "Error: couldn't find VCF files")

    command = 'SURVIVOR merge ' + structural_variant_outdir + \
            '/vcf_list 1000 2 1 1 0 ' + str(args.minlen_sv) + ' ' \
            + structural_variant_outdir + '/combined_sv_variants.vcf'
    print(command, file=output)
    run_process(command, "Error: SURVIVOR failed")