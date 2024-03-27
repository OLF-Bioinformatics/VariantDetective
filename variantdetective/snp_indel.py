"""
This file contains code for the snp_indel subcommand.

Copyright (C) 2024 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/OLF-Bioinformatics/VariantDetective
"""

import os
import sys
import io
import pandas as pd
import datetime
import psutil
import subprocess
from .tools import get_new_filename, run_process, read_vcf, generate_tab_csv_snp_summary

def snp_indel(args, snp_input, output=sys.stderr):
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tStarting snp_indel tool', file=output)
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
        if args.assembler == 'minimap2':
            print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning minimap2...', file=output)
            command = 'minimap2 -t ' + str(args.threads) + ' -ax sr '
        elif args.assembler == 'bwa':
            print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning bwa...', file=output)
            command = 'bwa index ' + reference
            run_process(command)
            command = 'bwa mem -t ' + str(args.threads) + ' '
        command += reference + ' ' + snp_input[0] + ' ' + snp_input[1] + \
        ' | samtools view -Sb - -@ ' + str(args.threads) + \
        ' | samtools sort -n - -@ ' +  str(args.threads) + \
        ' | samtools fixmate -m - - -@ ' + str(args.threads) + \
        ' | samtools sort - -@ ' + str(args.threads) + \
        ' | samtools markdup -r - -@ ' + str(args.threads) + ' ' + \
        bam_file_dir + '/alignment.sorted.bam'
        run_process(command)
    elif args.subparser_name == 'snp_indel':
        bam_file_dir = snp_indel_outdir
        rgpl = 'ONT'
        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning minimap2...', file=output)
        command = 'minimap2 -t ' + str(args.threads) + ' -ax map-ont ' + \
        reference + ' ' + snp_input + \
        ' | samtools view -Sb - -@ ' + str(args.threads) + \
        ' | samtools sort - -@ ' + str(args.threads) + \
        ' -o ' + bam_file_dir + '/alignment.sorted.bam'
        run_process(command)
    else:
        bam_file_dir = os.path.join(args.out, 'structural_variant')
        rgpl = 'ONT'

    # Run Picard 
    reference_base = os.path.splitext(reference)[0]
    dict_file = reference_base + ".dict"
    if os.path.exists(dict_file):
        os.remove(dict_file)
    command = 'picard CreateSequenceDictionary R=' + reference
    run_process(command)
    command = 'picard AddOrReplaceReadGroups I=' + bam_file_dir + '/alignment.sorted.bam O=' + \
        snp_indel_outdir + '/alignment.rg.sorted.bam RGID=1 RGLB=SAMPLE RGSM=SAMPLE RGPU=SAMPLE RGPL=' + rgpl
    run_process(command)
    command = 'samtools index ' + snp_indel_outdir + '/alignment.rg.sorted.bam'
    run_process(command)
    try:
        command = 'rm ' + snp_indel_outdir + '/alignment.sorted.bam'
        run_process(command)
    except:
        pass
    try:
        command = 'samtools faidx ' + reference
        run_process(command)
    except:
        pass


    # Run Freebayes
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning Freebayes...', file=output)
    command = 'freebayes-parallel <(fasta_generate_regions.py ' + reference + '.fai 100000) ' + \
            str(args.threads) + ' -f ' + reference + ' ' + \
            snp_indel_outdir + '/alignment.rg.sorted.bam -p 1 > ' + \
            freebayes_outdir + '/freebayes.vcf'
    run_process(command)

    command = 'vcffilter -f "QUAL > ' + str(args.minqual_snp) + '" ' + freebayes_outdir + '/freebayes.vcf > ' + \
        freebayes_outdir + '/freebayes.filt.vcf'
    run_process(command)

    command = 'bgzip -c -@ ' + str(args.threads) + ' ' + \
        freebayes_outdir + '/freebayes.filt.vcf > ' + \
        freebayes_outdir + '/freebayes.filt.vcf.gz'
    run_process(command)

    command = 'tabix -p vcf ' + freebayes_outdir + '/freebayes.filt.vcf.gz'
    run_process(command)

    # Run HaplotypeCaller
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning HaplotypeCaller...', file=output)
    total_memory_gb = psutil.virtual_memory().total / (1024 ** 3)
    ninety_percent_memory_gb = 0.9 * total_memory_gb
    xmx_value = f"-Xmx{int(ninety_percent_memory_gb)}G"
    xss_value = f"-Xss{int(ninety_percent_memory_gb // 10)}M"  # Assuming we use 10% of Xmx for Xss, and convert to MB
    if xss_value == "-Xss0M":
        xss_value = "-Xss1M"
    command = f'gatk HaplotypeCaller --java-options "{xmx_value} {xss_value}" -R ' + reference + \
        ' -I ' + snp_indel_outdir + '/alignment.rg.sorted.bam' + \
        ' -O ' + haplotypecaller_outdir + '/haplotypecaller.vcf' + \
        ' -ploidy 1'
    run_process(command)
    command = 'vcffilter -f "QD > ' + str(args.minqual_snp) + '" ' + haplotypecaller_outdir + '/haplotypecaller.vcf > ' + \
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf'
    run_process(command)

    command = 'bgzip -c -@ ' + str(args.threads) + ' ' + \
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf > ' + \
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz'
    run_process(command)

    command = 'tabix -p vcf ' + haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz'
    run_process(command)

    # Run Clair3
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tRunning Clair3...', file=output)
    #model_path = pkg_resources.resource_filename('variantdetective', 'clair3_models/ilmn')
    
    if args.custom_clair3_model is not None:
        model_path = args.custom_clair3_model

    else:
        command = "dirname $(which run_clair3.sh)"
        # Run the command and capture the output
        try:
            bin_output = subprocess.check_output(command, shell=True, universal_newlines=True).strip()
        except subprocess.CalledProcessError as e:
            print("Error:", e)
        model_path = bin_output + "/models/ilmn"

    command = 'run_clair3.sh -f ' + reference + \
        ' -b ' + snp_indel_outdir + '/alignment.rg.sorted.bam' + \
        ' -o ' + clair3_outdir + \
        ' -p "ilmn" -m ' + model_path + ' --include_all_ctgs ' + \
        ' --no_phasing_for_fa --haploid_precise -t ' + str(args.threads)
    run_process(command)
    
    command = 'mv ' + clair3_outdir + '/merge_output.vcf.gz ' + clair3_outdir + '/clair3.vcf.gz'
    run_process(command)

    command = 'gunzip -f ' + clair3_outdir + '/clair3.vcf.gz'
    run_process(command)

    command = 'vcffilter -f "QUAL > ' + str(args.minqual_snp) + ' & FILTER = PASS" ' + clair3_outdir + '/clair3.vcf > ' + \
        clair3_outdir + '/clair3.filt.vcf'
    run_process(command)

    command = 'bgzip -c -@ ' + str(args.threads) + ' ' + \
        clair3_outdir + '/clair3.filt.vcf > ' + \
        clair3_outdir + '/clair3.filt.vcf.gz'
    run_process(command)

    command = 'tabix -p vcf ' + clair3_outdir + '/clair3.filt.vcf.gz'
    run_process(command)

    # Combine Variants
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tCombining variants...', file=output)
    def has_variants(vcf_file):
        """Check if a VCF file has variants."""
        command = f"zcat {vcf_file} | grep -v '#' | wc -l"
        try:
            result = subprocess.check_output(command, shell=True, universal_newlines=True).strip()
            return int(result) > 0
        except subprocess.CalledProcessError as e:
            print("Error:", e)
            return False

    def create_vcf_if_not_exists(source_file, dest_file, snp_indel_outdir):
        if not os.path.exists(dest_file):
            command = f"zcat {snp_indel_outdir}/{source_file} | grep '#' | bgzip -c > {snp_indel_outdir}/{dest_file}"
            run_process(command)
    
    def run_tabix(vcf_file, snp_indel_outdir):
        command = 'tabix -p vcf ' + snp_indel_outdir + '/' + vcf_file
        run_process(command)

    vcf_files = [
        freebayes_outdir + '/freebayes.filt.vcf.gz',
        haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz',
        clair3_outdir + '/clair3.filt.vcf.gz'
    ]
    
    # Filter out files with no variants
    valid_vcf_files = [file for file in vcf_files if has_variants(file)]

    if len(valid_vcf_files) > 0:
        command = 'vcf-isec -p ' + snp_indel_outdir + '/snp_ ' + ' '.join(valid_vcf_files)
        run_process(command)
        if len(valid_vcf_files) == 3:
            source_vcf = 'snp_0_1_2.vcf.gz'
            dest_vcfs = ['snp_0_1.vcf.gz', 'snp_0_2.vcf.gz', 'snp_1_2.vcf.gz',
                        'snp_0.vcf.gz', 'snp_1.vcf.gz', 'snp_2.vcf.gz']
            run_tabix(source_vcf, snp_indel_outdir)
            for dest_vcf in dest_vcfs:
                create_vcf_if_not_exists(source_vcf, dest_vcf, snp_indel_outdir)
                run_tabix(dest_vcf, snp_indel_outdir)
        elif len(valid_vcf_files) == 2:
            source_vcf = 'snp_0_1.vcf.gz'
            dest_vcfs = ['snp_0_1_2.vcf.gz', 'snp_0_2.vcf.gz', 'snp_1_2.vcf.gz',
                        'snp_0.vcf.gz', 'snp_1.vcf.gz', 'snp_2.vcf.gz']
            run_tabix(source_vcf, snp_indel_outdir)
            for dest_vcf in dest_vcfs:
                create_vcf_if_not_exists(source_vcf, dest_vcf, snp_indel_outdir)
                run_tabix(dest_vcf, snp_indel_outdir) 
        elif len(valid_vcf_files) == 1:
            source_vcf = 'snp_0.vcf.gz'
            dest_vcfs = ['snp_0_1_2.vcf.gz', 'snp_0_1.vcf.gz', 'snp_0_2.vcf.gz',
                        'snp_1_2.vcf.gz', 'snp_1.vcf.gz', 'snp_2.vcf.gz']
            run_tabix(source_vcf, snp_indel_outdir)
            for dest_vcf in dest_vcfs:
                create_vcf_if_not_exists(source_vcf, dest_vcf, snp_indel_outdir)
                run_tabix(dest_vcf, snp_indel_outdir)
        if args.snp_consensus == 3:
            command = 'gunzip -c ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz > ' + \
                snp_indel_outdir + '/snp_final.vcf' 
            run_process(command)
        elif args.snp_consensus == 2:
            command = 'bcftools concat -a ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + \
                snp_indel_outdir + '/snp_0_1.vcf.gz ' + \
                snp_indel_outdir + '/snp_0_2.vcf.gz ' + \
                snp_indel_outdir + '/snp_1_2.vcf.gz ' + \
                '-o ' + snp_indel_outdir + '/snp_final.vcf'
            run_process(command)
        elif args.snp_consensus == 1:
            command = 'bcftools concat -a ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + \
                snp_indel_outdir + '/snp_0_1.vcf.gz ' + \
                snp_indel_outdir + '/snp_0_2.vcf.gz ' + \
                snp_indel_outdir + '/snp_1_2.vcf.gz ' + \
                snp_indel_outdir + '/snp_0.vcf.gz ' + \
                snp_indel_outdir + '/snp_1.vcf.gz ' + \
                snp_indel_outdir + '/snp_2.vcf.gz ' + \
                '-o ' + snp_indel_outdir + '/snp_final.vcf'
            run_process(command)
        if len(valid_vcf_files) == 3:    
            command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
                    'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
        elif len(valid_vcf_files) == 2:
            missing_vcf = list(set(vcf_files) - set(valid_vcf_files))
            if freebayes_outdir + '/freebayes.filt.vcf.gz' in missing_vcf:
                command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
                    'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
            elif haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz' in missing_vcf:
                command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
                    'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
            elif clair3_outdir + '/clair3.filt.vcf.gz'in missing_vcf:
                command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
                    'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
        elif len(valid_vcf_files) == 1:
            if freebayes_outdir + '/freebayes.filt.vcf.gz' in valid_vcf_files:
                command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
                    'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
            elif haplotypecaller_outdir + '/haplotypecaller.filt.vcf.gz' in valid_vcf_files:
                command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
                    'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
            elif clair3_outdir + '/clair3.filt.vcf.gz'in valid_vcf_files:
                command = 'mv ' + snp_indel_outdir + '/snp_0.vcf.gz ' + snp_indel_outdir + '/clair3.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1.vcf.gz ' + snp_indel_outdir + '/freebayes.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.unique.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1.vcf.gz ' + snp_indel_outdir + '/freebayes.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_2.vcf.gz ' + snp_indel_outdir + '/haplotypecaller.clair3.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.vcf.gz ; ' + \
                    'mv ' + snp_indel_outdir + '/snp_0_1_2.vcf.gz ' + snp_indel_outdir + '/freebayes.haplotypecaller.clair3.vcf.gz ; ' + \
                    'rm ' + snp_indel_outdir + '/*.tbi ' + snp_indel_outdir + '/snp__README'
        run_process(command)
    else:
        command = 'gunzip -c ' + vcf_files[0] + ' > ' + snp_indel_outdir + '/snp_final.vcf' 
        run_process(command)
    
    generate_tab_csv_snp_summary(read_vcf(snp_indel_outdir + '/snp_final.vcf'), snp_indel_outdir)