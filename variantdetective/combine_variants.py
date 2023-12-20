"""
This file contains code for the combine_variants subcommand.

Copyright (C) 2024 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/OLF-Bioinformatics/VariantDetective
"""

import os
import sys
import datetime
from .tools import run_process, read_vcf, generate_tab_csv_snp_summary, generate_tab_csv_sv_summary


def combine_variants(args, vcf_lists, output=sys.stderr):
    if len(vcf_lists[0]) != 0:
        snp_indel_outdir = os.path.join(args.out, 'snp_indel')
        if not os.path.isdir(snp_indel_outdir):
            os.makedirs(snp_indel_outdir)

        out_vcf_list = []
        for vcf_file in vcf_lists[0]:
            basename = os.path.basename(vcf_file)
            out_vcf_file = snp_indel_outdir + '/' + basename + '.gz'
            out_vcf_list.append(out_vcf_file)
            command = 'bgzip -c -@ ' + str(args.threads) + ' ' + \
                vcf_file + ' > ' + out_vcf_file
            run_process(command)

            command = 'tabix -p vcf ' + out_vcf_file
            run_process(command)
        out_vcf_string = ' '.join(out_vcf_list)
        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tCombining SNP VCF files...', file=output)
        command = 'vcf-isec -p ' + snp_indel_outdir + '/snp_ ' + out_vcf_string
        run_process(command)
        
        snp_output_files = []
        for path, currentDirectory, files in os.walk(snp_indel_outdir):
            for file in files:
                if file.startswith("snp_") and file.endswith(".vcf.gz"):
                    snp_output_files.append(file) 


        for snp_file in snp_output_files:
            command = 'tabix -f -p vcf ' + snp_indel_outdir + '/' + snp_file
            run_process(command)

        filtered_list = [filename for filename in snp_output_files if filename.count('_') >= args.snp_consensus]
        updated_list = [snp_indel_outdir + "/" + filename for filename in filtered_list]

        filtered_string = ' '.join(updated_list)
        if len(updated_list) > 0:
            command = 'bcftools concat -a ' + filtered_string + \
            ' -o ' + snp_indel_outdir + '/snp_final.vcf'
            run_process(command)

        command = 'rm ' + snp_indel_outdir + '/*.tbi '
        run_process(command)

        generate_tab_csv_snp_summary(read_vcf(snp_indel_outdir + '/snp_final.vcf'), snp_indel_outdir)

    if len(vcf_lists[1]) != 0:
        structural_variant_outdir = os.path.join(args.out, 'structural_variant')
        if not os.path.isdir(structural_variant_outdir):
            os.makedirs(structural_variant_outdir)

        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tCombining SV VCF files...', file=output)
        vcf_string = ' '.join(vcf_lists[1])
        command = 'ls ' + vcf_string + ' > ' + structural_variant_outdir + '/vcf_list'
        run_process(command)

        command = 'SURVIVOR merge ' + structural_variant_outdir + \
            '/vcf_list 1000 ' + str(args.sv_consensus) + ' 1 1 0 ' + str(args.minlen_sv) + ' ' \
            + structural_variant_outdir + '/combined_sv.vcf'
        run_process(command)

        command = 'rm ' + structural_variant_outdir + '/vcf_list'
        run_process(command)

        generate_tab_csv_sv_summary(read_vcf(structural_variant_outdir + '/combined_sv.vcf'), structural_variant_outdir)