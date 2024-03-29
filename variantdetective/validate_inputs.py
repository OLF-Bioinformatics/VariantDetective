"""
This file contains code needed to validate the inputs used
to ensure they are in the right format for VariantDetective.

Copyright (C) 2024 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/OLF-Bioinformatics/VariantDetective
"""

import datetime
import os
import sys

from .combine_variants import combine_variants
from .simulate import simulate
from .snp_indel import snp_indel
from .structural_variant import structural_variant
from .tools import get_fasta_info, get_fastq_info, get_input_type, get_new_filename, get_open_function

def validate_inputs(args, output=sys.stderr):
    input_file_type = []
    snp_vcf_list = []
    sv_vcf_list = []
    #actual_file_type = []
   
    if 'genome' in args and args.genome is not None:
        genome_file = get_new_filename(args.genome, args.out)
        input_file_type.append("Genomic FASTA")
        #actual_file_type.append(check_input(genome_file, output))
        #actual_file_type.append("Genomic FASTA")
    if 'long' in args and args.long is not None:
        long_file = get_new_filename(args.long, args.out)
        input_file_type.append("Long-read FASTQ")
        #actual_file_type.append(check_input(long_file, output))
        #actual_file_type.append("Long-read FASTQ")
    if 'short1' in args and args.short1 is not None:
        short1_file = get_new_filename(args.short1, args.out)
        input_file_type.append("Short-read FASTQ")
        #actual_file_type.append(check_input(short1_file, output))
        #actual_file_type.append("Short-read FASTQ")
    if 'short2' in args and args.short2 is not None:
        short2_file = get_new_filename(args.short2, args.out)
        input_file_type.append("Short-read FASTQ")
        #actual_file_type.append(check_input(short2_file, output))
        #actual_file_type.append("Short-read FASTQ")
    if 'snp_vcf' in args and args.snp_vcf is not None:
        for vcf_file in args.snp_vcf:
            snp_vcf_list.append(get_new_filename(vcf_file, args.out))
        input_file_type.append("SNP VCF")
    if 'sv_vcf' in args and args.sv_vcf is not None:
        for vcf_file in args.sv_vcf:
            sv_vcf_list.append(get_new_filename(vcf_file, args.out))
        input_file_type.append("SV VCF")


    #if input_file_type == actual_file_type:
    if 'Genomic FASTA' in input_file_type:
        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tStarting genome pipeline')
        sim_file = simulate(args, genome_file, output=sys.stderr)
        if args.subparser_name == "structural_variant":
            long_bam_file = structural_variant(args, sim_file, output=sys.stderr)
        if args.subparser_name == "snp_indel":
            snp_indel(args, sim_file, output=sys.stderr)
        if args.subparser_name == "all_variants":
            long_bam_file = structural_variant(args, sim_file, output=sys.stderr)
            snp_indel(args, long_bam_file, output=sys.stderr)
    elif 'Long-read FASTQ' in input_file_type and 'Short-read FASTQ' in input_file_type:
        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tStarting short and long read pipeline')
        short_inputs = [short1_file, short2_file]
        long_bam_file = structural_variant(args, long_file, output=sys.stderr)
        snp_indel(args, short_inputs, output=sys.stderr)
    elif 'Long-read FASTQ' in input_file_type:
        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tStarting long read pipeline')
        long_bam_file = structural_variant(args, long_file, output=sys.stderr)
    elif 'Short-read FASTQ' in input_file_type:
        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tStarting short read pipeline')
        short_inputs = [short1_file, short2_file]
        snp_indel(args, short_inputs, output=sys.stderr)
    elif 'SNP VCF' or 'SV VCF' in input_file_type:
        print(str(datetime.datetime.now().replace(microsecond=0)) + '\tStarting combine variants tool')
        vcf_lists = [snp_vcf_list, sv_vcf_list]
        combine_variants(args, vcf_lists, output=sys.stderr)
    
    #else:
    #    for i in range(len(input_file_type)):
    #        if input_file_type[i] == "Long-read FASTQ":
    #            if actual_file_type[i] == "Genomic FASTA":
    #                message = 'Input file was supposed to be long-read FASTQ but genomic FASTA was detected.'
    #            elif actual_file_type[i] == "Short-read FASTQ":
    #                message = 'Input file was supposed to be long-read FASTQ but short-read FASTQ was detected.'
    #        elif input_file_type[i] == "Genomic FASTA":
    #            if actual_file_type[i] == "Long-read FASTQ":
    #                message = 'Input file was supposed to be genomic FASTA but long-read FASTQ was detected.'
    #            elif actual_file_type[i] == "Short-read FASTQ":
    #                message = 'Input file was supposed to be genomic FASTA but short-read FASTQ was detected.'
    #        elif input_file_type[i] == "Short-read FASTQ":
    #            if actual_file_type[i] == "Long-read FASTQ":
    #                message = 'Input file was supposed to be short-read FASTQ but long-read FASTQ was detected.'
    #            elif actual_file_type[i] == "Genomic FASTA":
    #                message = 'Input file was supposed to be short-read FASTQ but genomic FASTA was detected.'
    #        message = message + ' Please verify inputs or use appropriate tool and parameters.'
    #        raise Exception(message)

def check_input(file, output=sys.stderr):
    file_extension = os.path.splitext(file)
    open_func = get_open_function(file_extension[1])
    file_type = get_input_type(open_func, file)

    if file_type == 'FASTQ':
        count, min_value, max_value, median, average = get_fastq_info(open_func, file)
        if average > 301:
            actual_file_type = "Long-read FASTQ"
            print("Input file type:\tLong-read FASTQ", file=output)
        elif average > 0:
            actual_file_type = "Short-read FASTQ"
            print("Input file type:\tShort-read FASTQ", file=output)
        else:
            raise Exception('Average length of reads is 0')
        print("Number of reads:\t{}".format(count), file=output)

    elif file_type == 'FASTA':
        count = get_fasta_info(open_func, file)
        actual_file_type = "Genomic FASTA"
    
    return actual_file_type






