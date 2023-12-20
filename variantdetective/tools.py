"""
This file contains miscellaneous tools that are used in various
components of VariantDetective.

Copyright (C) 2024 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/OLF-Bioinformatics/VariantDetective
"""

import gzip
import io
import os
import pandas as pd
import statistics
from subprocess import Popen, PIPE, STDOUT

def get_fasta_info(open_func, file):
    count = 0 
    with open_func(file, 'rt') as seq_file:
        for line in seq_file:
            if line.startswith('>'):
                count += 1
    return count


def get_fastq_info(open_func, file):
    num_lines = sum(1 for i in open_func(file, 'rb'))
    if num_lines % 4 != 0:
        raise ValueError('File might be corrupted, unexpected number of lines was found')
    else:
        with open_func(file, 'rt') as seq_file:
            length_list = []
            for i in range(int(num_lines/4)):
                next(seq_file)
                length_list.append(len(next(seq_file).strip()))
                next(seq_file)
                next(seq_file)
        count = int(num_lines/4)
        min_value = min(length_list)
        max_value = max(length_list)
        median = int(statistics.median(length_list))
        average = int(statistics.mean(length_list))
        return count, min_value, max_value, median, average


def get_input_type(open_func, file):
    with open_func(file, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''
    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        raise ValueError('File is not FASTA or FASTQ')

def get_new_filename(file, out):
    basename = os.path.basename(file)
    new_filename = os.path.join(out, basename)
    return new_filename


def get_open_function(file_extension):
    if file_extension == ".gz":
        open_func = gzip.open
    else:
        open_func = open
    return open_func

def run_process(command):
    process = Popen([command],
                    universal_newlines=True, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
    output, error = process.communicate()
    
    if process.returncode != 0:
        raise Exception(error)

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def generate_tab_csv_snp_summary(vcf, output_dir):
    CHROM = vcf.iloc[:,0]
    POS = vcf.iloc[:,1]
    FORMAT_ID = vcf.iloc[:,8].str.split(':', expand=True)
    FORMAT = vcf.iloc[:,9].str.split(':', expand=True)
    INFO = vcf.iloc[:,7].str.split(';', expand=True)
    REF = vcf.iloc[:,3]
    ALT = vcf.iloc[:,4]
    SUPPORT = pd.Series([None] * len(vcf), name='SUPPORT')
    TYPE = INFO.iloc[:,40].str.split('=',expand=True).iloc[:,1]
    TYPE.name = 'TYPE'
    for i, v in TYPE.items():
        try:
            AD_index =list(FORMAT_ID.iloc[i]).index('AD')
            RAW_SUPPORT = FORMAT.iloc[i,AD_index].split(",")
            SUPPORT[i] = "REF=" + RAW_SUPPORT[0] + ";ALT=" + RAW_SUPPORT[1]
        except ValueError:
            AF_index =list(FORMAT_ID.iloc[i]).index('AF')
            DP_index =list(FORMAT_ID.iloc[i]).index('DP')
            REF_COUNT = round(int(FORMAT.iloc[i,DP_index]) - (float(FORMAT.iloc[i,AF_index]) * int(FORMAT.iloc[i,DP_index])))
            ALT_COUNT = round(float(FORMAT.iloc[i,AF_index]) * int(FORMAT.iloc[i,DP_index]))
            SUPPORT[i] = "REF=" + str(REF_COUNT) + ";ALT=" + str(ALT_COUNT)
        if (v == None):
            if (len(REF[i]) < len(ALT[i])):
                TYPE[i] = 'ins'
            elif (len(REF[i]) > len(ALT[i])):
                TYPE[i] = 'del'
            elif (len(REF[i]) == 1):
                TYPE[i] = 'snp'
            else:
                TYPE[i] = 'complex'

    TAB_DATA = pd.concat([CHROM, POS, TYPE, REF, ALT, SUPPORT], axis=1)

    VARIANT_TYPES = pd.Series(['SNP', 'DEL', 'INS', 'MNP', 'COMPLEX','TOTAL'], name = 'TYPE')
    VARIANT_DATA = pd.Series([(TYPE=='snp').sum(), (TYPE=='del').sum(), (TYPE=='ins').sum(), (TYPE=='mnp').sum(), (TYPE=='complex').sum(), len(TYPE)], name = 'COUNT')
    SUMMARY_DATA = pd.concat([VARIANT_TYPES, VARIANT_DATA], axis=1)

    TAB_DATA.to_csv(output_dir + '/snp_final.csv', index=False)
    TAB_DATA.to_csv(output_dir + '/snp_final.tab', sep='\t', index=False)
    SUMMARY_DATA.to_csv(output_dir + '/snp_final_summary.txt', sep='\t', index=False)

def generate_tab_csv_sv_summary(vcf, output_dir):
    CHROM = vcf.iloc[:,0]
    CHROM.name = 'REF_CHROM'
    START = vcf.iloc[:,1]
    START.name = 'REF_START'
    END = vcf.iloc[:,7].str.split(';', expand=True).iloc[:,6].str[4:]
    END.name = 'REF_STOP'
    SIZE = vcf.iloc[:,7].str.split(';', expand=True).iloc[:,2].str[6:]
    SIZE.name = 'SIZE'
    TYPE = vcf.iloc[:,7].str.split(';', expand=True).iloc[:,3].str[7:]
    TYPE.name = 'TYPE'
    INFO = vcf.iloc[:,4].copy()
    INFO.name = 'INFO'
    for i, v in INFO.items():
        if (TYPE[i] == "TRA"):
            if ("]" not in v and "[" not in v):
                INFO[i] = ''
            elif ("]" in v):
                INFO[i] = ']'+(v.split(']'))[1].split(']')[0]+']'
            elif ("[" in v):
                INFO[i] = '['+(v.split('['))[1].split('[')[0]+'['
            else:
                INFO[i] = ''
        else:
            INFO[i] = ''

    TAB_DATA = pd.concat([CHROM, START, END, SIZE, TYPE, INFO], axis=1)

    VARIANT_TYPES = pd.Series(['TRANSLOCATION', 'INVERSION', 'DELETION', 'INSERTION', 'DUPLICATION','TOTAL'], name = 'TYPE')
    VARIANT_DATA = pd.Series([(TYPE=='TRA').sum(), (TYPE=='INV').sum(), (TYPE=='DEL').sum(), (TYPE=='INS').sum(), (TYPE=='DUP').sum(), len(TYPE)], name = 'COUNT')
    SUMMARY_DATA = pd.concat([VARIANT_TYPES, VARIANT_DATA], axis=1)
    
    TAB_DATA.to_csv(output_dir + '/combined_sv.csv', index=False)
    TAB_DATA.to_csv(output_dir + '/combined_sv.tab', sep='\t', index=False)
    SUMMARY_DATA.to_csv(output_dir + '/combined_sv_summary.txt', sep='\t', index=False)