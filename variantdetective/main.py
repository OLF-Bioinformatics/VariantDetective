"""
This file contains the main (entry-point) function into VariantDetective.
It can be run using the variantdetective.py script. 

Copyright (C) 2022 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/philcharron-cfia/VariantDetective
"""

import argparse
import os
import pathlib
import shutil
import sys
import datetime

from variantdetective.tools import get_new_filename
from .validate_inputs import validate_inputs
from .version import __version__

def main(output=sys.stderr):
    check_python_version()
    args = parse_args(sys.argv[1:])

    if args.subparser_name == 'structural_variant':
        check_structural_variant_args(args)
    
    elif args.subparser_name == 'snp_indel':
        check_snp_indel_args(args)
    
    elif args.subparser_name == 'all_variants':
        check_all_variants_args(args)

    create_outdir(args)
    copy_inputs(args)
    print(str(datetime.datetime.now().replace(microsecond=0)) + '\tValidating input files...', file=output)
    validate_inputs(args, output=output)    

def parse_args(args):
    parser = argparse.ArgumentParser(description='VariantDetective: Identify single nucleotide variants (SNV), '
                                        'insertions/deletions (indel) and/or structural variants (SV) from '
                                        'FASTQ reads or FASTA genomic sequences.',
                      formatter_class=NoSubparsersMetavarFormatter,
                      add_help=False)

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help',
                           default=argparse.SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('-v', '--version', action='version',
                           version='VariantDetective v' + __version__,
                           help="Show program version number and exit")

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name',
                                       metavar=None)
    all_variants_subparser(subparsers)
    structural_variant_subparser(subparsers)
    snp_indel_subparser(subparsers)




    # If no arguments were used, print the base-level help.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def structural_variant_subparser(subparsers):
    help = 'Identify structural variants (SV) from long reads (FASTQ) or genome sequence (FASTA). \
                 If input is FASTA, long reads will be simulated to detect SVs.'
    definition = 'Identify structural variants (SV) from long reads (FASTQ) or genome sequence (FASTA). \
                 If input is FASTA, long reads will be simulated to detect SVs.'

    group = subparsers.add_parser('structural_variant', description=definition,
                                  help=help, 
                                  formatter_class=argparse.HelpFormatter,
                                  add_help=False)

    help_args = group.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    help_args.add_argument('-v', '--version', action='version',
                            version='V v' + __version__,
                            help="Show program version number and exit")

    input_args = group.add_argument_group('Input')
    input_args.add_argument('-l', '--long', type=str, metavar='[FASTQ]',
                               help="Path to long reads FASTQ file. Can't be combined with -g")
    input_args.add_argument('-g', '--genome', type=str, metavar='[FASTA]',
                               help="Path to query genomic FASTA file. Can't be combined with -l")
    input_args.add_argument('-r', '--reference', type=str, required=True, metavar='[FASTA]',
                               help='Path to reference genome in FASTA. Required')

    simulate_args = group.add_argument_group('Simulate')
    simulate_args.add_argument("--readcov", type=str, default='50x',
                                help='Either an absolute value (e.g. 250M) or a relative depth (e.g. 50x) (default: %(default)s)')
    simulate_args.add_argument("--readlen", type=str, default='15000,13000',
                                help='Fragment length distribution (mean,stdev) (default: %(default)s)')
    
    nanovar_args = group.add_argument_group('Structural Variant Call')
    nanovar_args.add_argument("--data_type_sv", type=str, default='ont', choices=['ont', 'pacbio'],
                                help="Type of long-read data (ont or pacbio) (default: %(default)s)")
    nanovar_args.add_argument("--mincov_sv", type=int, default=2,
                                help='Minimum number of reads required to call variant (default: %(default)i)')
    nanovar_args.add_argument("--minlen_sv", type=int, default=25,
                                help='Minimum length of SV to be detected (default: %(default)i)')
    nanovar_args.add_argument("--minqual_sv", type=int, default=15,
                                help='Minimum quality of SV to be filtered out from SVIM (default: %(default)i)')
    nanovar_args.add_argument("--sv_consensus", type=int, default=3,
                                help='Specifies the minimum number of tools required to detect an SV to include it in the consensus list (default: %(default)i)')


    other_args = group.add_argument_group('Other')
    other_args.add_argument('-o', "--out", type=str, default='./',
                                help='Output directory. Will be created if it does not exist')
    other_args.add_argument('-t', '--threads', type=int, default=1,
                                help='Number of threads used for job (default: %(default)i)')                            


def all_variants_subparser(subparsers):
    help = 'Identify structural variants (SV) from long reads (FASTQ) and SNPs/indels from short reads (FASTQ). \
        If genome sequence (FASTA) is provided instead, simulate reads and predict SV, SNPs and indels.'
    definition = 'Identify structural variants (SV) from long reads (FASTQ) and SNPs/indels from short reads (FASTQ). \
        If genome sequence (FASTA) is provided instead, simulate reads and predict SV, SNPs and indels.'

    group = subparsers.add_parser('all_variants', description=definition,
                                  help=help, 
                                  formatter_class=argparse.HelpFormatter,
                                  add_help=False)

    help_args = group.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    help_args.add_argument('-v', '--version', action='version',
                            version='VariantDetective v' + __version__,
                            help="Show program version number and exit")

    input_args = group.add_argument_group('Input')
    input_args.add_argument('-l', '--long', type=str, metavar='[FASTQ]',
                               help="Path to long reads FASTQ file. Must be combined with -1 and -2")
    input_args.add_argument('-1', '--short1', type=str, metavar='[FASTQ]',
                               help="Path to pair 1 of short reads FASTQ file. Must be combined with -l and -2")
    input_args.add_argument('-2', '--short2', type=str, metavar='[FASTQ]',
                               help="Path to pair 2 of short reads FASTQ file. Must be combined with -l and -1")
    input_args.add_argument('-g', '--genome', type=str, metavar='[FASTA]',
                               help="Path to query genomic FASTA file. Can't be combined with -l, -1 or -2")
    input_args.add_argument('-r', '--reference', type=str, required=True, metavar='[FASTA]',
                               help='Path to reference genome in FASTA. Required')

    simulate_args = group.add_argument_group('Simulate')
    simulate_args.add_argument("--readcov", type=str, default='50x',
                                help='Either an absolute value (e.g. 250M) or a relative depth (e.g. 50x) (default: %(default)s)')
    simulate_args.add_argument("--readlen", type=str, default='15000,13000',
                                help='Fragment length distribution (mean,stdev) (default: %(default)s)')
    
    nanovar_args = group.add_argument_group('Structural Variant Call')
    nanovar_args.add_argument("--data_type_sv", type=str, default='ont', choices=['ont', 'pacbio'],
                                help="Type of long-read data (ont or pacbio) (default: %(default)s)")
    nanovar_args.add_argument("--mincov_sv", type=int, default=2,
                                help='Minimum number of reads required to call SV (default: %(default)i)')
    nanovar_args.add_argument("--minlen_sv", type=int, default=25,
                                help='Minimum length of SV to be detected (default: %(default)i)')
    nanovar_args.add_argument("--minqual_sv", type=int, default=15,
                                help='Minimum quality of SV to be filtered out with SVIM (default: %(default)i)')
    nanovar_args.add_argument("--sv_consensus", type=int, default=3,
                                help='Specifies the minimum number of tools required to detect an SV to include it in the consensus list (default: %(default)i)')

    snp_args = group.add_argument_group('SNP/Indel Call')
    snp_args.add_argument("--mincov_snp", type=int, default=2,
                                help='Minimum number of reads required to call SNP/Indel (default: %(default)i)')
    snp_args.add_argument("--minqual_snp", type=int, default=20,
                                help='Minimum quality of SNP/Indel to be filtered out (default: %(default)i)')
    snp_args.add_argument("--assembler", type=str, default='bwa', choices=['bwa', 'minimap2'],
                                help='Choose which assembler (bwa or minimap2) to use when using paired-end short reads (default: %(default)s)')
    snp_args.add_argument("--snp_consensus", type=int, default=2,
                                help='Specifies the minimum number of tools required to detect an SNP or Indel to include it in the consensus list (default: %(default)i)')


    other_args = group.add_argument_group('Other')
    other_args.add_argument('-o', "--out", type=str, default='./',
                                help='Output directory. Will be created if it does not exist')
    other_args.add_argument('-t', '--threads', type=int, default=1,
                                help='Number of threads used for job (default: %(default)i)')

def snp_indel_subparser(subparsers):
    help = 'Identify SNPs/indels from short reads (FASTQ). \
        If genome sequence (FASTA) is provided instead, simulate reads and predict SNPs and indels.'
    definition = 'Identify SNPs/indels from short reads (FASTQ). \
        If genome sequence (FASTA) is provided instead, simulate reads and predict  SNPs and indels.'

    group = subparsers.add_parser('snp_indel', description=definition,
                                  help=help, 
                                  formatter_class=argparse.HelpFormatter,
                                  add_help=False)

    help_args = group.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    help_args.add_argument('-v', '--version', action='version',
                            version='VariantDetective v' + __version__,
                            help="Show program version number and exit")

    input_args = group.add_argument_group('Input')
    input_args.add_argument('-1', '--short1', type=str, metavar='[FASTQ]',
                               help="Path to pair 1 of short reads FASTQ file. Must be combined with -2")
    input_args.add_argument('-2', '--short2', type=str, metavar='[FASTQ]',
                               help="Path to pair 2 of short reads FASTQ file. Must be combined with -1")
    input_args.add_argument('-g', '--genome', type=str, metavar='[FASTA]',
                               help="Path to query genomic FASTA file. Can't be combined with -1 or -2")
    input_args.add_argument('-r', '--reference', type=str, required=True, metavar='[FASTA]',
                               help='Path to reference genome in FASTA. Required')

    simulate_args = group.add_argument_group('Simulate')
    simulate_args.add_argument("--readcov", type=str, default='50x',
                                help='Either an absolute value (e.g. 250M) or a relative depth (e.g. 50x) (default: %(default)s)')
    simulate_args.add_argument("--readlen", type=str, default='15000,13000',
                                help='Fragment length distribution (mean,stdev) (default: %(default)s)')
    
    snp_args = group.add_argument_group('SNP/Indel Call')
    snp_args.add_argument("--mincov_snp", type=int, default=2,
                                help='Minimum number of reads required to call SNP/Indel (default: %(default)i)')
    snp_args.add_argument("--minqual_snp", type=int, default=20,
                                help='Minimum quality of SNP/Indel to be filtered out (default: %(default)i)')
    snp_args.add_argument("--assembler", type=str, default='bwa', choices=['bwa', 'minimap2'],
                                help='Choose which assembler (bwa or minimap2) to use when using paired-end short reads (default: %(default)s)')
    snp_args.add_argument("--snp_consensus", type=int, default=2,
                                help='Specifies the minimum number of tools required to detect an SNP or Indel to include it in the consensus list (default: %(default)i)')

    other_args = group.add_argument_group('Other')
    other_args.add_argument('-o', "--out", type=str, default='./',
                                help='Output directory. Will be created if it does not exist')
    other_args.add_argument('-t', '--threads', type=int, default=1,
                                help='Number of threads used for job (default: %(default)i)')


def check_all_variants_args(args):
    if args.long is not None and not pathlib.Path(args.long).is_file():
        sys.exit(f'Error: input file {args.long} does not exist')
    if args.short1 is not None and not pathlib.Path(args.short1).is_file():
        sys.exit(f'Error: input file {args.short1} does not exist')
    if args.short2 is not None and not pathlib.Path(args.short2).is_file():
        sys.exit(f'Error: input file {args.short2} does not exist')        
    if args.genome is not None and not pathlib.Path(args.genome).is_file():
        sys.exit(f'Error: input file {args.genome} does not exist')
    if not pathlib.Path(args.reference).is_file():
        sys.exit(f'Error: reference file {args.reference} does not exist')
      
    if args.long is None and args.short1 is None and args.short2 is None and args.genome is None:
        sys.exit("At least one input must be specified. Must use genomic FASTA (-g) or long read FASTQ (-l), short read pair 1 FASTQ (-1) and short read pair 2 FASTQ (-2).")
    if args.genome is not None and (args.long is not None or args.short1 is not None or args.short2 is not None):
        sys.exit("Cannot use FASTA (-g) with other inputs.")
    if args.genome is None and (args.long is None or args.short1 is None or args.short2 is None):
        sys.exit("Must use long read FASTQ (-l), short read pair 1 FASTQ (-1) and short read pair 2 FASTQ (-2) when calling all variants.")
    if args.mincov_sv < 1:
        sys.exit(f'Error: minimum coverage must be over 1')
    if args.minlen_sv < 1:
        sys.exit('Error: minimum length of SV must be over 1')
    if args.mincov_snp < 1:
        sys.exit(f'Error: minimum coverage must be over 1')
    if args.snp_consensus < 1 or args.snp_consensus > 3:
        sys.exit(f'Error: snp_consensus must be between 1 and 3')
    if args.sv_consensus < 1 or args.sv_consensus > 4:
        sys.exit(f'Error: sv_consensus must be between 1 and 4')
    if args.minqual_sv < 0:
        sys.exit(f'Error: minimum quality of SV must be over 0')
    if args.minqual_snp < 0:
        sys.exit(f'Error: minimum quality of SNP must be over 0')
    try:
        length_parameters = [float(x) for x in args.readlen.split(',')]
        args.mean_frag_length = length_parameters[0]
        args.frag_length_stdev = length_parameters[1]
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --length values')
    if args.mean_frag_length <= 100:
        sys.exit(f'Error: mean read length must be at least 100')
    if args.frag_length_stdev < 0:
        sys.exit('Error: read length stdev cannot be negative')


def check_structural_variant_args(args):
    if args.long is not None and not pathlib.Path(args.long).is_file():
        sys.exit(f'Error: input file {args.long} does not exist')
    if args.genome is not None and not pathlib.Path(args.genome).is_file():
        sys.exit(f'Error: input file {args.genome} does not exist')
    if not pathlib.Path(args.reference).is_file():
        sys.exit(f'Error: reference file {args.reference} does not exist')
    if args.long is None and args.genome is None:
        sys.exit("At least one input must be specified. Must use long-read FASTQ (-l) or genomic FASTA (-g).")
    if args.long is not None and args.genome is not None:
        sys.exit("Only one input can be specified. Can't use FASTQ (-l) and FASTA (-g) together.")
    if args.mincov_sv < 1:
        sys.exit(f'Error: minimum coverage must be over 1')
    if args.minlen_sv < 1:
        sys.exit('Error: minimum length of SV must be over 1')
    if args.sv_consensus < 1 or args.sv_consensus > 4:
        sys.exit(f'Error: sv_consensus must be between 1 and 4')
    if args.minqual_sv < 0:
        sys.exit(f'Error: minimum quality of SV must be over 0')
      
    try:
        length_parameters = [float(x) for x in args.readlen.split(',')]
        args.mean_frag_length = length_parameters[0]
        args.frag_length_stdev = length_parameters[1]
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --length values')
    if args.mean_frag_length <= 100:
        sys.exit(f'Error: mean read length must be at least 100')
    if args.frag_length_stdev < 0:
        sys.exit('Error: read length stdev cannot be negative')

def check_snp_indel_args(args):
    if args.short1 is not None and not pathlib.Path(args.short1).is_file():
        sys.exit(f'Error: input file {args.short1} does not exist')
    if args.short2 is not None and not pathlib.Path(args.short2).is_file():
        sys.exit(f'Error: input file {args.short2} does not exist')        
    if args.genome is not None and not pathlib.Path(args.genome).is_file():
        sys.exit(f'Error: input file {args.genome} does not exist')
    if not pathlib.Path(args.reference).is_file():
        sys.exit(f'Error: reference file {args.reference} does not exist')
      
    if args.short1 is None and args.short2 is None and args.genome is None:
        sys.exit("At least one input must be specified. Must use genomic FASTA (-g) or short read pair 1 FASTQ (-1) and short read pair 2 FASTQ (-2).")
    if args.genome is not None and (args.short1 is not None or args.short2 is not None):
        sys.exit("Cannot use FASTA (-g) with other inputs.")
    if args.genome is None and (args.short1 is None or args.short2 is None):
        sys.exit("Must use short read pair 1 FASTQ (-1) and short read pair 2 FASTQ (-2) when calling SNPs and indels.")
    if args.mincov_snp < 1:
        sys.exit(f'Error: minimum coverage must be over 1')
    if args.snp_consensus < 1 or args.snp_consensus > 3:
        sys.exit(f'Error: snp_consensus must be between 1 and 3')
    if args.minqual_snp < 0:
        sys.exit(f'Error: minimum quality of SNP must be over 0')

    try:
        length_parameters = [float(x) for x in args.readlen.split(',')]
        args.mean_frag_length = length_parameters[0]
        args.frag_length_stdev = length_parameters[1]
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --length values')
    if args.mean_frag_length <= 100:
        sys.exit(f'Error: mean read length must be at least 100')
    if args.frag_length_stdev < 0:
        sys.exit('Error: read length stdev cannot be negative')

def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 6:
        sys.exit('Error: VariantDetective requires Python 3.6 or later')

def copy_file(file1, file2):
    try:
        shutil.copyfile(file1, file2)
    except shutil.SameFileError:
        pass


def copy_inputs(args):
    copy_file(args.reference, get_new_filename(args.reference, args.out))
    if args.genome is not None:
        copy_file(args.genome, get_new_filename(args.genome, args.out))
    if args.subparser_name == 'structural_variant':
        if args.long is not None:
            copy_file(args.long, get_new_filename(args.long, args.out))    
    elif args.subparser_name == 'snp_indel':
        if args.short1 is not None:
            copy_file(args.short1, get_new_filename(args.short1, args.out))
        if args.short2 is not None:
            copy_file(args.short2, get_new_filename(args.short2, args.out))
    else:
        if args.short1 is not None:
            copy_file(args.short1, get_new_filename(args.short1, args.out))
        if args.short2 is not None:
            copy_file(args.short2, get_new_filename(args.short2, args.out))
        if args.long is not None:
            copy_file(args.long, get_new_filename(args.long, args.out))
        
def create_outdir(args):
    if not os.path.isdir(args.out):
        os.makedirs(args.out)
    




class NoSubparsersMetavarFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom
    formatting, in particular for the help texts when dealing with subparsers
    action. It removes subparsers metavar and help line in subcommand argument
    group, and removes extra indentation of those subcommands.
    https://stackoverflow.com/questions/11070268/ 
    """

    def _format_action(self, action):
        result = super()._format_action(action)
        if isinstance(action, argparse._SubParsersAction):
            # fix indentation on first line
            return "%*s%s" % (self._current_indent, "", result.lstrip())
        return result
    def _format_action_invocation(self, action):
        if isinstance(action, argparse._SubParsersAction):
            # remove metavar and help line
            return ""
        return super()._format_action_invocation(action)
    def _iter_indented_subactions(self, action):
        if isinstance(action, argparse._SubParsersAction):
            try:
                get_subactions = action._get_subactions
            except AttributeError:
                pass
            else:
                # remove indentation
                yield from get_subactions()
        else:
            yield from super()._iter_indented_subactions(action)

if __name__ == '__main__':
    main()
