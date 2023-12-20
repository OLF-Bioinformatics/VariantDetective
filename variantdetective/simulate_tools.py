"""
This file contains functions that are used in the long read simulation
component of VariantDetective.

Copyright (C) 2024 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/OLF-Bioinformatics/VariantDetective
Portions Copyright (C) 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread
"""

import collections
import gzip
import os
import random
import re
import sys

def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'

REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a',
                 'g': 'c', 'c': 'g', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 
                 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', 
                 '.': '.', '-': '-', '?': '?'}

def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type

def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open

def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA.
    """
    if not os.path.isfile(filename):
        sys.exit('Error: could not find {}'.format(filename))
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(filename, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''
    if first_char == '>':
        return 'FASTA'
    else:
        raise ValueError('File is not FASTA')

def load_fasta(filename):
    if get_sequence_file_type(filename) != 'FASTA':
        sys.exit('Error: {} is not FASTA format'.format(filename))
    fasta_seqs = collections.OrderedDict()
    depths, circular = {}, {}
    p = re.compile(r'depth=([\d.]+)')
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs[name.split()[0]] = ''.join(sequence)
                    sequence = []
                name = line[1:]
                short_name = name.split()[0]
                if 'depth=' in name.lower():
                    try:
                        depths[short_name] = float(p.search(name.lower()).group(1))
                    except (ValueError, AttributeError):
                        depths[short_name] = 1.0
                else:
                    depths[short_name] = 1.0
                circular[short_name] = 'circular=true' in name.lower()
            else:
                sequence.append(line)
        if name:
            fasta_seqs[name.split()[0]] = ''.join(sequence)
    return fasta_seqs, depths, circular

def random_chance(chance):
    assert 0.0 <= chance <= 1.0
    return random.random() < chance

def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])
