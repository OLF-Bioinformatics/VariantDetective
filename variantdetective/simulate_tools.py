"""
This module contains miscellaneous functions that are used in various
components of LongReadGenerator.

Portions Copyright (C) 2022 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/philcharron-cfia/LongReadGenerator
Portions Copyright (C) 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread

This file is part of LongReadGenerator.

LongReadGenerator is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

LongReadGenerator is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
LongReadGenerator. If not, see <http://www.gnu.org/licenses/>.
"""

import collections
import contextlib
import gzip
import io
import os
import random
import re
import sys

def bold(text):
    return BOLD + text + END_FORMATTING

BOLD = '\033[1m'
END_FORMATTING = '\033[0m'

@contextlib.contextmanager
def captured_output():
    new_out, new_err = io.StringIO(), io.StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err



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

def float_to_str(v, decimals=1, trim_zeros=False):
    if float(int(v)) == v:
        return str(int(v))
    else:
        formatter = '%.' + str(decimals) + 'f'
        result = formatter % v
        if trim_zeros:
            while result.endswith('0'):
                result = result[:-1]
        return result
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

def get_random_base():
    """
    Returns a random base with 25% probability of each.
    """
    return RANDOM_SEQ_DICT[random.randint(0, 3)]

RANDOM_SEQ_DICT = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

def get_random_sequence(length):
    """
    Returns a random sequence of the given length.
    """
    return ''.join([get_random_base() for _ in range(length)])

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

def print_in_two_columns(l1p1, l2p1, l3p1, l1p2, l2p2, l3p2, output, space_between=6):
    part_1_len = max(len(l1p1), len(l2p1), len(l3p1)) + space_between
    format_str = '{:<' + str(part_1_len) + '}'
    l1p1 = format_str.format(l1p1)
    l2p1 = format_str.format(l2p1)
    l3p1 = format_str.format(l3p1)
    print(l1p1 + l1p2, file=output)
    print(l2p1 + l2p2, file=output)
    print(l3p1 + l3p2, file=output)

def random_chance(chance):
    assert 0.0 <= chance <= 1.0
    return random.random() < chance

def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])
