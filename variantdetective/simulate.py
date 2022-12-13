"""
This file contains code required to run the long read simulation of VariantDetective.

Portions Copyright (C) 2022 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/philcharron-cfia/VariantDetective
Portions Copyright (C) 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread
"""

import multiprocessing
import os 
import random
import sys
import uuid

from .simulate_tools import load_fasta, reverse_complement, random_chance
from .fragment_lengths import FragmentLengths
from .version import __version__

def simulate(args, input_fasta, output=sys.stderr):
    split_path = os.path.splitext(input_fasta)
    if split_path[1] == ".gz":
        output_name = os.path.splitext(split_path[0])[0] + ".fastq"
    else:
        output_name = split_path[0] + ".fastq"
    
    print(f'Simulating long-reads from genomic sequence...', end=' ', file=output)
    ref_seqs, ref_depths, ref_circular = load_reference(input_fasta, output)
    rev_comp_ref_seqs = {name: reverse_complement(seq) for name, seq in ref_seqs.items()}
    frag_lengths = FragmentLengths(args.mean_frag_length, args.frag_length_stdev, output)
    adjust_depths(ref_seqs, ref_depths, ref_circular, frag_lengths, args)   
    ref_contigs, ref_contig_weights = get_ref_contig_weights(ref_seqs, ref_depths)    
    ref_size = sum(len(x) for x in ref_seqs.values())
    target_size = get_target_size(ref_size, args.readcov)
    #process_size = target_size / args.threads
    process_size = target_size 
    #print('', file=output)
    #print(f'Target read set size: {target_size:,} bp', file=output)
    #print(f'Number of threads: {args.threads}', file=output)
    #print('', file=output)
   
    output_file = open(output_name, 'w')

    process_list = []
    for i in range(1):
        p =  multiprocessing.Process(target= generate_reads,
                                     args = [process_size, frag_lengths, ref_seqs,
                                     rev_comp_ref_seqs, ref_contigs,
                                     ref_contig_weights, ref_circular, output_file])
        p.start()
        process_list.append(p)
    for process in process_list:
        process.join()

    output_file.close()
    print('Complete', file=output)
    return output_name

def generate_reads(target_size, frag_lengths, ref_seqs, rev_comp_ref_seqs,
                   ref_contigs, ref_contig_weights, ref_circular, output_file):
    total_size = 0
    count = 0

    while total_size < target_size:
        fragment, info = build_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs, ref_contigs,
                                        ref_contig_weights, ref_circular)     
        quals = 'S'*len(fragment)

        if len(fragment) == 0:
            continue

        info.append(f'length={len(fragment)}')

        read_name = uuid.UUID(int=random.getrandbits(128))
        info = ' '.join(info)
        print(f'@{read_name} {info}\n'+ fragment + '\n+\n' + quals, file=output_file)

     
        total_size += len(fragment)
        count += 1  

def adjust_depths(ref_seqs, ref_depths, ref_circular, frag_lengths, args):
    sampled_lengths = [frag_lengths.get_fragment_length() for x in range(100000)]
    total = sum(sampled_lengths)
    for ref_name, ref_seq in ref_seqs.items():
        ref_len = len(ref_seq)
        ref_circ = ref_circular[ref_name]

        # Circular plasmids may have to have their depth increased due compensate for misses.
        if ref_circ:
            passing_total = sum(length for length in sampled_lengths if length <= ref_len)
            if passing_total == 0:
                sys.exit('Error: fragment length distribution incompatible with reference lengths.')
            adjustment = total / passing_total
            ref_depths[ref_name] *= adjustment

        # Linear plasmids may have to have their depth increased due compensate for truncations.
        if not ref_circ:
            passing_total = sum(min(ref_len, length) for length in sampled_lengths)
            adjustment = total / passing_total
            ref_depths[ref_name] *= adjustment

def build_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs, ref_contigs, ref_contig_weights,
                   ref_circular):
    info = []
    frag_seq, frag_info = get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs,
                                       ref_contigs, ref_contig_weights, ref_circular)
    info.append(','.join(frag_info))
    return frag_seq, info

def get_fragment(frag_lengths, ref_seqs, rev_comp_ref_seqs, ref_contigs, ref_contig_weights,
                 ref_circular):
    fragment_length = frag_lengths.get_fragment_length()
 
    # The get_real_fragment function can return nothing so we try repeatedly
    # until we get a result.
    for _ in range(1000):
        seq, info = get_real_fragment(fragment_length, ref_seqs, rev_comp_ref_seqs, ref_contigs,
                                      ref_contig_weights, ref_circular)
        if seq != '':
            return seq, info
    sys.exit('Error: failed to generate any sequence fragments - are your read lengths '
             'incompatible with your reference contig lengths?')

def get_real_fragment(fragment_length, ref_seqs, rev_comp_ref_seqs, ref_contigs,
                      ref_contig_weights, ref_circular):
    if len(ref_contigs) == 1:
        contig = ref_contigs[0]
    else:
        contig = random.choices(ref_contigs, weights=ref_contig_weights)[0]
    info = [contig]
    if random_chance(0.5):
        seq = ref_seqs[contig]
        info.append('+strand')
    else:
        seq = rev_comp_ref_seqs[contig]
        info.append('-strand')

    # If the reference contig is linear and the fragment length is long enough, then we just
    # return the entire fragment, start to end.
    if fragment_length >= len(seq) and not ref_circular[contig]:
        info.append('0-' + str(len(seq)))
        return seq, info

    # If the reference contig is circular and the fragment length is too long, then we fail to get
    # the read.
    if fragment_length > len(seq) and ref_circular[contig]:
        return '', ''

    start_pos = random.randint(0, len(seq)-1)
    end_pos = start_pos + fragment_length

    info.append(f'{start_pos}-{end_pos}')

    # For circular contigs, we may have to loop the read around the contig.
    if ref_circular[contig]:
        if end_pos <= len(seq):
            return seq[start_pos:end_pos], info
        else:
            looped_end_pos = end_pos - len(seq)
            assert looped_end_pos > 0
        return seq[start_pos:] + seq[:looped_end_pos], info

    # For linear contigs, we don't care if the ending position is off the end - that will just
    # result in the read ending at the sequence end (and being shorter than the fragment
    # length).
    else:
        return seq[start_pos:end_pos], info

def get_ref_contig_weights(ref_seqs, ref_depths):
    ref_contigs = [x[0] for x in ref_depths.items()]
    ref_contig_weights = [x[1] * len(ref_seqs[x[0]]) for x in ref_depths.items()]
    return ref_contigs, ref_contig_weights

def get_target_size(ref_size, coverage):
    try:
        return int(coverage)
    except ValueError:
        pass
    coverage = coverage.lower()
    try:
        last_char = coverage[-1]
        value = float(coverage[:-1])
        if last_char == 'x':
            return int(round(value * ref_size))
        elif last_char == 'g':
            return int(round(value * 1000000000))
        elif last_char == 'm':
            return int(round(value * 1000000))
        elif last_char == 'k':
            return int(round(value * 1000))
    except (ValueError, IndexError):
        pass
    sys.exit('Error: could not parse coverage\n'
             '--coverage must be either an absolute value (e.g. 250M) or a relative depth '
             '(e.g. 25x)')

def load_reference(reference, output):
    #print('', file=output)
    #print(f'Loading reference from {reference}', file=output)
    ref_seqs, ref_depths, ref_circular = load_fasta(reference)
    plural = '' if len(ref_seqs) == 1 else 's'
    #print(f'  {len(ref_seqs):,} contig{plural}:', file=output)
    for contig in ref_seqs:
        circular_linear = 'circular' if ref_circular[contig] else 'linear'
        #print(f'    {contig}: {len(ref_seqs[contig]):,} bp, {circular_linear}, '
        #      f'{ref_depths[contig]:.2f}x depth', file=output)
    if len(ref_seqs) > 1:
        total_size = sum(len(s) for s in ref_seqs.values())
        #print(f'  total size: {total_size:,} bp', file=output)
    return ref_seqs, ref_depths, ref_circular

def print_intro(output):
    #print('', file=output)
    #print(f'LongReadGenerator v{__version__}', file=output)
    #print(f'Generate fake long reads from genomic sequence', file=output)
    print(f'Simulating long-reads from genomic sequence...', end=' ', file=output)
    
def print_progress(count, bp, target, output):
    plural = ' ' if count == 1 else 's'
    percent = int(1000.0 * bp / target) / 10
    if percent > 100.0:
        percent = 100.0
    #print(f'\rSimulating: {count:,} read{plural}  {bp:,} bp  {percent:.1f}%',
    #      file=output, flush=True, end='')