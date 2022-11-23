import gzip
import os
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

def run_process(command, error_message):
    process = Popen([command],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT, shell=True, executable='/bin/bash')
    exitcode = process.wait()
    if exitcode != 0:
        raise Exception(error_message)
    