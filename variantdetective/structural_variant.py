import os
import sys
from subprocess import Popen, PIPE, STDOUT

from .tools import get_new_filename

def structural_variant(args, input_reads, output=sys.stderr):
    print('Running structural_variant tool...', file=output)
    reference = get_new_filename(args.reference, args.out)   
    structural_variant_outdir = os.path.join(args.out, 'structural_variant')
    if not os.path.isdir(structural_variant_outdir):
        os.makedirs(structural_variant_outdir)

    # Run Nanovar 
    command = 'nanovar -t ' + str(args.threads) + ' ' + input_reads  + ' ' + reference + ' ' + structural_variant_outdir
    print(command)
    process = Popen([command],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT, shell=True, executable='/bin/bash')
    exitcode = process.wait()
    if exitcode != 0:
        raise Exception("Error: Nanovar failed")
