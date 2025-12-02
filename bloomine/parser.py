import sys
import argparse

from .utilities import *
from ._version import __version__
from ._splash import __splash__

"""
A Bloom filter driven read mining tool.

BlooMine is a tool for identifying and quantifying the presence of specific target sequences within a metagenomic dataset.
It leverages a combination of Bloom filters and sequence alignment to achieve high sensitivity and specificity,
allowing for the accurate detection of even low-abundance targets.
"""


def run(args):
    ## level 0 run function
    from .run import RunManager

    print(__splash__)

    ## check that the BlooMine executable binary exists
    checkBMexec(args, args.bloomine_exec_bin)

    ## create the base output directory
    mkchdir(args.outdir)

    indir = path.abspath(args.indir)
    fq_box = groupReads(indir, args.suffix)
    
    flank_dict = splitMultiFasta(args.target_fasta, args.outdir)
    
    ## print some info
    print(
        f""" 
        BLOOMINE RUN STARTED

        ## Run Information
        Run directory       : {args.outdir}
        Num. flanks         : {len(flank_dict)}
        Num. FASTQs         : {len(fq_box)}

        ## Arguments        
        Kmer size : {args.kmer}
        Anchor kmer size : {args.min_kmer}
        Skip MOI : {args.skip_moi}
        False positive rate : {args.false_positive}
        FP-screen threshold : {args.FP_sim}
        SP-screen threshold : {args.SP_error}
        Threads             : {args.threads}
        On disk             : {args.on_disk}\n            
        """
        )

    runner = RunManager(args, fq_box, flank_dict, args.outdir)
    runner.run()


## bloomine args parser
parser_bloomine = argparse.ArgumentParser(add_help=True, description=__doc__)

parser_bloomine.add_argument('-v', '--version', action='version',
    version='%(prog)s {version}'.format(version=__version__))

# parser_bloomine.add_argument('script_path',
#     action='store',
#     help=argparse.SUPPRESS)

parser_bloomine.add_argument('target_fasta',
    type=isFile,
    action='store',
    help='A FASTA/MULTIFASTA file containing target sequence flanks.')

parser_bloomine.add_argument('indir',
    action='store',
    type=isDir,
    help='Directory containing paired FASTQ files to screen for the target sequence.')

parser_bloomine.add_argument('-k', '--kmer',
    action='store',
    type=int,
    default=7,
    help='Kmer size. Default=7')

parser_bloomine.add_argument('-f', '--false_positive',
    action='store',
    type=float,
    default=0.0001,
    help='False positive rate for building the Bloom filter. Range=0-1. Default=0.0001')

parser_bloomine.add_argument('-s', '--FP_sim',
    action='store',
    type=float,
    default=35.0,
    help='FP-screen threshold for gene orthology inference as a percentage of kmer array identity. Range=0-100. Default=35.0')

parser_bloomine.add_argument('-e', '--SP_error',
    action='store',
    type=float,
    default=4.0,
    help='SP-screen screening error threshold for alignment, where the maximum number of errors to return a read as a hit is 1/n given some n. Default=4.0')

parser_bloomine.add_argument('-m', '--min_kmer',
    action='store',
    type=int,
    default=11,
    help='The minimum kmer size for use as an anchor during target region extraction. Default=11')

parser_bloomine.add_argument('--skip_moi',
    action='store_true',
    default=False,
    help='Skip MOI analysis. Default=False')

parser_bloomine.add_argument('--polyfamily',
    action='store_true',
    default=False,
    help='Bin reads by highest-scoring probe family and summarise STR variants to JSON. Default=False')

parser_bloomine.add_argument('-t', '--threads',
    action='store',
    default=4,
    help='Number of threads to use. Must be an even <int>. Default=4.')

parser_bloomine.add_argument('-o', '--outdir',
    action='store',
    type=path.abspath,
    default="./BlooMine_RUNDIR",
    help='Directory for writing output files. Default=./BlooMine_RUNDIR.')

parser_bloomine.add_argument('--suffix',
    action='store',
    type=expandSuffix,
    default="_{1,2}.fastq.gz",
    help='The fastq suffix to use for pairing read files. Default=_{1,2}.fastq.gz')

parser_bloomine.add_argument('--on_disk',
    action='store_true',
    default=False,
    help='Write read partitions to disk and stream from them. Low memory usage but slower. Default=False.')

parser_bloomine.add_argument('--bloomine_exec_bin',
    action='store',
    default='BlooMine_exec',
    help='Path to a BlooMine_exec binary for use in this run. By default it will use the executable in $PATH.')

parser_bloomine.set_defaults(func=run)
