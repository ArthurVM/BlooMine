"""
Afanc high-resolution Metagenomics disambiguator.
"""
import sys
import argparse
from ._version import __version__

from Afanc.utilities.generalUtils import isDir, isFile, checkDate, reformat_mapping_arg

"""
Parse arguments for BlooMine
"""

def run(args):
    ## level 0 run function
    runBlooMine(args)


base_parser = argparse.ArgumentParser(add_help=True, description=__doc__)

base_parser.add_argument('-v', '--version', action='version',
    version='%(prog)s {version}'.format(version=__version__))

subparsers = base_parser.add_subparsers(
    title="[sub-commands]", dest="command"
)

## bloomine args parser
parser_bloomine = subparsers.add_parser()

parser_getDataset.add_argument('target_fasta',
    type=isFile,
    action='store',
    help='file containing target sequences. Supported formats: FASTA')

parser_getDataset.add_argument('fastq',
    type=isFile,
    action='store',
    help='file containing reads to be screened for containing taret sequences. Supported formats: FASTQ')

parser_getDataset.add_argument('-p', '--prefix',
    action='store',
    default='BM',
    help='prefix for output files. Default=BM')

parser_getDataset.add_argument('-k', '--kmer',
    type=int,
    default=7,
    action='store',
    help='kmer size. Default = 7')

parser_getDataset.add_argument('-f', '--false_positive',
    type=float,
    default=0.0001,
    action='store',
    help='false positive rate for building the Bloom filter. Range = 0-1. Default = 0.0001')

parser_getDataset.add_argument('-s', '--FP_sim',
    type=float,
    default=50.0,
    action='store',
    help='FP-screen threshold for gene orthology inference as a percentage of kmer array identity. Range = 0-100. Default = 50.0')

parser_getDataset.add_argument('-e', '--SP_error',
    type=float,
    default=4.0,
    action='store',
    help='SP-screen screening error threshold for alignment, where the maximum number of errors to return a read as a hit is 1/n given n. Default = 4.0')

parser_bloomine.set_defaults(func=run)
