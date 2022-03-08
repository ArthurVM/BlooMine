#!/usr/bin/python3

import sys, os
import argparse
import subprocess
from Bio import SeqIO

def split_multi(fasta):
    """ Takes a multifasta and splits it into individual fasta files.
    """

    seqrecs = list(SeqIO.parse(fasta, "fasta"))
    if len(seqrecs) > 2:
        print(f"More than 2 sequences in {fasta}, exiting...")
        sys.exit(1)

    elif len(seqrecs) == 1:
        f1_out = "./tmp.f1.fa"
        SeqIO.write(seqrecs, f1_out, "fasta")

    elif len(seqrecs) == 2:
        f1_out = "./tmp.f1.fa"
        f2_out = "./tmp.f2.fa"

        SeqIO.write(seqrecs[0], f1_out, "fasta")
        SeqIO.write(seqrecs[1], f2_out, "fasta")

def runBlooMine(args):
    """ Run BlooMine using both flanks
    """
    split_multi(args.target_fasta)

    f1_runline = f"{os.path.dirname(args.script_path)}/bin/BlooMine\
     -t ./tmp.f1.fa \
     -q {','.join(args.fastq)} \
     --prefix BM_F1_screen \
     --kmer {args.kmer} \
     --false_positive {args.false_positive} \
     --FP-sim {args.FP_sim} \
     --SP-error {args.SP_error} \
     --threads {args.threads}"

    subprocess.run(f1_runline, shell=True)

    f2_runline = f"{os.path.dirname(args.script_path)}/bin/BlooMine\
      -t ./tmp.f2.fa \
      -q BM_F1_screen_BMfiltered.fq \
      --prefix {args.prefix} \
      --kmer {args.kmer} \
      --false_positive {args.false_positive} \
      --FP-sim {args.FP_sim} \
      --SP-error {args.SP_error} \
      --threads {args.threads}"

    subprocess.run(f2_runline, shell=True)

    moi_runline = f"python3 {os.path.dirname(args.script_path)}/scripts/MOIscreen.py {args.target_fasta} {args.prefix}_BMfiltered.fq"
    subprocess.run(moi_runline, shell=True)

def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('target_fasta', action='store', help='A FASTA file containing target sequence flanks.')
    parser.add_argument('fastq', action='store', nargs="*", help='Paired or unpaired FASTQ files to screen for the target sequence.')

    parser.add_argument('-p', '--prefix', action='store', default="BM", help='Prefix for output files. Default=BM.')
    parser.add_argument('-k', '--kmer', action='store', default=7, help='Kmer size. Default=7')
    parser.add_argument('-f', '--false_positive', action='store', default=0.0001, help='False positive rate for building the Bloom filter. Range=0-1. Default=0.0001')
    parser.add_argument('-s', '--FP_sim', action='store', default=50.0, help='FP-screen threshold for gene orthology inference as a percentage of kmer array identity. Range=0-100. Default=50.0')
    parser.add_argument('-e', '--SP_error', action='store', default=4.0, help='SP-screen screening error threshold foralignment, where the maximum number of errors to return a read as a hit is 1/n given some n. Default=4.0.')
    parser.add_argument('-t', '--threads', action='store', default=4, help='Number of threads to use. Must be an even <int>. Default=4.')


    args = parser.parse_args(argv)
    return args

def main(argv):
    """ BlooMine main function
    """
    args = parseArgs(argv)
    runBlooMine(args)

if __name__=="__main__":
    main(sys.argv)
