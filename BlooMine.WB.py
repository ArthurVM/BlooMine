#!/usr/bin/python3

import sys, os
import argparse
import subprocess
from collections import defaultdict
from Bio import SeqIO

def splitMulti(fasta, single_flank):
    """ Takes a multifasta and splits it into individual fasta files.
    """

    target_fastas = []
    seqrecs = list(SeqIO.parse(fasta, "fasta"))

    if not single_flank:
        fasta_pair_dict = pairRecords(seqrecs)
        print(f"{len(fasta_pair_dict)} flank pairs detected.")

        for id, recs in fasta_pair_dict.items():

            f1_out = f"./tmp.{recs[0].id}_f1.fa"
            SeqIO.write(f1_out, recs[0], "fasta")

            f2_out = f"./tmp.{recs[1].id}_f2.fa"
            SeqIO.write(f2_out, recs[1], "fasta")

            target_fastas.append([id, f1_out, f2_out])

        return target_fastas

    else:
        print(f"{len(seqrecs)} single end flanks detected.")

        for rec in seqrecs:

            f_out = f"./tmp.{rec.id}.fa"
            SeqIO.write(f_out, rec, "fasta")

            f_out.append([rec.id, f_out])

        return target_fastas

def pairRecords(seqrecs):
    """ Pair records in a multifasta according to their ID.
    Records should take the form of:

    >interesting_target_1 | flank_1
    TTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAAAAAAAGGGGGGGGGG
    >interesting_target_1 | flank_2
    TTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAAAAAAAGGGGGGGGGG
    >interesting_target_2 | flank_1
    TTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAAAAAAAGGGGGGGGGG
    >interesting_target_2 | flank_2
    TTTTTTTTTTTTCCCCCCCCCCCCCCAAAAAAAAAAAAAAAGGGGGGGGGG
    ...
    """

    fasta_pair_dict = defaultdict(list)
    for rec in seqrecs:
        fasta_pair_dict[rec.id].append(rec)

    return fasta_pair_dict

def run(args):
    """ Run BlooMine
    """

    target_fastas = splitMulti(args.target_fasta, args.single_flanks)

    ## if single flank mode is not selected, as by default, run on the second flank
    if not args.single_flank:

        for id, f1, f2 in target_fastas:
            BlooMine(args, f1, ",".join(args.fastq), f"{prefix}.{id}_F1screen")
            BlooMine(args, f2, f"{prefix}.{id}_FPscreen_BMfiltered.fq", f"{prefix}.{id}_SPscreen")

            ## MOI detection can only be done in instances where there are multiple flanks

            moi_runline = f"python3 {os.path.dirname(args.script_path)}/scripts/MOIscreen.py {args.target_fasta} {args.prefix}_BMfiltered.fq"
            subprocess.run(moi_runline, shell=True)

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

def BlooMine(args, target, query, prefix):
    """ Function for running BlooMine
    """

    runline  = f"{os.path.dirname(args.script_path)}/bin/BlooMine\
     -t ./{target} \
     -q {query} \
     --prefix {prefix} \
     --kmer {args.kmer} \
     --false_positive {args.false_positive} \
     --FP-sim {args.FP_sim} \
     --SP-error {args.SP_error} \
     --threads {args.threads}"

     subprocess.run(runline, shell=True)

def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('target_fasta', action='store', help='A FASTA/MULTIFASTA file containing target sequence flanks.')
    parser.add_argument('fastq', action='store', nargs="*", help='Paired or unpaired FASTQ files to screen for the target sequence.')

    parser.add_argument('-p', '--prefix', action='store', default="BM", help='Prefix for output files. Default=BM.')
    parser.add_argument('-k', '--kmer', action='store', default=7, help='Kmer size. Default=7')
    parser.add_argument('-f', '--false_positive', action='store', default=0.0001, help='False positive rate for building the Bloom filter. Range=0-1. Default=0.0001')
    parser.add_argument('-s', '--FP_sim', action='store', default=50.0, help='FP-screen threshold for gene orthology inference as a percentage of kmer array identity. Range=0-100. Default=50.0')
    parser.add_argument('-e', '--SP_error', action='store', default=4.0, help='SP-screen screening error threshold foralignment, where the maximum number of errors to return a read as a hit is 1/n given some n. Default=4.0.')
    parser.add_argument('--single_flank', action='store_true', default=False, help='Run in single flank mode, where only 1 flanking region is provided. MOI detection is not run. Default=False.')
    parser.add_argument('-t', '--threads', action='store', default=4, help='Number of threads to use. Must be an even <int>. Default=4.')


    args = parser.parse_args(argv)
    return args

def main(argv):
    """ BlooMine main function
    """
    args = parseArgs(argv)
    run(args)

if __name__=="__main__":
    main(sys.argv)
