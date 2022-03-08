#!/usr/bin/python3

"""
Split multifasta using the biopython SeqIO module.
"""

import sys, os
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

if __name__=="__main__":

    if len(sys.argv) != 2:
        print("Arguments malformed, please use -h to see usage.")
        sys.exit(1)

    elif "-h" in sys.argv:
        print(__doc__)
        print("\tUSAGE: splitMultifasta.py <fasta>")
        sys.exit(0)

    split_multi(sys.argv[1])
