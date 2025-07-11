"""
BlooMine MOI screening module.
"""

import sys
import os
import re
import argparse
import datetime
import collections

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def isolate_target(fastq, flanks_fasta, min_kmer):
	""" Isolate the target sequence from the read using the flanking regions.

	INPUT:
	fastq    <string> : path to a file containing reads in FASTQ format
	flanks_fasta   <string> : path to a file containing target sequence flanking regions in FASTA format
	min_kmer <int>    : the minimum kmer size to use in generating kmer lists

	OUTPUT:
	target_dict <dict> : a dictionary containing flanks_fasta sequences that have been isolated from filtered reads
	"""

	target_dict = {}
	flank_coords = {}

	fq_record = SeqIO.parse(fastq, "fastq")

	head_rec, tail_rec = SeqIO.parse(flanks_fasta, "fasta")

	## generate a kascade array for both flanks
	head_kas_array = make_kascade(str(head_rec.seq), min_kmer)
	tail_kas_array = make_kascade(str(tail_rec.seq), min_kmer)

	for record in fq_record:
		## iterate through each read and generate start and end coordinates for the flanks_fasta in the read
		head_pos_kmer_tup = kmer_hit(head_kas_array, record.seq, "head", min_kmer)
		tail_pos_kmer_tup = kmer_hit(tail_kas_array, record.seq, "tail", min_kmer)

		## TODO: work this try catch container out
		head_pos, kmer_h, orientation_h = head_pos_kmer_tup
		tail_pos, kmer_t, orientation_t = tail_pos_kmer_tup

		## if either flank is not found, or if they are found on different strands, skip the read
		if head_pos is None or tail_pos is None or orientation_h != orientation_t:
			continue
		orientation = orientation_h

		## Deal with reverse complemented flanking sequences
		flankRCflag = False
		if head_pos > tail_pos:
			head_pos = len(record.seq) - head_pos + len(head_rec.seq)+1
			tail_pos = len(record.seq) - tail_pos - len(tail_rec.seq)
			flankRCflag = True

		if orientation == "+":
			if flankRCflag:
				target_dict[record.description] = record.seq.reverse_complement()[head_pos:tail_pos]
			else:
				target_dict[record.description] = record.seq[head_pos+1:tail_pos].reverse_complement()

		elif orientation == "-":
			if flankRCflag:
				target_dict[record.description] = record.seq[head_pos:tail_pos]
			else:
				target_dict[record.description] = record.seq.reverse_complement()[head_pos+1:tail_pos].reverse_complement()
				## detects that this is a reverse read, and therefore compliments the sequence to improve ease of analysis
				## then RCs this to ensure comparison against other flanks_fasta regions is on the same strand as the flanking sequences

	return target_dict


def kmer_hit(kas_array, read, flank_flag, min_kmer):
	""" Takes a kascade array and a read and identify the start and end of the flanks_fasta sequence based on
	the hit position of the longest kmer.

	INPUT:
	kmer_array <list>   : an uncounted list of kmers in the given sequence
	read       <string> : an input read
	flank_flag <string> : a flag defining the flank (head/tail)
	min_kmer   <int>    : the minimum kmer size to use in generating kmer lists

	OUTPUT:
	{head/tail}_pos <int>    : the coordinates for the interface between the flank and the flanks_fasta sequence
	kmer            <string> : the longest kmer mapping to the read sequence
	orientation     <string> : the orientation of the mapped kmer
	"""

	len_flank = len((kas_array[0])[0])

	for k_array in kas_array:
		k = len(k_array[0])

		read_array = make_kmer_array(str(read), k)
		comp_read_array = make_kmer_array(str(read.reverse_complement()), k)


		for i, kmer in enumerate(k_array):
			if kmer in read_array:
				orientation = "+"

				if flank_flag == "head":
					head_pos = read_array.index(kmer) + len_flank - i-1
					return head_pos, kmer, orientation

				elif flank_flag == "tail":
					tail_pos = read_array.index(kmer) - i
					return tail_pos, kmer, orientation

			elif kmer in comp_read_array:
				orientation = "-"

				if flank_flag == "head":
					head_pos = comp_read_array.index(kmer) + len_flank - i-1
					return head_pos, kmer, orientation

				elif flank_flag == "tail":
					tail_pos = comp_read_array.index(kmer) - i
					return tail_pos, kmer, orientation
		
	return None, None, None


def report_out(target_dict, fastq, flanks_fasta, out_file="subpop_report.txt"):
	""" Generates a report of flanks_fasta variation in the read set.

	INPUT:
	target_dict <dict>   : a dictionary containing flanks_fasta sequences that have been isolated from filtered reads
	fastq       <string> : path to a file containing reads in FASTQ format
	flanks_fasta      <string> : path to a file containing flanks_fasta sequence flanking regions in FASTA format
	out_file    <string> : path to the output report file
	"""

	date_time = datetime.datetime.now().strftime("%d-%m-%y %H:%M:%S")

	seqvar_counter = collections.Counter([str(value) for key, value in target_dict.items()])
	lenvar_counter = collections.Counter([len(str(value)) for key, value  in target_dict.items()])

	with open(out_file, "w") as rp:

		print(f"Subpop report generated {date_time}\n", file=rp)
		print(f"fastq:\t{fastq}", file=rp)
		print(f"flanks_fasta flanks:\t{flanks_fasta}\n", file=rp)

		print("\nSequence variants:", file=rp)
		for seq, count in seqvar_counter.items():
			print(f"{seq}\t{count}", file=rp)

		print("\nLength variants:", file=rp)
		for size, count in lenvar_counter.items():
			print(f"{size}\t{count}", file=rp)

		print("", file=rp)


def make_kmer_array(seq, k):
	""" Generate a kmer array of seq

	INPUT:
	seq <string> : a DNA sequence
	k   <int>    : a kmer sizes

	OUTPUT:
	kmer_array <list> : an uncounted list of kmers in the given sequence
	"""

	kmer_array = []
	for i in range(len(seq)-k+1):
		kmer = str(seq[i:i+k])
		kmer_array.append(kmer)
	return kmer_array


def make_kascade(seq, min_kmer):
	""" Generates a 'kascade' array for a given sequence, using min_kmer as a lowe bound for k.

	A kascade array is effectively a set of all kmer arrays of a sequence S using a range of k between
	min_kmer and the length of S.

	Given a sequence, S, and a minimum kmer size, m, a kascade array, C, is:
	C^m_S = { K^k_S | m <= k <= |S| }
	where
	K^k_S = { S_(i,j+k) | 0 <= i <= |S|-k+1 }

	INPUT:
	seq <string>   : a DNA sequence
	min_kmer <int> : the minimum kmer size to use in generating kmer lists

	OUTPUT:
	kas_array <list> : a list of kmer lists
	"""

	kas_array = []
	k = len(seq)
	while k >= min_kmer:
		kmer_array = make_kmer_array(seq, k)
		kas_array.append(kmer_array)
		k-=1
	return kas_array