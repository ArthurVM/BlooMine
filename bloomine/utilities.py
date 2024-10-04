import sys
import re
import glob
import subprocess
from collections import defaultdict
from Bio import SeqIO
from os import path, mkdir, chdir, popen, system


def isFile(filename):
    """ Checks if a path is an existing file """

    if not path.isfile(filename):
        print(f"No file found at {filename}")
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))


def isDir(dirname):
    """ Checks if a path is an existing directory """

    if not path.isdir(dirname):
        print(f"No directory found at {dirname}")
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(dirname)))


def mkchdir(dir, ch=True):
    """ make and change to a dir
    """
    if not path.isdir(dir): mkdir(dir)
    if ch: chdir(dir)


def expandSuffix(suffix):
    """Expands a bash wildcard expression like _{1,2}.fastq.gz to individual strings.

    INPUT:
        suffix (str): The bash wildcard expression to expand.

    OUPUT:
        list: A list of individual strings captured by the expression.
    """

    ## check bash expansion syntax is used
    if "{" not in suffix:
        return [suffix]

    tmp = re.split(r"{|}|,", suffix)
    parts = tmp[1:-1]
    s1 = tmp[0]
    s2 = tmp[-1]

    expanded_suffixes = []
    for part in parts:
        expanded_suffixes.append(s1 + part + s2)
    return expanded_suffixes


def splitMultiFasta(fasta, outdir):
    """ Takes a multifasta containing multiple flanking sequences, and pairs then, producing a fasta containing individual flank pairs.
    """

    flankdict = defaultdict(list)
    flank_pairs = {}

    for rec in SeqIO.parse(fasta, "fasta"):
        flankdict[rec.id].append(rec)

    # print(f"{len(flankdict.keys())} flank headers found.\n")
        
    for id, flanks in flankdict.items():
        if len(flanks) != 2:
            print("""Flank headers malformed! Please ensure flank headers take the form of:

            >target_1 | flank1
            ...
            >target_1 | flank2
            ...
            >target_2 | flank1
            ...
            >target_2 | flank2
            ...""")
            sys.exit(1)

        f1, f2 = flanks

        f1f2_fasta = f"{outdir}/tmp.{f1.id}_flanks.fa"
        with open(f1f2_fasta, "w") as f1f2out:
            f1f2out.write(f">{f1.description}\n{f1.seq}\n>{f2.description}\n{f2.seq}")

        flank_pairs[id] = f1f2_fasta

    return flank_pairs


def groupReads(indir, suffix_list):
    """ Groups reads based on their prefix and suffix.

    INPUT:
        indir (str): The directory containing the reads.
        suffix_list (list): A list of suffixes to group the reads by.

    OUTPUT:
        list: A list of lists, where each inner list contains the prefix and the paths to the reads with the corresponding suffixes.
              like [[prefix, fq1, fq2, ...], ... ]
    """

    base_suffix = suffix_list[0]
    base_box = [f for f in glob.glob(f'{indir}/*{base_suffix}')]
    fq_box = []

    for fq in base_box:

        prefix = path.basename(fq).split(base_suffix)[0]
        tmp = [prefix]

        for end in suffix_list:
            fqpath =  f'{indir}/{prefix}{end}'

            if path.exists(fqpath):
                tmp.append(f'{indir}/{prefix}{end}')
            else:
                print(f"Cannot locate {fqpath}. Please check the input directory and suffix arguments are correct.")
                exit(1)

        fq_box.append(tmp)
    
    if len(fq_box) == 0:
        print(f"Cannot find reads in {indir}! Exiting...")
        exit(1)

    return fq_box


def prepareFQ(fastq, rundir):
    """ Takes a fastq file and checks if it needs decompressing.
    """

    if fastq.endswith(".gz"):
        # print(f"Decompressing {fastq}...")
        fqdc = path.join(rundir, f"tmp.{path.basename(fastq[:-3])}")

        dcline = f"gunzip -c {fastq} > {fqdc}"
        subprocess.run(dcline, shell=True)

        return fqdc
    else:
        return fastq
    

def checkBMexec(args, BLOOMINE_EXEC):
    """ Checks if the BlooMine executable exists and prints its version.

    INPUT:
        BLOOMINE_EXEC (str): The path to the BlooMine executable.

    OUTPUT:
        bool: True if the executable exists and its version is printed, False otherwise.
    """

    try:
        # Try running BlooMine_exec -h to check for basic functionality
        subprocess.check_output([BLOOMINE_EXEC, "-h"])

    except FileNotFoundError:
        print(f"Cannot locate BlooMine executable: {BLOOMINE_EXEC}. Please check the path to the executable or reinstall BlooMine.")
        sys.exit(101)
    except subprocess.CalledProcessError:
        print(f"Error running BlooMine executable: {BLOOMINE_EXEC}. Please check the executable or reinstall BlooMine.")
        sys.exit(101)


def checkRun(prefix, rundir):
    """ Checks whether the MOI report exists.
    """
    return path.exists(f"{rundir}/{prefix}.report.txt")