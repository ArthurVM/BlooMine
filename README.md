# BlooMine
A Bloom filter driven read filtering tool for investigating variation of a target locus within read sets.

## Quick Install

BlooMine requires python 3.8 or 3.9, cmake and boost. Installing these is pretty straightforward.

Install cmake and boost

```
sudo apt-get install cmake libboost-all-dev
```

Compile BlooMine binaries

```
bash install.sh
```

Install the BlooMine package

```
pip install ./
```

## Quick Start
```
BlooMine ./my_targets.fa ./my_fastq_dir/ -t 10 --suffix _{1,2}.fastq.gz
```
This runs BlooMine using targets found in `my_targets.fa` on FASTQ file pairs found in `./my_fastq_dir/` with the pairing suffix `_1.fastq.gz` and `_2.fastq.gz`.

```
usage: BlooMine [-h] [-v] [-k KMER] [-f FALSE_POSITIVE] [-s FP_SIM] [-e SP_ERROR] [-t THREADS]
                [-o OUTDIR] [--suffix SUFFIX] [--on_disk] [--bloomine_exec_bin BLOOMINE_EXEC_BIN]
                target_fasta indir

positional arguments:
  target_fasta          A FASTA/MULTIFASTA file containing target sequence flanks.
  indir                 Directory containing paired FASTQ files to screen for the target sequence.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -k KMER, --kmer KMER  Kmer size. Default=7
  -f FALSE_POSITIVE, --false_positive FALSE_POSITIVE
                        False positive rate for building the Bloom filter. Range=0-1. Default=0.0001
  -s FP_SIM, --FP_sim FP_SIM
                        FP-screen threshold for gene orthology inference as a percentage of kmer array
                        identity. Range=0-100. Default=50.0
  -e SP_ERROR, --SP_error SP_ERROR
                        SP-screen screening error threshold for alignment, where the maximum number of
                        errors to return a read as a hit is 1/n given some n. Default=4.0
  -t THREADS, --threads THREADS
                        Number of threads to use. Must be an even <int>. Default=4.
  -o OUTDIR, --outdir OUTDIR
                        Directory for writing output files. Default=./BlooMine_RUNDIR.
  --suffix SUFFIX       The fastq suffix to use for pairing read files. Default=_{1,2}.fastq.gz
  --on_disk             Write read partitions to disk and stream from them. Low memory usage but
                        slower. Default=False.
  --bloomine_exec_bin BLOOMINE_EXEC_BIN
                        Path to a BlooMine_exec binary for use in this run. By default it will use the
                        executable in $PATH.
```
BlooMine takes a multifasta containing target read sequences formatted like:
```
>target_1 | flank1
...TACG...
>target_1 | flank2
...TACG...
>target_2 | flank1
...TACG...
>target_2 | flank2
...TACG...
```
Sequences must be idenfiable as being flank1 or flank2 of a target sequence.

A directory containing Illumina reads to be screened should be specified, along with a wildcard suffix for pairing read files. For example in the directory
```
.
└── mydata
    ├── cp1_1.fq.gz
    ├── cp1_2.fq.gz
    ├── cp2_1.fq.gz
    ├── cp2_2.fq.gz
    ├── cp3_1.fq.gz
    └── cp3_2.fq.gz
```
the suffix would be `_{1,2}.fq.gz`.

By default, BlooMine runs on physical memory. However, this may require a supercomputing cluster for larger datasets. Consequently, the `--on_disk` flag partitions read files to temporary files and streams them from disk. This is much slower, but requires substantially less memory.
