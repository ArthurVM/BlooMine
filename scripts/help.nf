def printHelp() {
  log.info"""
  Usage:
    nextflow run BlooMine.nf target fastq --k [int] --threads [int] --outdir [str] (--MOI)

  Description:
    Screen a fastq read set using a pair of target flanking sequences to identify all reads which contain the target sequence.
      - MOI: Run the Multiplicity of Infection pipeline to generate a report of variation at the target locus exhibited in the
        read set.

  Mandatory Arguments:
    --target            Sequences which flank the target region to identify within the read set. Alternatively this may be the target
                        sequence itself.
    --fastq             Paired or unpaired reads in fastq format. Paired reads must be split into 2 files, representing the forward and backward
                        read sets.
    --prefix            A string prefix for output files.

  Optional Arguments:
    --k                 The kmer size to use for screening. Default=7.
    --FPsim             First pass screen threshold for sequence similarity, as a percentage kmer array intersection. Default=50.0.
    --SPerror           Second pass screen error threshold for alignment, where the max number of acceptable errors is 1/n. Default n=4.
    --threads           The number of threads. Default=1.
    --outdir            Output directory for results. Default=./results.
  """.stripIndent()
}
