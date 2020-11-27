#!/usr/bin/env nextflow

// nextflow pipeline for BlooMine MOI screening

// enable DSL2
nextflow.enable.dsl=2

// include modules
include {printHelp} from './scripts/help.nf'

params.help = false
params.target = false
params.fastq = false
params.prefix = false
params.FPsim = 50.0
params.SPerror = 4
params.k = 7
params.threads = 4
params.outdir = "results"
params.MOI = false

if (params.help){
    printHelp()
    exit 0
}

if ( !params.target ) {
    println("Please supply a multifasta file containing target flanking sequences.")
    println("Use --help to print help")
    System.exit(1)
}

if ( !params.fastq ) {
    println("Please supply a FASTQ file containing reads to be screened.")
    println("Use --help to print help")
    System.exit(1)
}

if ( !params.prefix ) {
    println("Please supply a prefix for your output files with --prefix.")
    println("Use --help to print help")
    System.exit(1)
} else {
     if ( params.prefix =~ /\// ){
         println("The --prefix that you supplied must not contain the character: \"/\".")
         System.exit(1)
     }
}

println()
println("Target : ${params.target}")
println("FASTQ  : ${params.fastq}")
println("Prefix : ${params.prefix}")
println()

runOutDir = "$PWD/${params.outdir}/${params.prefix}"

process splitMultifasta {
  // make directory structure and split the flank multifasta into temporary files for screening

  output:
  path "f1.tmp.fa", emit: f1
  path "f2.tmp.fa", emit: f2

  """
  if [ ! -d $PWD/${params.outdir} ] ; then
    mkdir $PWD/${params.outdir}
  fi

  if [ ! -d ${runOutDir} ] ; then
    mkdir ${runOutDir}
  fi

  $PWD/scripts/splitMultifasta.py ${params.target}
  """
}

process FirstPassScreen {
  // run first pass screen, using flank 1 as target and the unscreened reads as input fastq

  publishDir "${runOutDir}", pattern: "${params.prefix}*"

  input:
  path flank1

  output:
  path "${params.prefix}_BMfiltered.fq", emit: FPfiltered

  script:
  """
  $PWD/bin/BlooMine -t ${flank1} -q ${params.fastq} --kmer ${params.k} --prefix ${params.prefix} --threads ${params.threads} --FP-sim ${params.FPsim} --SP-error ${params.SPerror} > ${params.prefix}_BM1.log
  mv ${params.prefix}*.log ${runOutDir}/${params.prefix}_FP.log
  """
}

process SecondPassScreen {
  // run second pass screen, using flank 2 as target and the output from first pass screening as input fastq

  publishDir "${runOutDir}", pattern: "${params.prefix}*"

  input:
  path flank2
  path FPout

  output:
  path "${params.prefix}_SP_BMfiltered.fq", emit: SPfiltered

  script:
  """
  $PWD/bin/BlooMine -t ${flank2} -q ${FPout} --prefix ${params.prefix}_SP --threads ${params.threads} --FP-sim ${params.FPsim} --SP-error ${params.SPerror} > ${params.prefix}_BM2.log
  mv ${params.prefix}*.log ${runOutDir}/${params.prefix}_SP.log
  """
}

process MOIscreen {
  // run MOI sceening script to quantify the amount of target variation exhibited in the read set

  publishDir "${runOutDir}", pattern: "${params.prefix}*"

  input:
  path SPout

  output:
  path "${params.prefix}_SP_BMfiltered.report.txt", emit: MOIreport

  script:
  """
  python3 $PWD/scripts/MOIscreen.py ${params.target} ${SPout}
  """
}

workflow {
  // define the workflow

  splitMultifasta()
  FirstPassScreen( splitMultifasta.out.f1.collect() )
  SecondPassScreen( splitMultifasta.out.f2.collect(), FirstPassScreen.out.FPfiltered.collect() )

  // execute MOI screening script if the MOI flag is present
  if (params.MOI){
    MOIscreen( SecondPassScreen.out.SPfiltered.collect() )
  }
}
