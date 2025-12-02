#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP

#include "BlooMineUtils.hpp"
#include "version.hpp"
#include <iostream>
#include <vector>
#include <limits>


/****************************************************************************************************
 * Parses arguments for the BlooMine pipeline.
 *
 * INPUT:
 *  argc <int> : the number of arguments provided
 *  argv <char**> : a ptr-ptr to the arguments
 *
 * OUTPUT:
 *  vm <variables_map> : a boost program options map containing argument variables
 ****************************************************************************************************/
po::variables_map parseArgs( int argc, char** argv ) {

  int threads;
  int k = 0;
  bool disk = 0;
  double f = 0.0;
  double s = 0.0;
  double e = 0.0;
  checkArgc( argc, argv );     // check the number of args provided

  // Declare the supported options.
  po::options_description desc("Allowed options");

  // boost::po args parse
  try {
    desc.add_options()
      ("help,h",          "Show this message and exit.")
      ("target_fasta,t",  po::value<std::string>()->required(),                       "file containing target sequences. Supported formats: FASTA")
      ("fastq,q",         po::value<std::string>()->required(),                       "file containing reads to be screened for containing taret sequences. Supported formats: FASTQ")
      ("prefix",          po::value<std::string>()->default_value("BM")->required(),  "prefix for output files. Default=BM")
      ("kmer",            po::value<int>(&k)->default_value(7)->required(),           "set kmer size. Default = 7")
      ("false_positive",  po::value<double>(&f)->default_value(0.0001)->required(),   "false positive rate for building the Bloom filter. Range = 0-1. Default = 0.0001")
      ("FP-sim",          po::value<double>(&s)->default_value(50.0)->required(),     "FP-screen threshold for gene orthology inference as a percentage of kmer array identity. Range = 0-100. Default = 50.0")
      ("SP-error",        po::value<double>(&e)->default_value(4.0)->required(),      "SP-screen screening error threshold for alignment, where the maximum number of errors to return a read as a hit is 1/n given n. Default = 4.0.")
      ("flank_number",    po::value<int>()->default_value(0),                         "flank identifier for logging (e.g. 1 or 2). Default = 0")
      ("threads",         po::value<int>(&threads)->default_value(4)->required(),     "number of threads to use. Must be an even <int>. Default = 4")
      ("on-disk",                                                                     "write temporary read partitions to disk and stream from them. Slower but low memory consumption. Default = Off");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if ( vm.count("help") ) {
      cout << desc << "\n";
      exit(0);
    }

    po::notify(vm);

    // argscheck chunk
    if ( threads != 1 ) { checkThreadCount( threads ); }  // just returns if running in single core mode
    checkFloats( f, s );  // check the args of the false positive and FP-sim floats
    checkKsize( k );     // check the given kmer size

    // display input arguments
    vprint("FP-ALN", "Started BlooMine with parameters:", "g");
    printParams( vm );

    return vm;
  }

  catch ( const po::error& e ) {
    cerr << "ARG-ERROR: " << e.what() << endl << endl;
    cerr << desc << endl;
    exit (1);
  }
}


/****************************************************************************************************
 * Check the number of arguments provided to BlooMine.
 * Throws an exception (kind of) if too few are provided.
 *
 * INPUT:
 *  argc <int> : the number of arguments provided
 *  argv <char**> : a ptr-ptr to the arguments
 ****************************************************************************************************/
void checkArgc(int argc, char** argv) {
  try {
    if (argc < 2) {
      throw argc;
    }
  } catch (int argc) {
    std::cerr << "Too few argument provided! Please use -h or --help to see arguments." << std::endl;
    exit(1);
  }
}


/****************************************************************************************************
 * Checks the given threadcount is suitable.
 * i.e. threads <= number of available threads, threads = 1, or threads%2 != 0
 * Throws an exception (kind of) if one or more of the above are false.
 *
 * INPUT:
 *  threads <int> : the number of threads requested in arguments
 ****************************************************************************************************/
void checkThreadCount(int threads) {
  const auto processor_count = std::thread::hardware_concurrency();
  try {
    if (threads % 2 != 0 || threads > processor_count) {
      throw threads;
    }
  } catch (int threads) {
    std::string error_msg = "Inappropriate thread count specified. Must be an even <int> of between 1 and " + std::to_string(processor_count) + ". Exiting...";
    vprint("ARG-ERROR", error_msg, "r");
    exit(1);
  }
}


/****************************************************************************************************
 * Checks that the false postitive rate and FP-sim scores are within acceptable range.
 *
 * INPUT:
 *  f <float> : the false positive rate
 *  s <float> : the minimum acceptable similarity score for first pass screening
 ****************************************************************************************************/
void checkFloats(float f, float s) {
  /* Checks that the false postitive rate and FP-sim scores are within acceptable range
  */
  try {
    if (f < 0 || f > 1.0) {
      throw f;
    }
  } catch (float f) {
    vprint("ARG-ERROR", "False-positive rate must be between 0-1. Exiting...", "r");
    exit(1);
  }

  try {
    if (s < 0 || s > 100) {
      throw s;
    }
  } catch (float s) {
    vprint("ARG-ERROR", "Kmer array FP-sim must be between 0-100. Exiting...", "r");
    exit(1);
  }
}


/****************************************************************************************************
 * Checks given kmer size isnt silly.
 *
 * INPUT:
 *  k <int> : the given kmer size
 ****************************************************************************************************/
void checkKsize(int k) {
  /* Checks given kmer size isnt silly
  */
  try {
    if (k < 3) {
      throw k;
    } else if (k >= 3 && k <= 5) {
      vprint("ARG-WARNING", "Kmer size is small. This may give unexpected results. Proceeding with due pessimism...", "y");
    }
  } catch (int k) {
    vprint("ARG-ERROR", "Kmer size too small. Please specific a larger value. Exiting...", "r");
    exit(1);
  }
}


/****************************************************************************************************
 * Prints parameters to stdout during runtime.
 *
 * INPUT:
 *  vm <variables_map> : a boost program options map containing argument variables
 ****************************************************************************************************/
void printParams( po::variables_map vm ) {

  cout << "\t\ttarget_fasta : " << vm["target_fasta"].as<std::string>() << "\n"
  << "\t\tfastq : " << vm["fastq"].as<std::string>() << "\n"
  << "\t\tkmer : " << vm["kmer"].as<int>() << "\n"
  << "\t\tFP : " << vm["false_positive"].as<double>() << "\n"
  << "\t\tFP-sim : " << vm["FP-sim"].as<double>() << "\n"
  << "\t\tSP-error : " << vm["SP-error"].as<double>() << "\n"
  << "\t\tthreads : " << vm["threads"].as<int>() << "\n\n";
}


#endif /* ARGPARSE_HPP */
