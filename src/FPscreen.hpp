#ifndef FPSCREEN
#define FPSCREEN

#include "utilities.hpp"
#include "BloomFilter.hpp"
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

std::unordered_set<std::string> genkmerset( std::string, int );

class BlooMineFP
{
  /* Class to control exectuion of an instance of the BlooMineFP pipeline for easier multithreading
  */
  const po::variables_map vm;     // boost po variables map
  const std::string fq_part;      // input fastq file string

  public:
    std::string _tfasta, _fastq;                  // input target fasta and fastq path strings
    float _sim, _fp;                              // similarity and false positive rate
    int _k, _t;                                   // kmer size and threshold score
    std::vector<std::string> _FPhits;             // vector containing read hits from first pass screening
    BlooMineFP( po::variables_map, const std::string ); // BlooMineFP constructor
    ~BlooMineFP();                                // BlooMineFP destructor
    int runFPscreen();                            // run first pass screen

  private:
    std::unordered_set<std::string> fasta2kmerarray();
    void screen_fastq( BloomFilter );
    bool screen_read( std::string, BloomFilter );
};

BlooMineFP::BlooMineFP( const po::variables_map vm, const std::string fq_part )
{
  /* BlooMineFP constructor.
  where:
    vm (variables_map) : a boost program options variables map
  */
  _tfasta = vm["target_fasta"].as<std::string>();
  _fastq = fq_part;
  _sim = vm["similarity"].as<float>();
  _fp = vm["false_positive"].as<float>();
  _k = vm["kmer"].as<int>();
}

BlooMineFP::~BlooMineFP(void) { /* BlooMineFP class destructor */ }

int BlooMineFP::runFPscreen()
{
  /* Public function which controls the execution of the BlooMineFP first-pass screening step
  */
  std::unordered_set<std::string> kmer_array;

  try {
    kmer_array = BlooMineFP::fasta2kmerarray();
  }
  catch ( const std::out_of_range& oor ) {
    vprint( "FA-ERROR", "Generation of kmer array from " + _tfasta + " failed. Exiting...", "r" );
    exit(1);
  }

  _t = kmer_array.size() * _sim/100;     // define the minimum number of kmers that much match to return a read hit
  // cout << "threshold : " << threshold << endl;         // uncomment for debugging

  BloomFilter BF(kmer_array.size(), _fp);

  for ( auto & kmer : kmer_array ) { BF.push(kmer); } // iterates through the target kmer array and pushes each kmer into the Bloom Filters

  // start first pass screening
  screen_fastq( BF );
  return 0;
}

std::unordered_set<std::string> BlooMineFP::fasta2kmerarray()
{
  // N.B. This will only return the final sequence, needs fixing.

  std::unordered_set<std::string> kmer_array;
  std::ifstream infile(_tfasta);

  try {
    if ( !infile.good() ) {
      throw _tfasta;
    }
  }
  catch ( std::string _tfasta ) {
    vprint("FA-ERROR", "Cannot open FASTA file : " + _tfasta, "r" );
    exit(1);
  }

 std::string line, id, sequence;

  while ( std::getline(infile, line) ) {
    if ( line.empty() ) { continue; }

    if ( line[0] == '>' ) {
      // output previous line before overwriting id if it is non-empty
      if ( !id.empty() ) { kmer_array = genkmerset( sequence, _k ); }

      id = line.substr(1);
      sequence.clear();
    }

    else { sequence += line; }
  }

  // output final entry if id is non-empty
  if ( !id.empty() ) { kmer_array = genkmerset( sequence, _k ); }

  infile.close();

  return kmer_array;
}

void BlooMineFP::screen_fastq( BloomFilter BF )
{
  std::int32_t linecounter = 0;

  std::ifstream infile(_fastq);

  if ( !infile.good() ) {
    vprint( "FQ-ERROR", "There was a problem opening the target FASTQ file. Exiting...", "r" );
    exit (2);
  }

  std::string line, id, seq, qual, desc;   // declaire fastq line data
  std::vector<std::string> read;           // read vector

  while ( getline(infile, line) ) {
    linecounter++;    // increment the line counter

    if ( linecounter % 4 != 0 ) { read.push_back(line); }   // check if the linecounter has reaches a new block

    else {
      read.push_back(line);    // capture the description line

      // std::transform(read[1].begin(), read[1].end(), read[1].begin(), std::ptr_fun<int, int>(std::toupper));

      if ( screen_read( read[1], BF ) ) {
        // if the sequence is a hit, concatenate into fastq format and push to the FPhits vector
        std::string read_concat;
        for ( const auto & l : read ) { read_concat += l + "\n"; }
        _FPhits.push_back(read_concat);
      }

      read.clear(); // redefine the vector as empty
    }
  }
  infile.close();
}

bool BlooMineFP::screen_read( std::string seq, BloomFilter BF )
{
  /* Screen the read against the target sequence Bloom filter
  */
  int hit_count = 0;
  std::unordered_set<std::string> kset = genkmerset( seq, _k );

  for ( auto & kmer : kset ) {
    if ( BF.check(kmer) ) {
      hit_count++;
    }
  }
  if ( hit_count >= _t ) { return true; }
  else { return false; }
}

#endif /* FPSCREEN */
