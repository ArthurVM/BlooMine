#ifndef SPSCREEN
#define SPSCREEN

#include "utilities.hpp"
#include "BloomFilter.hpp"
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

std::unordered_set<std::string> genkmerset( std::string, int );

class BlooMineSP
{
  /* Class to control execution of an instance of the BlooMineSP pipeline for easier multithreading
  */



  public:
    BlooMineSP( po::variables_map, std::string ); // BlooMineSP constructor
    ~BlooMineSP();                                // BlooMineSP destructor

  private:
    std::string _tfasta;                          // input target fasta and fastq path strings
    float _min;                                   // minimum fraction of unique kmers to call a hit
    int _k;                                       // kmer size and threshold score
    std::vector<std::string> _FPhits;             // vector containing reads hits from first-pass screening
    std::vector<std::string*> _SPhits;            // vector containing reads hits from second pass screening as pointers to elements of _FPhits

};

BlooMineSP::BlooMineSP( po::variables_map vm, std::string fq_part )
{
  /* BlooMineSP constructor.
  where:
    vm (variables_map) : a boost program options variables map
  */
}

BlooMineSP::~BlooMineSP(void) { /* BlooMineSP class destructor */ }


#endif /* SPSCREEN */
