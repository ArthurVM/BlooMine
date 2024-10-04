#ifndef BLOOMINE_UTILS_HPP
#define BLOOMINE_UTILS_HPP

#include <vector>
#include <future>
#include <iostream>
#include <fstream>
#include "utilities.hpp"
#include "constants.hpp"
#include "BloomFilter.hpp"
#include "FQread.hpp"
#include "FastQ.hpp"
#include "SPscreenutils.hpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;


// Function Declarations
po::variables_map parseArgs( int argc, char** argv );
void checkArgc( int argc, char** argv );
void checkThreadCount( int threads );
void checkFloats( float f, float s );
void checkKsize( int k );
void printParams( po::variables_map args );

// Bloom Filter and Threshold Generation
std::pair<BloomFilter, int> generateBloomFilter( std::string target_fasta, 
                                                 int kmer, 
                                                 float first_pass_similarity, 
                                                 float fp_rate );

// Minimum Score Threshold Calculation
double calculateMinimumScoreThreshold( float second_pass_error,
                                       int kmer, 
                                       std::unordered_set<std::string> target_kset );

// Data Partitioning into meory
std::vector<std::vector<std::string>> partitionData( std::string fastq, 
                                                     int threads );

// Data Partitioning into disk
FastQ partitionData( std::string fastq, 
                     int threads,
                     bool on_disk );

// FASTQ File Reading
std::vector<std::string> readFQ( std::string fastq );

// Thread Spawning from memory
std::vector<std::future<std::vector<std::string>>> spawnThreads( int threads,
                                                                 int kmer,
                                                                 std::string prefix,
                                                                 std::vector<std::vector<std::string>> read_data_parts, 
                                                                 BloomFilter BF,
                                                                 std::unordered_set<std::string> target_kset,
                                                                 int fp_threshold,
                                                                 double spMST );

// Thread Spawning from disk
std::vector<std::future<std::vector<std::string>>> spawnThreads( int threads,
                                                                 int kmer,
                                                                 std::string prefix,
                                                                 FastQ FQ, 
                                                                 BloomFilter BF,
                                                                 std::unordered_set<std::string> target_kset,
                                                                 int fp_threshold,
                                                                 double spMST );

// Run BlooMine from memory
std::vector<std::string> runBM( int kmer,
                                std::vector<std::string>,
                                BloomFilter,
                                std::unordered_set<std::string>,
                                int,
                                double );

// Run BlooMine from disk
std::vector<std::string> runBMdisk( int kmer,
                                    std::string fq_file,
                                    BloomFilter BF,
                                    std::unordered_set<std::string> target_kset,
                                    int fp_threshold,
                                    double spMST );

// Result Collection
// void collectResults( VariablesMap, std::vector<std::future<std::vector<std::string>>> );

#endif /* BLOOMINE_UTILS_HPP */
