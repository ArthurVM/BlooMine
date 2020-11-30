#ifndef BLOOMINE_HPP
#define BLOOMINE_HPP

#include <thread>
#include <future>
#include "program_options.hpp"
#include "utilities.hpp"
#include "BloomFilter.hpp"
#include "FQread.hpp"
#include "FastQ.hpp"
#include "SPscreenutils.hpp"

namespace po = boost::program_options;

po::variables_map parse_args( int, char** );
void check_argc( int, char**);
void check_threadcount( int );
void check_floats( float, float );
void check_ksize( int );
BloomFilter genBF( po::variables_map, int& );
std::vector<std::string> runBM( po::variables_map,
                                std::string,
                                BloomFilter,
                                std::unordered_set<std::string>,
                                const int,
                                const double,
                                const double,
                                const double,
                                const double );
void print_params( po::variables_map);
std::vector<std::string> partition_file( std::string, int );

#endif /* BLOOMINE_HPP */
