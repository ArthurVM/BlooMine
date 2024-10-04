#ifndef BLOOMINE_HPP
#define BLOOMINE_HPP

#include <vector>
#include <future>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include "utilities.hpp"
#include "BloomFilter.hpp"
#include "FQread.hpp"
#include "SPscreenutils.hpp"

namespace po = boost::program_options;

po::variables_map parseArgs( int, char** );
void checkArgc( int, char**);
void checkThreadCount( int );
void checkFloats( float, float );
void checkKsize( int );
std::vector<std::string> readFQ( std::string fq_file );
BloomFilter genBF( po::variables_map, int& );
std::vector<std::string> runBM( po::variables_map,
                                std::vector<std::string>,
                                BloomFilter,
                                std::unordered_set<std::string>,
                                const int,
                                const double,
                                const double,
                                const double,
                                const double );
void printParams( po::variables_map);

#endif /* BLOOMINE_HPP */
