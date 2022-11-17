/* A set of generic utility functions for BlooMine.
A. V. Morris 2020
*/

#ifndef BMUTILITIES_HPP
#define BMUTILITIES_HPP

#include <iostream>
#include <functional>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <chrono>

#include "colours.hpp"

using std::cout;
using std::cerr;
using std::endl;
using namespace std::chrono;

std::unordered_set<std::string> genkmerset( std::string, const int );
std::unordered_set<std::string> fasta2kmerset( std::string, int );
int vprint( std::string, std::string, std::string );

void printNtimes( char c, int n )
{
  // character c will be printed n times
  cout << std::string(n, c);
}

std::unordered_set<std::string> genkmerset( std::string seq, const int k )
{
  /* takes a sequence and a value of k and generates a kmer set
  N.B. this is different to a kmer array, in that repeated elements are collapsed and no counter is provided
  */

  std::unordered_set<std::string> kmer_set;

  for ( long unsigned int i = 0; i < seq.size()-k+1; i++ ) {
    std::string kmer = seq.substr(i, k);

    kmer_set.insert(kmer);
  }
  return kmer_set;
}

std::unordered_set<std::string> fasta2kmerset( std::string fasta, int k )
{
  // N.B. This will only return the final sequence, needs fixing.

  std::unordered_set<std::string> target_kset;
  std::ifstream infile(fasta);

  try {
    if ( !infile.good() ) {
      throw fasta;
    }
  }
  catch ( std::string fasta ) {
    vprint("FA-ERROR", "Cannot open FASTA file : " + fasta, "r" );
    exit(1);
  }

 std::string line, id, sequence;

  while ( std::getline(infile, line) ) {
    if ( line.empty() ) { continue; }

    if ( line[0] == '>' ) {
      // output previous line before overwriting id if it is non-empty
      if ( !id.empty() ) { target_kset = genkmerset( sequence, k ); }

      id = line.substr(1);
      sequence.clear();
    } else { sequence += line; }
  }

  // output final entry if id is non-empty
  if ( !id.empty() ) { target_kset = genkmerset( sequence, k ); }

  try {
    if ( target_kset.size() < (sequence.size()-k+1)*0.5 ) {
      throw target_kset.size();
    }
  }
  catch ( size_t ka_size ) {
    std::string warning_msg = "The target sequence you have selected is low complexity (less than hald kmers unique) :\n\t\tpattern length: " + std::to_string(sequence.size()) +
      "\n\t\tunique kmers: " + std::to_string(ka_size) +
      "\n\t\tThis may result in false positive results.\n";
    vprint("WARNING", warning_msg, "y");
  }

  infile.close();

  return target_kset;
}

int vprint( std::string subprocess, std::string message, std::string col )
{
  /* prints process output in a nice pretty (useful) way
  */
  time_t my_time = time(0);
  tm* dt = localtime(&my_time);

  // ctime() used to give the present time
  cout << " " << dt->tm_hour << ":" << dt->tm_min << ":" << dt->tm_sec << "\t";

  if ( col == "r") {
    cout << "\x1B[31m" + subprocess;
  } else if ( col == "b") {
    cout << "\x1B[33m" + subprocess;
  } else if ( col == "g") {
    cout << "\x1B[32m" + subprocess;
  } else {
    cout << subprocess;
  }

  std::cout << "\x1B[0m  ::  " << message << "\n\n";

  return 0;
}

char complement(char n)
{
    switch(n)
    {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    default:
      return 'N';
    }
}

std::string reverse_complement( std::string seq )
{
  /* Generates the reverse complement of a DNA sequence
  */
  std::string cseq;

  for ( auto it = seq.crbegin(); it != seq.crend(); it++ ) {
    try { cseq += complement(*it); }   // complements the base and adds to cseq
    catch ( std::out_of_range& e) { cseq += "N"; }          // if it is an unknown base, add N to cseq
  }
  return cseq;
}

#endif /* BMUTILITIES_HPP */
