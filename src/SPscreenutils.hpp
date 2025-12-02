/* Header file for utilities associated with the BlooMine second-pass screening subprocess.
A. V. Morris 2020
*/

#ifndef BMSPSCREEN_HPP
#define BMSPSCREEN_HPP

#include "utilities.hpp"
#include "BloomFilter.hpp"


/****************************************************************************************************
 * Simple subalignment chunk struct.
 *
 * INPUT:
 *  aln <string> : the alignment chunk string
 *  gap <int>    : the gap count directly downstream to the alignment for aggregating subalignments
 *  score <int>  : the score of the subalignment
 *  start <int>  : start position of the subalignment in the read
 *  end <int>    : end position of the subalignment in the read
 ****************************************************************************************************/
struct SubAln
{
  std::string aln;
  int gap, score, start, end;
};


/****************************************************************************************************
 * Simple read alignment struct.
 *
 * INPUT:
 *  readseq <string>   : the read sequence
 *  alnstring <string> : the maximum scoring subalignment string
 *  score <int>        : the score of the aforementinoed alignment
 ****************************************************************************************************/
struct ReadAln
{
  std::string readseq;
  std::string alnstring;
  int score;
};


/****************************************************************************************************
 * Takes a vector and begin and end indices, and returns the subvector.
 *
 * INPUT:
 *  vec <vector<T>> : mem address of a vector
 *  m <int>         : start index of the subvector
 *  n <int>         : stop index of the subvector
 *
 * OUTPUT:
 *  vector <vector> : an m-n subvector of the vector
 ****************************************************************************************************/
template<typename T>
std::vector<T> subvector(std::vector<T> const &vec, int m, int n) {
   auto first = vec.begin() + m;
   auto last = vec.begin() + n + 1;
   std::vector<T> vector(first, last);
   return vector;
}

/****************************************************************************************************
 * Calculates the minimum scoring threshold for second pass alignment.
 *
 * Min Score calculation:
 *  The minimum score is calculated to attempt to approximate a score of a sequence with an error
 *  frequency of 1/error_rate, so 1 error in every 4 bases where error_rate=4. The obs_match uses
 *  the premise that errors will not occur in a consistent manner, as: (match, match, match, gap)*,
 *  since if k>=error_rate, here will be 0 kmers hitting to the read array and consequently will yield
 *  a score of 0. Therefore, the approximate observed rate of error (mismatch) events is resolved by
 *  assuming that you will see error_rate-1 kmer hits followed by a sequence of error events which will
 *  resolve to 1/error_rate errors over a given sequence.
 *
 *  The obs_window refers to the number of observation windows, within each of which you will observe
 *  obs_match match events and obs_error error events. The min_score is therefore calculated using the
 *  following equation:
 *
 *  H-(wg + n(e(w-1)))
 *
 *  where:
 *   H = the score of a perfect match to the query sequence (repeats collapsed)
 *   w = number of observation windows
 *   g = gap open penalty
 *   n = error/mismatch/gap extension penalty
 *   e = expected error events in an observation window
 *
 * INPUT:
 *  kslen <sizt_t>    : the length of the kmer set generated from the target sequence
 *  er <double>       : the maximum nucleotide error rate across the alignment
 *  k <int>           : the kmer size
 *  hit <double>      : the value to increment the score by for a hit event
 *  gap_open <double> : the value to penalise the score by for a gap open event
 *  neg <double>      : the value to penalise the score by for a gap extension event
 *
 * OUTPUT:
 *  minScoreThreshold <double> : a minimum scoring threshold for second pass screening
 ****************************************************************************************************/
double minscore(  size_t kslen,
                  double er,
                  int k,
                  double hit,
                  double gap_open,
                  double neg )
{
  double obs_match = k+er-1;            // a window which will be a length of kmer hits followed by an arbitrarily long string of mismatches which will resolve to 1/F over the seq
  double obs_error = obs_match/er;      // number of kmers which would need to fail to map to result in 1/error_rate base mismatches within the window
  double obs_windows = kslen/obs_match; // len_pattern = kmer array length as a set to account for repetetive elements
  double max_score = kslen*hit;         // the maximum alignment score assuming 100% match

  double minScoreThreshold = max_score-((obs_windows*gap_open) + (neg*(obs_error*(obs_windows-1))));  // minimum scoring threshold assuming errors aggregate

  return minScoreThreshold;
}

#endif /* BMSPSCREEN_HPP */
