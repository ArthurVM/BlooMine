#ifndef FQREAD_HPP
#define FQREAD_HPP

#include "SPscreenutils.hpp"


class FQread
{
  public:
    std::string _id;
    std::string _readseq;
    int _k;
    size_t _readlen;
    double _h, _go, _ge;

    FQread( const std::string& seq, const int k, const double hit, const double gap_open, const double gap_extend, const std::string& id )
      : _id(id), _readseq(seq), _k(k), _readlen(seq.size()), _h(hit), _go(gap_open), _ge(gap_extend) {}
    ~FQread() { /* FQread class destructor */ }

    bool FPscreen( const BloomFilter&, const int );
    bool SPscreen( const std::unordered_set<std::string>&, double, SubAln* out_aln = nullptr );


  private:
    std::unordered_set<std::string> _readkset;

    // generates a kmer set from a sequence
    std::unordered_set<std::string> genKmerSet( std::string, const int );

    // generates a map of kmers and the positions they map to in the read
    std::map<std::string, std::vector<int>> genKmerPosMap( std::string, const int );

    // performs a kmer alignment by mapping from a target sequence to a read
    SubAln kmerAlign( const std::unordered_set<std::string>& );

    // splits the kmer alignment into subalignments for maximum subalignment scoring
    std::vector<SubAln> splitSubalignments( SubAln);

    // finds the maximum scoring subalignment
    SubAln findMaxSubalignment( std::vector<SubAln> );

    // concatenate a subalingment vector
    SubAln concatSubalnChunk( std::vector<SubAln> );

    // takes a subalignment chunk and returns the alignment score of that chunk
    int scoreSubAlignment( std::string );

    // removes trailing gap ragions from an alignment chunk
    SubAln removeTrailing( std::string );
};


/****************************************************************************************************
 * FIRST-PASS SCREEN
 *
 * Perform first-pass screening of _readseq using the given BloomFilter and threshold.
 *
 * INPUT:
 *  BF <BloomFilter> : an instance of the BloomFilter class, generated using a target sequence
 *  threshold <int>  : a threshold of kmer hits to call a read hit
 *
 * OUTPUT:
 *  <bool> : true if the read has similar kmer contents to the target sequence, else false
 *
 ****************************************************************************************************/
bool FQread::FPscreen( const BloomFilter& BF, const int threshold  )
{
  // Fast early-exit counter; only unique kmers count toward the threshold.
  if (threshold <= 0) { return true; }

  int hit_count = 0;
  const int limit = static_cast<int>(_readseq.size()) - _k + 1;
  if (limit <= 0) { return false; }

  std::unordered_set<std::string> seen;
  seen.reserve(static_cast<size_t>(limit));

  for ( int i = 0; i < limit; i++ ) {
    std::string kmer = _readseq.substr(i, _k);

    // skip duplicate kmers
    if ( !seen.insert(kmer).second ) { continue; }

    if ( BF.check(kmer) && ++hit_count >= threshold ) {
      return true;  // early exit when threshold reached
    }
  }

  return false;
}


/****************************************************************************************************
 * Takes a sequence and a value of k and generates a kmer set.
 * N.B. this is different to a kmer array, in that repeated elements are collapsed and no counter is
 * provided.
 *
 * INPUT:
 *  seq <string> : DNA sequence to digest into kmers
 *  k <int>      : kmer size
 *
 * OUTPUT:
 *  kmer_set <unordered_set<string>> : an unordered set of kmers contained within the given sequence
 ****************************************************************************************************/
std::unordered_set<std::string> FQread::genKmerSet( std::string seq, const int k )
{
  std::unordered_set<std::string> kmer_set;

  for ( int i = 0; i < seq.size()-k+1; i++ ) {
    std::string kmer = seq.substr(i, k);

    kmer_set.insert(kmer);
  }
  return kmer_set;
}


/****************************************************************************************************
 * SECOND-PASS SCREEN
 *
 * Perform second-pass scoring of the target sequence. Calls alignment of the target against the
 * read and returns true if the score of the maximum scoring subalignment exceeds the minimum score
 * threshold.
 *
 * DEVELOPERS NOTE: THIS ALGORITHM IS DESIGNED TO BE USED IN BLOOMINE ALONE AS A SECOND-PASS FILTER.
 * BECAUSE OF ITS HIGH TOLERANCE TO GAPS AND MISMATCHES, WORKING ON THE ASSUMPTION THAT THE MAJORITY
 * OF FALSE POSITIVE READS WILL BE DISMISSED DURING THE FIRST PASS SCREEN, IT HAS A LOW SPECIFICITY
 * AND THEREFORE PRONE TO REPORTING LARGE NUMBERS OF FALSE POSITIVE READS IF USED ON AN UNFILTERED
 * .FASTQ FILE. IT HAS BEEN DESIGNED WITH THE INTENTION OF ALLOWING THE USER TO TUNE THE FIRST-PASS
 * FILTER TO BE FAR MORE LIBERAL WITH WHAT IT REPORTS, TO REDUCE LOSS OF INFORMATION, AND TO PERFORM
 * THIS MORE ALIGNMENT DRIVEN APPROACH TO IDENTIFY SEQUENCE HOMOLOGY IN A TIME AND MEMORY EFFICIENT
 * MANNER. WE HAVE FOUND THIS TO BE A FAR MORE EFFECTIVE APPROACH THAN FINE TUNING THE FIRST PASS
 * FILTER TO ATTEMPT TO DO ALL  THE WORK, WHICH IS BOTH TIME CONSUMING AND CAN EASILY LEAD TO
 * INFORMATION LOSS.
 *
 * Given a reference kmer array (set), a read, and a minimum scoring threshold, it will generate a
 * positional kmer map for the read, capturing the relative positions of each kmer and identify
 * sequences of kmers within the read array that intersect with those in the reference array.
 *
 * Alignment is inferred by identifying sequences of kmer hits (referred to as subalignments),
 * and scoring them in a manner similar to a standard pairwise aligner, where scores are improved by
 * the presence of longer, or multiple subalignments, and penalised by gaps.
 *
 * There are a number of features specific to this function which make it differ from a conventional
 * 'mapping by alignment' tool such as bowtie. Notably it is that it is far less strict than the 2
 * mismatches in every 10 bases rule that is emplyed by the aforementioned.
 *
 * This algorithm can be represented as follows:
 *
 *  +--> | parse next read & kmerise |
 *  |               |
 *  |               v
 *  |   | perform positional intersection |
 *  |   | between reference and read kmer |
 *  |   | arrays                          |
 *  |               |
 *  |               v
 *  |   | search for kontiguities within |
 *  |   | the kmer hits array            |
 *  |               |
 *  |               v
 *  |   | gap > gap threshold |  yes  | split kontiguity   |
 *  |   | opened?             |  ---> | into subalignments |
 *  |               ^                           |
 *  |               |                           |
 *  |   | continue scoring |                    |
 *  |               ^                           |
 *  |               | no                        |
 *  |       | End of array? | <-----------------+
 *  |               | yes
 *  |               v
 *  |   | Score sequential subalignments   |
 *  |   | and return max scoring alignment |
 *  |               |
 *  |   no          v
 *  +------- | end of file? |
 *                 | yes
 *                 v
 *     | write all reads with score |
 *     | above homology threshold   |
 *     | to a .fastq                |
 *
 * INPUT:
 *  seq <string> : a target sequence to screen the reads against
 *  mst <double> : a minimum scoring threshold to score a read hit
 *
 * OUTPUT:
 *  <bool> : true if the read contains the target sequence, else false
 ****************************************************************************************************/
bool FQread::SPscreen( const std::unordered_set<std::string>& target_kset, double mst, SubAln* out_aln )
{
  SubAln max_aln = kmerAlign( target_kset );
  if ( out_aln ) {
    *out_aln = max_aln;
  }
  return ( max_aln.score >= mst ? true : false );
}


/****************************************************************************************************
 * In essence, maps kmers to the read and finds the largest scoring aggragation of mapped kmers.
 *
 * In reality, it initialises a blank array using the read length and fills it with kmers based on
 * where there are intersections between the target and read sequences.
 *
 * INPUT:
 *  target_kset <unordered_set<string>> : an unordered kmer set
 *
 * OUTPUT:
 *  maxmum_scoring_subaln <SubAln> : the maximum scoring subalignment
 *
 ****************************************************************************************************/
SubAln FQread::kmerAlign( const std::unordered_set<std::string>& target_kset )
{
  std::vector<std::string> ctarget_kset;

  for (const auto & kmer : target_kset ) {
    std::string cseq = reverseCompliment( kmer );
    ctarget_kset.push_back(cseq);
  }

  std::unordered_set<std::string> unique_kmer_hits;      // container to store kmers which intersect the read and target arrays

  // generate a kmer posmap to store kmers and positions where they belong in
  std::map<std::string, std::vector<int>> read_posmap = genKmerPosMap( _readseq, _k );

  // generate a zero array of the read to contain intersecting kmers
  // char zero_array[_readseq.size()];
  std::string zero_array;
  for ( int i=0; i<_readseq.size(); i++ ) { zero_array += '-'; }

  // iterate through the target kmer set and adjust the zero array according to where regions intersect within the read
  for (const auto & kmer : target_kset ) {
    std::vector<int> kmer_posvec = read_posmap[kmer];

    for (const auto & pos : kmer_posvec ) {
      for ( int i=0; i<kmer.size(); i++ ) {
        zero_array[pos+i] = kmer[i];
      }
    }
  }

  cout << _id << "\n";
  cout << _readseq << "\n";

  // Check that the zero_array has been filled, if not returns NULL.
  // This is necessary if the false-positive for the Bloom filter screening step is too high, leading
  // to a high number of false positive kmer hits which will not be found during direct comparison in
  // SP screening.
  // These cases can result in no kmers being places into the zero array, resulting in a seg fault
  // when trailing gap regions are being stripped.
  if ( std::count(zero_array.begin(), zero_array.end(), '-') == zero_array.size() ) {
    cout << "CIA brainwashed read at : " << zero_array << endl;

    // initiliase a blank subaln read
    SubAln ciaread{"", 0, 0, 0, 0};

    return ciaread;
  }

  // TODO : adjust these functions to create SubAln structs so that I can capture the start and end positions of each subalingment in the read
  // Done
  SubAln stripped_subaln = removeTrailing( zero_array );

  // split the alignment and find the max scoring subalignment
  std::vector<SubAln> chunk_vec = splitSubalignments( stripped_subaln );

  // TODO: this function is a bit flimsy and causes an error without the if catch
  SubAln max_aln = findMaxSubalignment( chunk_vec );

  std::string leading_blank = std::string(max_aln.start, '-');
  std::string tailing_blank = std::string(_readlen-max_aln.end, '-');

  cout << leading_blank << max_aln.aln << tailing_blank << " SCORE :  " << max_aln.score << "\n\n";

  return max_aln;
}


/****************************************************************************************************
 * Splits an alignment chunk into subalignments according to a gap threshold.
 *
 * gap_threshold = ceil(((hit*k)-go)/neg)
 *
 * where;
 *  hit : hit score
 *  k   : kmer size
 *  go  : gap open penalty
 *  neg : gap extend penalty
 *  ceil: ceiling function
 *
 * This gap threshold is calibrated such that a single kmer hit followed a gap equal to the gap
 * threshold will result in a score of approximately 0. The alignment is broken at each position
 * where a gap exists of >= gap threshold by the splitSubalignments function. Every permutation
 * of sequential subalignments are then scored, totalling n(n/2+0.5) where n is the number of
 * subalignments.
 *
 * E.g. with an alignment of ACT---GACT----ACTG-----GACTGA and a gap threshold of 4, 3 subalignments
 * will be yielded:

 *  ACT---GACT, ACTG & GACTGA.
 *
 * These will form 6 sequential subalignments:

 *  ACT---GACT
 *  ACTG
 *  GACTGA
 *  ACT---GACT----ACTG
 *  ACTG-----GACTGA
 *  ACT---GACT----ACTG-----GACTGA
 *
 * Each of these subalignments will be scored and the maximum scoring alignment returned.
 *
 * INPUT:
 *  aln_chunk <string> : a trailing stripped kmer alignment chunk
 *
 * OUTPUT:
 *  chunk_vec <vector<SubAln>> : a vector containing SubAln structures
 ****************************************************************************************************/
std::vector<SubAln> FQread::splitSubalignments( SubAln subaln_chunk )
{
  // TODO : I think the substr out of range bug is here...
  std::vector<SubAln> chunk_vec;

  int gap_threshold = ceil(((_h*_k)-_go)/_ge);

  std::string aln_fragment;   // declare an align fragment to store the current subalignment chunk for processing

  int c0=0;           // initialise where the start of the subalignment chunk is
  int cn=0;           // initialise where the end of the subalignment chunk is
  int gap_count=0;    // initialise a gap counter to keep track of gap size

  std::string aln_chunk = subaln_chunk.aln;

  for ( int i=0; i<aln_chunk.size(); i++ ) {
    if ( aln_chunk[i]  == '-' ) {
      gap_count++;    // if a - character is encountered, increment the gap counter
    } else {
      gap_count = 0;  // otherwise reset gap counter
    }

    if ( aln_chunk[i] == '-' && gap_count >= gap_threshold && aln_chunk[i+1] != '-' ) {
      cn = i+1-(c0+gap_count);                  // the end of the subaln chunk
      aln_fragment = aln_chunk.substr(c0, cn);  // define the new subaln chunk

      SubAln AlnChunk{ aln_fragment, gap_count, scoreSubAlignment( aln_fragment ), subaln_chunk.start+c0, subaln_chunk.start+c0+cn };   // define a new aln chunk
      // cout << "SPLIT-SUBALN :: " << AlnChunk.aln << " " << AlnChunk.start << ":" << AlnChunk.end << " " << c0 << ":" << cn <<  endl;

      chunk_vec.push_back(AlnChunk);

      c0 = i+1;
      gap_count = 0;
    }
  }

  // capture the last subalignment chunk
  cn = aln_chunk.size()-c0;
  aln_fragment = aln_chunk.substr(c0, cn);

  SubAln AlnChunk{ aln_fragment, 0, scoreSubAlignment( aln_fragment ), subaln_chunk.start+c0, subaln_chunk.start+c0+cn };
  // cout << "SPLIT-SUBALN :: " << AlnChunk.aln << " " << AlnChunk.start << ":" << AlnChunk.end << " " << c0 << ":" << cn <<  endl;

  chunk_vec.push_back(AlnChunk);

  return chunk_vec;
}


/****************************************************************************************************
 * Takes a vector of subalignment structures and calculates the maximum scoring subalignment.
 *
 * This works by calculating the score of concatenated subalignments and finding the maximum
 * scoring chunk.
 *
 * INPUT:
 *  subaln_vec <vector<SubAln>> : a vector of subalignment structures generated from kmer mapping
 *
 * OUTPUT:
 *  none
 ****************************************************************************************************/
SubAln FQread::findMaxSubalignment( std::vector<SubAln> chunk_vec )
{
  SubAln max_scoring_subaln{"", 0, 0, 0, 0};  // initialises a blank subalignment
  SubAln concat_aln{"", 0, 0, 0, 0}; // initialises a concat subalignment

  size_t vecsize = chunk_vec.size();

  // Handle the case of an empty vector
  if (chunk_vec.empty()) {
    return max_scoring_subaln;
  }
  
  max_scoring_subaln = chunk_vec[0];  // reinitialse the first subalignment as max score

  // Handle the case of a single subalignment
  if (chunk_vec.size() == 1) {
    return max_scoring_subaln;
  }


  // Iterates through the subalignment twice to generate all possible concatenations, and scores them
  // This is O(horrible) but works fine in practice because NGS reads are small...
  for ( int i=0; i<vecsize; i++ ) {
    for ( int j=i; j<vecsize; j++ ) {
      std::vector<SubAln> subvec = subvector( chunk_vec, i, j );
      if ( subvec.size() > 1 ) {
        concat_aln = concatSubalnChunk( subvec );
      } else {
        concat_aln = subvec[0];
      }
      
      if ( concat_aln.score >  max_scoring_subaln.score ) {
        max_scoring_subaln = concat_aln;
      }
    }
  }

  return max_scoring_subaln;
}


/****************************************************************************************************
 * Takes an alignment vector and returns a SubAln struct of the constituent SubAlns concatenated.
 *
 * INPUT:
 *  subaln_vec <vector<SubAln>> : a vector of subalignment structures generated from kmer mapping
 *
 * OUTPUT:
 *  concat_sa <SubAln> : a subalignment structure formed from the concatenation of all SubAln's in
 *                       the input vector
 ****************************************************************************************************/
SubAln FQread::concatSubalnChunk( std::vector<SubAln> subaln_vec )
{
  std::string seq;
  int score = 0;
  size_t vecsize = subaln_vec.size();

  int concat_start = subaln_vec[0].start;
  int concat_end = subaln_vec[vecsize-1].end;

  for ( int i=0; i<vecsize; i++ ) {
    seq+=subaln_vec[i].aln;
    score+=subaln_vec[i].score;

    // cout << "CAT-SUBALN :: "<< subaln_vec[i].aln << " " << subaln_vec[i].start << ":" << subaln_vec[i].end << endl;

    if ( i<subaln_vec.size()-1 ) {
      seq+=std::string(subaln_vec[i].gap, '-');
      score-=_go;
      score-=_ge*subaln_vec[i].gap-1;
    }
  }
  int gap = subaln_vec[vecsize-1].gap;
  SubAln concat_sa{ seq, gap, score, concat_start, concat_end };

  return concat_sa;
}


/****************************************************************************************************
 * Takes a subalignment chunk and calculates the score.
 *
 * The scoring system follows an affine convention, whereby gap open events are penalised more
 * heavily than gap extension events.
 *
 * INPUT:
 *  aln_chunk <string> : the subalignment chunk generated by the splitSubalignments function
 *
 * OUPUT:
 *  score <int> : the score of the individual subalignment
 ****************************************************************************************************/
int FQread::scoreSubAlignment( std::string aln_chunk )
{
  int score = 0;
  int gap = 0;
  for ( auto & pos : aln_chunk ) {
    if ( pos != '-' ) {
      score+= _h;
      gap=0;
    } else if ( pos == '-' && gap == 0 ) {
      score-=_go;
      gap+=1;
    } else if ( pos == '-' && gap > 0 ) {
      score-=_ge;
      gap+=1;
    }
  }
  return score;
}


/****************************************************************************************************
 * Takes a sequence and a value of k and generates a kmer array map of kmer:map_position pairs.
 *
 * INPUT:
 *  seq <string> : DNA sequence to digest into kmers
 *  k <int>      : kmer size
 *
 * OUTPUT:
 *  kmer_posmap <map<string,int>> : map containing kmers and positions they map to in the read
 ****************************************************************************************************/
std::map<std::string, std::vector<int>> FQread::genKmerPosMap( std::string seq, const int k )
{
  std::map<std::string, std::vector<int>> kmer_posmap;

  for ( int i = 0; i < seq.size()-k+1; i++ ) {
    std::string kmer = seq.substr(i, k);

    kmer_posmap[kmer].push_back(i);
  }
  return kmer_posmap;
}


/****************************************************************************************************
 * Remove trailing non-aligned regions from a raw alingment chunk.
 *
 * INPUT:
 *  raw_chunk <string> : a raw alignment chunk with trailing gap characters
 *
 * OUTPUT:
 *  stripped_chunk <string> : the alignment chunk with trailing gap characters stripped
 ****************************************************************************************************/
SubAln FQread::removeTrailing( std::string raw_chunk )
{
  int chunk_start = 0;
  int chunk_end = 0;

  for ( int fi=0; fi<_readseq.size(); fi++ ) {
    if ( raw_chunk[fi] != '-' ) {
      chunk_start = fi;
      break;
    }
  }

  for ( int bi=_readseq.size()-1; bi>chunk_start; bi-- ) {
    if ( raw_chunk[bi] != '-' ) {
      chunk_end = bi;
      break;
    }
  }

  // cout << chunk_start << " " << chunk_end << " " << chunk_end-chunk_start+1 << endl;
  std::string stripped_chunk = raw_chunk.substr(chunk_start, chunk_end-chunk_start+1);
  // std::string stripped_chunk;
  // for ( int i=chunk_start, j=0; i<chunk_end; i++, j++ ) {
  //   stripped_chunk += raw_chunk[i];
  // }

  SubAln stripped_subaln{ stripped_chunk, chunk_start, 0, chunk_start, chunk_end };

  return stripped_subaln;
}

#endif /* FQREAD_HPP */
