#ifndef BLOOMFILTER_HPP
#define BLOOMFILTER_HPP

#include <vector>
#include <math.h>

class BloomFilter
{
  unsigned int _elementsN;
  float _fp;

  public:
    unsigned int _baN, _hashN;
    std::vector<bool> _bitarray;

    // BloomFilter constructor and destructor
    BloomFilter( int, float );
    ~BloomFilter();

    // calculate the size of the Bloom Filter using the number of elements to be added and the desired FP rate
    int getsize( int, float );

    // calculate the number of hash functions to use using the FP rate and the size of the Bloom Filter
    int gethashN( int, int );

    // expose an element to the Bloom Filter for storage
    void push( std::string );

    // check an element against the Bloom Filter to infer membership
    bool check( std::string );
};

BloomFilter::BloomFilter( int els, float fp )
{
  /* BloomFilter constructor.
  where:
    els (int) : the number of elements to be stored in the Bloom Filter. In the case of BlooMine this is the size of the kmer set of the target sequence
    fp (float) : the desired false positive rate for the Bloom Filter
  */
  _elementsN = els;
  _fp = fp;
  _baN = getsize( _elementsN, _fp );
  _hashN = gethashN( _baN, _elementsN );

  // Please dont judge me for this, I spent 5 hours trying to find a way of initialising an empty bit array of length N without iteration, but couldnt get it to work on a previously defined bool vector...
  // TODO : fix this horrible way of doing things. It's a good thing BF bit arrays are supposed to be short.
  for (int i=0; i<=_baN; i++) {
    _bitarray.push_back(false);
  }
}

BloomFilter::~BloomFilter(void){ /* BloomFilter class destructon */ }

int BloomFilter::getsize( int _elementsN, float _fp )
{
  /* Calculates the size of the bit array as:
    m = -(n * log(p) / log(2)^2)
  where:
    n _elementsN (int) : the number of elements to be stored in the Bloom BloomFilter
    p _fp (float) : the desired false positive rate

  returns m
  */
  int m = -(_elementsN * log(_fp) / pow(log(2), 2));
  return m;
}

int BloomFilter::gethashN( int _baN, int _elementsN )
{
  /* Calculates the number of hash functions to use as:
    k = (m/n) * log(2)
  where:
    m _baN (int) : the size of the bit array as calculates by getsize
    n _elementsN (int) : the number of elements to be stored in the Bloom Filter

  returns k
  */
  int k = (_baN/_elementsN) * log(2);
  return k;
}

void BloomFilter::push( std::string element )
{
  /* Expose an element to the hash functions and adds it to the Bloom Filter
  where:
    element (string) : the element to check against the filter
  */
  for (int i=0; i<_hashN; i++) {
    std::hash<std::string> str_hash;                      // declaire the hash
    std::string unique_el = element + std::to_string(i);  // hacks a unique element using the hash number
    int hashval = str_hash(unique_el) % _baN;   // generate the hash value

    _bitarray[hashval] = 1;                  // adjust the index corresponding to the hash value in the Bloom Filter to True (1)
  }
}

bool BloomFilter::check( std::string element )
{
  /* Check an element against the Bloom Filter
  where:
    element (string) : the element to check against the filter

  returns true if the element putiatively exists within the BF, else false
  */
  for (int i=0; i<_hashN; i++) {
    std::hash<std::string> str_hash;                      // declaire the hash
    std::string unique_el = element + std::to_string(i);  // hacks a unique element using the hash number
    int hashval = str_hash(unique_el) % _baN;   // generate the hash value

    if (_bitarray[hashval] == false) {
      return false;
    }
  }
  return true;
}

#endif /* BLOOMFILTER_HPP */
