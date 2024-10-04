#ifndef FASTQ_HPP
#define FASTQ_HPP

#include "utilities.hpp"

// class for handling FASTQ data on disk
class FastQ
{
  const std::string _files;

  public:
    std::string _fq1;
    std::string _fq2;
    bool _paired;
    std::vector<std::string> _parts;
    void partition( int );

    FastQ( std::string fastq )
    {
      /* FastQ class constructor
      */
      try {
        std::size_t splitpos = fastq.find(",");

        if ( splitpos != std::string::npos ) {
          _paired = true;
          _fq1 = fastq.substr( 0, splitpos );
          _fq2 = fastq.substr( splitpos+1, fastq.size() );
        } else {
          _paired = false;
          _fq1 = fastq;
        }
      }
      catch ( std::string fastq ) {
        vprint( "FQ-ERROR", "Cannot parse FQ files: " + fastq, "r" );
        exit(1);
      }
    }

    ~FastQ(void) { /* Destructor */ }

  private:
    void runPart( std::string, int, std::string );
    std::int32_t countReads( std::ifstream& inFQ );
};


/****************************************************************************************************
 * Partition the fastq read set.
 *
 * INPUT:
 *  parts <int> : the number of partitions to split the fastq read set into
 *
 * OUTPUT:
 *  read_hits <vector<string>> : a vector of reads which contain the target sequence
 ****************************************************************************************************/
void FastQ::partition( int parts )
{
  if ( _paired ) {
    if ( parts > 2 ) {
      // check to see if partitioning is necessary
      vprint( "PREPROCESSING", "Partitioning paired reads for multithreading...", "g" );
      runPart( _fq1, parts/2, "_1" );
      runPart( _fq2, parts/2, "_2" );
    } else {
      // otherwise just populate _parts with the paired read files
      _parts.push_back( _fq1 );
      _parts.push_back( _fq2 );
    }
  } else {
    if ( parts > 1 ) {
      // check to see if partitioning is necessary
      vprint( "PREPROCESSING", "Partitioning unpaired reads for multithreading...", "g" );
      runPart( _fq1, parts, "_UNP" );
    } else {
      // otherwise just populate _parts with the single unpaired read file
      _parts.push_back( _fq1 );
    }
  }
}

void FastQ::runPart( std::string file, int p, std::string fq_ID )
{
  /*  Partitions a file into p partitions by reads entry (4 lines per read)
  */

  // open the file in binary, check, and store in string buffer
  std::ifstream inFQ(file, std::ios::binary);

  try {
    if ( !inFQ.good() ) {
      throw file;
    }
  }
  catch ( std::string fa_file ) {
    vprint("FQ-ERROR", "Cannot open FASTQ file : " + file, "r" );
    exit(1);
  }

  // count the number of reads in the fastq file
  std::size_t readcount = countReads( inFQ );
  // cout << readcount << endl;
  
  // Reset the file pointer to the beginning
  inFQ.clear();
  inFQ.seekg(0, std::ios::beg);

  // Calculate and store the number of reads in each partition
  std::vector<std::size_t> partvec;
  std::size_t base_partsize = (int)readcount / p;
  int mod_p = readcount % p;
  int sum_p = 0;

  for ( int i=0; i<p; i++ ) {
    if ( i < mod_p ) {
      partvec.push_back(base_partsize+1);
    } else {
      partvec.push_back(base_partsize);
    }
  }

  std::vector<std::string> fqouts_vec;

  // populate the _parts vector with partition files
  for ( int i=0; i<p; i++ ) {
    std::string filepath = "./tmp.p" + std::to_string(i) + fq_ID + ".fq";
    fqouts_vec.push_back( filepath );
    _parts.push_back( filepath );
  }

  std::string line;
  std::int32_t partline = 0;
  int iptr = 0;
  partline += partvec[iptr]*4;

  // cout << _parts[iptr] << " : " << partline/4 << endl;
  std::ofstream fqout;
  fqout.open(fqouts_vec[iptr]);

  // Iterate through the FASTQ file line by line
  std::int32_t lc = 0;
  while ( std::getline(inFQ, line) ) {
    lc++; // Increment line counter

    // Write the line to the current partition file
    fqout << line << "\n";

    // Check if the end of the partition is approaching
    if ( lc == partline && iptr < p ) {
      // Close the current partition file
      fqout.close();

      // Increment the partition pointer
      iptr++;

      // Update the line counter for the next partition
      partline += partvec[iptr]*4;

      // Open the next partition file
      fqout.open(fqouts_vec[iptr]);
    }
  }

  // Close the last partition file
  fqout.close();

  // Close the input FASTQ file
  inFQ.close();
}

std::int32_t FastQ::countReads( std::ifstream& inFQ )
{
  // Count the number of reads in the file
  std::int32_t readcount;
  std::int32_t numlines=0;
  std::string line;
  while ( std::getline(inFQ, line) ) {
    numlines++;
  }

  try {
    if ( numlines % 4 != 0 ) {
      throw numlines;
    } else {
      readcount = numlines/4;    // calc the number of reads in the file
    }
  }
  catch ( std::int32_t numlines ) {
    vprint("PART-ERROR", "FASTQ malformed. linecount % 4 non 0. Please check that there are no malformed reads. Exiting...", "r");
    // cout << "Line count : " << linecount << endl;
    exit (1);
  }

  return readcount;
}

#endif /* FASTQ_HPP */
