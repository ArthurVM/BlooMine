#ifndef FASTQ_HPP
#define FASTQ_HPP

#include "utilities.hpp"

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
    void runpart( std::string, int, std::string );
    std::int32_t count_lines( std::ifstream& );
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
      runpart( _fq1, parts/2, "_1" );
      runpart( _fq2, parts/2, "_2" );
    } else {
      // otherwise just populate _parts with the paired read files
      _parts.push_back( _fq1 );
      _parts.push_back( _fq2 );
    }
  } else {
    if ( parts > 1 ) {
      // check to see if partitioning is necessary
      vprint( "PREPROCESSING", "Partitioning unpaired reads for multithreading...", "g" );
      runpart( _fq1, parts, "_UNP" );
    } else {
      // otherwise just populate _parts with the single unpaired read file
      _parts.push_back( _fq1 );
    }
  }
}

void FastQ::runpart( std::string file, int p, std::string fq_ID )
{
  /*  Partitions a file into p partitions by reads entry (4 lines per read)
  */

  std::ifstream fileRL(file);

  try {
    if ( !fileRL.good() ) {
      throw file;
    }
  }
  catch ( std::string fa_file ) {
    vprint("FQ-ERROR", "Cannot open FASTQ file : " + file, "r" );
    exit(1);
  }

  std::size_t linecount = count_lines( fileRL ); // store the number of lines in the fastq file
  std::size_t readcount;                         // store the number of reads in the fastq file

  fileRL.close(); std::ifstream inFQ(file);       // Gross close and reopen file for iterating through

  // cout << "linecount : " << linecount << endl;

  try {
    if ( linecount % 4 != 0 ) {
      throw linecount;
    } else {
      readcount = linecount/4;    // calc the number of reads in the file
    }
  }
  catch ( std::int32_t linecount ) {
    vprint("PART-ERROR", "FASTQ malformed. linecount % 4 non 0. Please check that there are no fragmented reads. Exiting...", "r");
    // cout << "Line count : " << linecount << endl;
    exit (1);
  }

  // cout << "readcount : " << readcount << endl;

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

  for ( std::int32_t lc = 0; lc < linecount && std::getline(inFQ, line); lc++ ) {
    if ( lc < partline-1 ) {                      // check if the end of the partition is approaching
      fqout << line << "\n";                      // write the line to the partition file
    } else if ( lc == partline-1 && iptr < p ) {   // detects end of partition approaching
      // cout << iptr << " " << partline <<  endl;
      fqout << line << "\n";                      // capture the last line into the partition file
      iptr++;                                     // increments the partition ptr to specify the next partition
      partline += partvec[iptr]*4;                // specify the number of lines in this partition

      fqout.close();                              // close the previous partition file

      // cout << _parts[iptr] << " : " << partline/4 << endl;
      if ( iptr < p ) {
        // cout << "opening new file" << endl;
        fqout.open(fqouts_vec[iptr]);              // open the new partition file
      }
    }
  }
  inFQ.close();
}

std::int32_t FastQ::count_lines(std::ifstream& inFile)
{
  // dirty but quick (ish) way of counting the number of lines in an open file
  std::int32_t numlines=0;
  std::string line;

  // int linecount = std::count(std::istreambuf_iterator<char>(inFile),
  //                            std::istreambuf_iterator<char>(), '\n');

  while ( std::getline(inFile, line) ) {
    numlines++;
  }
  return numlines;
}

#endif /* FASTQ_HPP */
