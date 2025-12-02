/*
Utility functions for the BlooMine module.
*/

#include "BlooMineUtils.hpp"
#include "FastQ.hpp"

#include <mutex>
#include <atomic>

struct ScoreLog {
  std::string read_id;
  int flank;
  int score;
  double threshold;
  bool reverse_complement;
  bool pass;
};

static std::mutex score_log_mutex;
static std::vector<ScoreLog> score_logs;
static std::atomic<int> g_flank_number{0};

void setFlankNumber(int flank) {
  g_flank_number.store(flank);
}

static void appendScoreLog(const std::string& read_id,
                           int score,
                           double threshold,
                           bool rc,
                           bool pass) {
  std::lock_guard<std::mutex> lock(score_log_mutex);
  score_logs.push_back(ScoreLog{read_id, g_flank_number.load(), score, threshold, rc, pass});
}

static void clearScoreLog() {
  std::lock_guard<std::mutex> lock(score_log_mutex);
  score_logs.clear();
}

static void writeScoreLog(const std::string& prefix) {
  std::lock_guard<std::mutex> lock(score_log_mutex);
  std::ofstream out(prefix + "_flank_scores.tsv", std::ios::out);
  if (!out.is_open()) {
    std::cerr << "LOG-ERROR: Cannot open flank score log for writing\n";
    return;
  }
  out << "read_id\tflank\tscore\tthreshold\treverse_complement\tpass\n";
  for (const auto& rec : score_logs) {
    out << rec.read_id << "\t"
        << rec.flank << "\t"
        << rec.score << "\t"
        << rec.threshold << "\t"
        << (rec.reverse_complement ? "1" : "0") << "\t"
        << (rec.pass ? "1" : "0") << "\n";
  }
  out.close();
}


/****************************************************************************************************
 * Generate a Bloom Filter from the target sequence.
 *
 * INPUT:
 *  vm <variables_map> : a boost program options map containing argument variables
 *  threshold <int&>   : mem addr of the threshold variable
 *
 * OUTPUT:
 *  BF <BloomFilter> : a Bloom Filter generated from the target sequence, using a given FP rate
 *
 * MODIFIES:
 *  threshold <int> main : the threshold for kmer hits to call a read hit
 ****************************************************************************************************/
std::pair<BloomFilter, int> generateBloomFilter( std::string target_fasta, 
                                                 int kmer, 
                                                 float first_pass_similarity, 
                                                 float fp_rate ) {

  std::unordered_set<std::string> kmer_array;

  try {
    kmer_array = fasta2kmerset( target_fasta, kmer );
  }
  catch ( const std::out_of_range& oor ) {
    vprint( "FA-ERROR", "Generation of kmer array from " + target_fasta + " failed. Exiting...", "r" );
    exit(1);
  }

  int threshold = kmer_array.size() * first_pass_similarity/100;     // define the minimum number of kmers that must match to return a read hit
  // cout << "threshold : " << threshold << endl;         // uncomment for debugging

  BloomFilter BF( kmer_array.size(), fp_rate );

  for ( auto & kmer : kmer_array ) { BF.push(kmer); } // iterates through the target kmer array and pushes each kmer into the Bloom Filters

  return std::make_pair(BF, threshold);;
}


/****************************************************************************************************
 * Calculate the Minimum Scoring Threshold (MST) for this run.
 *
 * INPUT:
 *  vm <variables_map>                    : a boost program options map containing argument variables
 *  target_kset <unordered_set<string>>   : a set of kmers generated from the target sequence
 *
 * OUTPUT:
 *  spMST <double> : the Minimum Scoring Threshold
 ****************************************************************************************************/
double calculateMinimumScoreThreshold( float second_pass_error,
                                       int kmer, 
                                       std::unordered_set<std::string> target_kset ) {

  const double spMST =  minscore( target_kset.size(), second_pass_error, kmer, hit, gap_open, neg );
  vprint("SP-SCREEN", "MINIMUM SCORE THRESHOLD SET TO : " + std::to_string(spMST), "g");

  return spMST;
}


/****************************************************************************************************
 * Partition the reads within the input FASTQ for multithreading
 *
 * INPUT:
 *  fastq <string>  : the input fastq file to partition
 *  threads <int>   : the number of threads for partitioning
 *
 * OUTPUT:
 *  read_data_parts <vector<vector<string>>> : a vector of read partitions
 ****************************************************************************************************/
std::vector<std::vector<std::string>> partitionData( std::string fastq, 
                                                     int threads ) {

  vprint( "PREPROCESSING", "Partitioning paired reads for multithreading...", "g" );

  // Read FASTQ file
  std::vector<std::string> read_data = readFQ( fastq );

  size_t chunk_size = (read_data.size()/4) / threads;

  std::vector<std::vector<std::string>> read_data_parts;

  for (int tnum = 0; tnum < threads; ++tnum) {
    size_t chunk_lines = chunk_size*4;
    size_t start = tnum * chunk_lines;
    size_t end = (tnum == threads - 1) ? read_data.size() : (tnum + 1) * chunk_lines;
    std::vector<std::string> chunk(read_data.begin() + start, read_data.begin() + end);

    read_data_parts.push_back(chunk);
  }

  return read_data_parts;
}


/****************************************************************************************************
 * Partition the reads within the input FASTQ for multithreading
 *
 * INPUT:
 *  fastq <string>  : the input fastq file to partition
 *  threads <int>   : the number of threads for partitioning
 *
 * OUTPUT:
 *  FQ <FastQ> : a FastQ object containing read partitions
 ****************************************************************************************************/
FastQ partitionData( std::string fastq, 
                     int threads,
                     bool on_disk ) {

  // Read FASTQ file
  FastQ FQ( fastq );
  FQ.partition( threads );

  return FQ;
}


/****************************************************************************************************
 * Checks and reads a FASTQ file.
 *
 * INPUT:
 *  fq_file <string> : the path to a FASTQ file
 *
 * OUTPUT:
 *  read_data <vector> : a vector containing all reads within the FASTQ file
 ****************************************************************************************************/
std::vector<std::string> readFQ( std::string fastq_in ) {

  std::vector<std::string> read_data;
  std::stringstream ss(fastq_in);
  std::string fq_file;

  while (std::getline(ss, fq_file, ',')) {
    std::ifstream inFQ(fq_file, std::ios::binary);

    if (!inFQ.good()) {
      std::cerr << "FQ-ERROR: Cannot open FASTQ file: " << fq_file << std::endl;
      exit(1);
    }

    std::string line;
    while (std::getline(inFQ, line)) {
      read_data.push_back(line);
    }
    inFQ.close();
  }

  return read_data;
}


/****************************************************************************************************
 * Spawn threads to run the BlooMine pipeline on each read partition in memory.
 *
 * INPUT:
 *  vm <po:variables_map>                       : a boost program options map containing argument 
 *                                                variables
 *  read_data_parts <vector<vector<string>>>    : a vector of read partitions
 *  BF <BloomFilter>                            : a Bloom filter generated using the target sequence
 *  fp_threshold <int>                          : the threshold for kmer hits to call a read hit 
 *                                                during first-pass screening
 *
 * OUTPUT:
 *  read_hits <vector<string>> : a vector of reads which contain the target sequence
 ****************************************************************************************************/
std::vector<std::future<std::vector<std::string>>> spawnThreads( int threads,
                                                                 int kmer,
                                                                 std::string prefix,
                                                                 std::vector<std::vector<std::string>> read_data_parts, 
                                                                 BloomFilter BF,
                                                                 std::unordered_set<std::string> target_kset,
                                                                 int fp_threshold,
                                                                double spMST ) {
  
  vprint( "FP-SCREEN", "Starting BlooMine screening pipeline...", "g" );

  clearScoreLog();

  int thread_num = threads;

  // ThreadPool pool(thread_num);

  // Create a vector of futures to store the results from each thread
  std::vector<std::future<std::vector<std::string>>> futures;

  // for ( int i = 0; i < thread_num; i++ ) {
  //   cout << "\t\tThread created : " << i << "\tfqpart : " << read_data_parts[i].size() << endl;
  //   if ( i == thread_num - 1 ) { cout << endl; }
  //   futures.push_back( pool.enqueue( runBM, vm, read_data_parts[i], BF, target_kset, fp_threshold, spMST ) );
  // }
  
  // spawn threads
  for (const auto& read_data_part : read_data_parts ) {
    std::future<std::vector<std::string>> fp = std::async( runBM, kmer, read_data_part, BF, target_kset, fp_threshold, spMST );
    futures.push_back(std::move(fp));
  }
  
  // initialise and open read hit outfile
  std::ofstream fp_outfile;
  fp_outfile.open( prefix + "_BMfiltered.fq", std::ios::out );

  // Collect the results from the futures
  for (auto& future : futures) {
    std::vector<std::string> hitsvec = future.get();
    for (auto& read : hitsvec) {
      fp_outfile << read;
    }
  }
  
  fp_outfile.close();

  writeScoreLog(prefix);

  return futures;
}


/****************************************************************************************************
 * PARALLEL FUNCTION - RUN BLOOMINE PIPELINE FROM MEMORY.
 * Reads through and screens each read within the fq_part fastq file.
 *
 * INPUT:
 *  kmer <int>                          : kmer size
 *  read_data_part <vector<string>>     : vector of partitioned reads
 *  BF <BloomFilter>                    : a Bloom filter generated using the target sequence
 *  target_kset <unordered_set<string>> : a set of kmers generated from the target sequence
 *  fp_threshold <int>                  : a threshold for number of kmer hits to call a read hit during
 *                                        first-pass screening
 *  spMST <double>                      : the Minimum Scoring Threshold
 * 
 * OUTPUT:
 *  read_hits <vector<string>> : a vector of reads which contain the target sequence
 ****************************************************************************************************/
std::vector<std::string> runBM( int kmer,
                                std::vector<std::string> read_data_part,
                                BloomFilter BF,
                                std::unordered_set<std::string> target_kset,
                                int fp_threshold,
                                double spMST ) {

  // Iterate through the FASTQ file line by line
  std::vector<std::string> read_hits;
  std::vector<std::string> read_box;
  std::string line;
  std::int32_t lc = 0;

  for (size_t i = 0; i < read_data_part.size(); i++) {
    lc++; // Increment line counter

    line = read_data_part[i];

    if ( lc % 4 != 0 ) {  // check if the linecounter has reached a new block
      read_box.push_back(line);
    } else {
      read_box.push_back(line);    // capture the description line
      
      // Check if the read is complete
      if (read_box.size() == 4) {

        // Screen block
        const std::string read_id = read_box[0];
        FQread Read( read_box[1], kmer, hit, gap_open, neg, read_id );

        // Screen reads
        if ( Read.FPscreen( BF, fp_threshold ) ) {            // if read passes first-pass screen
          SubAln aln{"",0,0,0,0};
          bool sp_pass = Read.SPscreen( target_kset, spMST, &aln );
          appendScoreLog(read_id, aln.score, spMST, false, sp_pass);
          if ( sp_pass ) {        // if read passes second-pass screen

            std::string read_concat;
            for ( const auto & l : read_box ) { read_concat += l + "\n"; }
            read_hits.push_back(read_concat);
            // std::cout << read_concat;
          }
        } else {
          // reverse complement the read sequence and rerun through the BF, adding to the counter on hits
          // Screen block for complementing the read
          std::string readseqRC = reverseCompliment( read_box[1] );
          FQread Read( readseqRC, kmer, hit, gap_open, neg, read_id );

          if ( Read.FPscreen( BF, fp_threshold ) ) {            // if reverse complemented read passes first-pass screen
            SubAln aln{"",0,0,0,0};
            bool sp_pass = Read.SPscreen( target_kset, spMST, &aln );
            appendScoreLog(read_id, aln.score, spMST, true, sp_pass);
            if ( sp_pass ) {        // if reverse complemented read passes second-pass screen
            
              std::string read_concat;
              for ( const auto & l : read_box ) { read_concat += l + "\n"; }
              read_hits.push_back(read_concat);
            }
          }
        }

        read_box.clear(); // redefine the vector as empty
      }
    }
  }

  return read_hits;
}


/****************************************************************************************************
 * Spawn threads to run the BlooMine pipeline on each read partition on disk.
 *
 * INPUT:
 *  vm <po:variables_map>          : a boost program options map containing argument 
 *                                   variables
 *  fq_file_parts <vector<string>> : a vector of partitioned read files
 *  BF <BloomFilter>               : a Bloom filter generated using the target sequence
 *  fp_threshold <int>             : the threshold for kmer hits to call a read hit 
 *                                   during first-pass screening
 *
 * OUTPUT:
 *  read_hits <vector<string>> : a vector of reads which contain the target sequence
 ****************************************************************************************************/
std::vector<std::future<std::vector<std::string>>> spawnThreads( int threads,
                                                                 int kmer,
                                                                 std::string prefix,
                                                                 FastQ FQ, 
                                                                 BloomFilter BF,
                                                                 std::unordered_set<std::string> target_kset,
                                                                 int fp_threshold,
                                                                 double spMST ) {
  
  vprint( "FP-SCREEN", "Starting BlooMine screening pipeline...", "g" );

  clearScoreLog();

  // ThreadPool pool(thread_num);

  // Create a vector of futures to store the results from each thread
  std::vector<std::future<std::vector<std::string>>> futures;

  // for ( int i = 0; i < thread_num; i++ ) {
  //   cout << "\t\tThread created : " << i << "\tfqpart : " << read_data_parts[i].size() << endl;
  //   if ( i == thread_num - 1 ) { cout << endl; }
  //   futures.push_back( pool.enqueue( runBM, vm, read_data_parts[i], BF, target_kset, fp_threshold, spMST ) );
  // }
  
  // spawn threads
  for ( int i=0; i < threads; i++ ) {
    std::future<std::vector<std::string>> fp = std::async( runBMdisk, kmer, FQ._parts[i], BF, target_kset, fp_threshold, spMST );
    futures.push_back(std::move(fp));
  }
  
  // initialise and open read hit outfile
  std::ofstream fp_outfile;
  fp_outfile.open( prefix + "_BMfiltered.fq", std::ios::out );

  // Collect the results from the futures
  for (auto& future : futures) {
    std::vector<std::string> hitsvec = future.get();
    for (auto& read : hitsvec) {
      fp_outfile << read;
    }
  }
  
  fp_outfile.close();

  writeScoreLog(prefix);

  return futures;
}


/****************************************************************************************************
 * PARALLEL FUNCTION - RUN BLOOMINE PIPELINE FROM DISK.
 * Reads through and screens each read within the fq_part fastq file.
 *
 * INPUT:
 *  kmer <int>                          : kmer size
 *  fq_file <string>                    : path to a read file (or a partition of one) in fastq format
 *  BF <BloomFilter>                    : a Bloom filter generated using the target sequence
 *  target_kset <unordered_set<string>> : a set of kmers generated from the target sequence
 *  fp_threshold <int>                  : a threshold for number of kmer hits to call a read hit during
 *                                        first-pass screening
 *  spMST <double>                      : the Minimum Scoring Threshold
 * 
 * OUTPUT:
 *  read_hits <vector<string>> : a vector of reads which contain the target sequence
 ****************************************************************************************************/
std::vector<std::string> runBMdisk( int kmer,
                                    std::string fq_file,
                                    BloomFilter BF,
                                    std::unordered_set<std::string> target_kset,
                                    int fp_threshold,
                                    double spMST ) {

  std::vector<std::string> read_hits;
  std::ifstream infile(fq_file);
  std::vector<std::string> read_box;
  std::string line;
  std::int32_t lc = 0;

  while (std::getline(infile, line)) {
    lc++; // Increment line counter

    if (lc % 4 != 0) {  // check if the linecounter has reached a new block
      read_box.push_back(line);
    } else {
      read_box.push_back(line);    // capture the description line

      // Check if the read is complete
      if (read_box.size() == 4) {

        // Screen block
        const std::string read_id = read_box[0];
        FQread Read(read_box[1], kmer, hit, gap_open, neg, read_id);

        // Screen reads
        if (Read.FPscreen(BF, fp_threshold)) {            // if read passes first-pass screen
          SubAln aln{"",0,0,0,0};
          bool sp_pass = Read.SPscreen(target_kset, spMST, &aln);
          appendScoreLog(read_id, aln.score, spMST, false, sp_pass);
          if (sp_pass) {        // if read passes second-pass screen

            std::string read_concat;
            for (const auto& l : read_box) { read_concat += l + "\n"; }
            read_hits.push_back(read_concat);
            // std::cout << read_concat;
          }
        } else {
          // reverse complement the read sequence and rerun through the BF, adding to the counter on hits
          // Screen block for complementing the read
          std::string readseqRC = reverseCompliment(read_box[1]);
          FQread Read(readseqRC, kmer, hit, gap_open, neg, read_id);

          if (Read.FPscreen(BF, fp_threshold)) {            // if reverse complemented read passes first-pass screen
            SubAln aln{"",0,0,0,0};
            bool sp_pass = Read.SPscreen(target_kset, spMST, &aln);
            appendScoreLog(read_id, aln.score, spMST, true, sp_pass);
            if (sp_pass) {        // if reverse complemented read passes second-pass screen

              std::string read_concat;
              for (const auto& l : read_box) { read_concat += l + "\n"; }
              read_hits.push_back(read_concat);
            }
          }
        }

        read_box.clear(); // redefine the vector as empty
      }
    }
  }

  return read_hits;
}


/****************************************************************************************************
 * Collect results from threads and write to output file.
 *
 * INPUT:
 *  vm <po:variables_map>                               : a boost program options map containing 
 *                                                        argument variables
 *  hits_collection <vector<future<vector<string>>>>    : a vector of futures containing read hits
 *                                                        from each thread
 *
 * OUTPUT:
 *  None
 ****************************************************************************************************/
// void collectResults( po::variables_map vm, std::vector<std::future<std::vector<std::string>>> hits_collection ) {
//   // initialise and open read hit outfile
//   std::ofstream fp_outfile;
//   fp_outfile.open( vm["prefix"].as<std::string>() + "_BMfiltered.fq", std::ios::out );

//   // Check if the file opened successfully
//   if (!fp_outfile.is_open()) {
//     std::cerr << "FQ-ERROR: Cannot open output file: " << vm["prefix"].as<std::string>() + "_BMfiltered.fq" << std::endl;
//     exit(1);
//   }

//   // Collect results from threads
//   for (auto& fp : hits_collection) {
//     std::vector<std::string> hitsvec = fp.get();
//     for (auto& read : hitsvec) {
//       fp_outfile << read;
//     }
//   }
//   fp_outfile.close();
// }
