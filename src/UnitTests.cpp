#include "gtest/gtest.h"
#include "BlooMineUtils.hpp"
#include <fstream>
#include <sstream>

// Test fixture for BlooMineUtils functions
class BlooMineTest : public ::testing::Test {
protected:
  // Setup for each test
  void SetUp() override {
    // Create a temporary FASTA file for testing
    std::ofstream temp_fasta("temp.fasta");
    temp_fasta << ">target_sequence\n"
              << "ACGTACGTACGT\n";
    temp_fasta.close();

    // Create a temporary FASTQ file for testing
    std::ofstream temp_fastq("temp.fastq");
    temp_fastq << "@read1\n"
              << "ACGTACGTACGT\n"
              << "+\n"
              << "IIIIIIIIII\n"
              << "@read2\n"
              << "ATCGATCGATCG\n"
              << "+\n"
              << "IIIIIIIIII\n";
    temp_fastq.close();

    // Create a boost program options map for testing
    vm.emplace("target_fasta", "temp.fasta");
    vm.emplace("fastq", "temp.fastq");
    vm.emplace("kmer", 7);
    vm.emplace("false_positive", 0.0001);
    vm.emplace("FP-sim", 50.0);
    vm.emplace("SP-error", 4.0);
  }

  // Teardown for each test
  void TearDown() override {
    // Remove temporary files
    std::remove("temp.fasta");
    std::remove("temp.fastq");
  }

  // Boost program options map for testing
  po::variables_map vm;
};

// Test for generateBloomFilter function
TEST_F(BlooMineTest, GenerateBloomFilter) {
  auto [BF, threshold] = generateBloomFilter(vm);

  // Check if the Bloom filter is generated correctly
  ASSERT_NE(BF.size(), 0);
  ASSERT_EQ(threshold, 3);

  // Check if the Bloom filter contains the expected kmers
  ASSERT_TRUE(BF.contains("ACGTACG"));
  ASSERT_TRUE(BF.contains("CGTACGT"));
  ASSERT_TRUE(BF.contains("GTACGTA"));
}

// Test for calculateMinimumScoreThreshold function
TEST_F(BlooMineTest, CalculateMinimumScoreThreshold) {
  std::unordered_set<std::string> target_kset = fasta2kmerset(vm["target_fasta"].as<std::string>(), vm["kmer"].as<int>());
  double spMST = calculateMinimumScoreThreshold(vm, target_kset);

  // Check if the minimum score threshold is calculated correctly
  ASSERT_EQ(spMST, 1.0);
}

// Test for readFQ function
TEST_F(BlooMineTest, ReadFQ) {
  std::vector<std::string> read_data = readFQ(vm["fastq"].as<std::string>());

  // Check if the read data is read correctly
  ASSERT_EQ(read_data.size(), 8);
  ASSERT_EQ(read_data[1], "ACGTACGTACGT");
  ASSERT_EQ(read_data[5], "ATCGATCGATCG");
}

// Test for partitionData function
TEST_F(BlooMineTest, PartitionData) {
  std::vector<std::string> read_data = readFQ(vm["fastq"].as<std::string>());
  std::vector<std::vector<std::string>> read_data_parts = partitionData(read_data, 2);

  // Check if the read data is partitioned correctly
  ASSERT_EQ(read_data_parts.size(), 2);
  ASSERT_EQ(read_data_parts[0].size(), 4);
  ASSERT_EQ(read_data_parts[1].size(), 4);
}

// Test for spawnThreads function (integration test)
TEST_F(BlooMineTest, SpawnThreads) {
  auto [BF, threshold] = generateBloomFilter(vm);
  std::unordered_set<std::string> target_kset = fasta2kmerset(vm["target_fasta"].as<std::string>(), vm["kmer"].as<int>());
  std::vector<std::string> read_data = readFQ(vm["fastq"].as<std::string>());
  std::vector<std::vector<std::string>> read_data_parts = partitionData(read_data, 2);

  // Spawn threads and collect results
  std::vector<std::future<std::vector<std::string>>> hits_collection = spawnThreads(vm, read_data_parts, BF, target_kset, threshold, 1.0);

  // Check if the output file is created and contains the expected reads
  std::ifstream outfile("BM_BMfiltered.fq");
  std::string line;
  std::getline(outfile, line);
  ASSERT_EQ(line, "@read1");
  outfile.close();
}

// Test for reverseCompliment function
TEST_F(BlooMineTest, ReverseCompliment) {
  std::string seq = "ACGTACGTACGT";
  std::string rc_seq = reverseCompliment(seq);

  // Check if the reverse complement is calculated correctly
  ASSERT_EQ(rc_seq, "ACGTACGTACGT");
}

// Test for minscore function
TEST_F(BlooMineTest, Minscore) {
  std::unordered_set<std::string> target_kset = fasta2kmerset(vm["target_fasta"].as<std::string>(), vm["kmer"].as<int>());
  double score = minscore(target_kset.size(), vm["SP-error"].as<float>(), vm["kmer"].as<int>(), hit, gap_open, neg);

  // Check if the minimum score is calculated correctly
  ASSERT_EQ(score, 1.0);
}

// Test for fasta2kmerset function
TEST_F(BlooMineTest, Fasta2Kmerset) {
  std::unordered_set<std::string> kmer_set = fasta2kmerset(vm["target_fasta"].as<std::string>(), vm["kmer"].as<int>());

  // Check if the kmer set is generated correctly
  ASSERT_EQ(kmer_set.size(), 3);
  ASSERT_TRUE(kmer_set.contains("ACGTACG"));
  ASSERT_TRUE(kmer_set.contains("CGTACGT"));
  ASSERT_TRUE(kmer_set.contains("GTACGTA"));
}
