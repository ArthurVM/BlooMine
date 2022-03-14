cmake_minimum_required(VERSION 3.13.4)
project(BlooMine)

## set compiler and flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

## find and include boost program_options library
find_package( Boost 1.40 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

## executables and source files
add_executable(BlooMine ./src/BlooMine.cpp)

target_link_libraries( BlooMine LINK_PUBLIC ${Boost_LIBRARIES} )

## install rules
install(TARGETS BlooMine DESTINATION bin)

## test module
enable_testing()

## run test
add_test(NAME Runs COMMAND BlooMine)
set_tests_properties(Runs
  PROPERTIES PASS_REGULAR_EXPRESSION "Too few argument provided!*")

## help test
add_test(NAME Help COMMAND BlooMine --help)
set_tests_properties(Help
  PROPERTIES PASS_REGULAR_EXPRESSION "Allowed options:*")

## test data run
add_test(NAME test_run COMMAND BlooMine --target_fasta ../test/targets/D614G.fa --fastq ../test/fqs/testFQ.fq --threads 2)
set_tests_properties(test_run
  PROPERTIES PASS_REGULAR_EXPRESSION "*Finished screening.*")
