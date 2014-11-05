#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main


#include "txtProc.h"

#include <boost/test/unit_test.hpp>

#include <sstream>


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_read_fasta_single)
{
  std::stringstream ss;
  ss << ">test" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;

  auto proteins = read_fasta(ss);

  BOOST_CHECK_EQUAL(proteins.size(), 1);
}

BOOST_AUTO_TEST_CASE(test_read_fasta_multiple)
{
  std::stringstream ss;
  ss << ">test1" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;
  ss << ">test2" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;
  ss << ">test3" << std::endl;
  ss << "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN" << std::endl;
  auto proteins = read_fasta(ss);

  BOOST_CHECK_EQUAL(proteins.size(), 3);
}

BOOST_AUTO_TEST_SUITE_END()
