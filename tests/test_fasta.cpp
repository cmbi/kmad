#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main

#include "fasta.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_parse_fasta)
{
  fasta::FastaData fd = fasta::parse_fasta("tests/TAU_SPECI.fasta.7c", 7);

  BOOST_CHECK_EQUAL(fd.sequences.size(), 19);
  BOOST_CHECK_EQUAL(fd.probabilities.size(), 97);
}

BOOST_AUTO_TEST_SUITE_END()
