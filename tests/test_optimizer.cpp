#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "optimizer.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>


BOOST_AUTO_TEST_SUITE(test_optimizer)

BOOST_AUTO_TEST_CASE(test_find_gap_end)
{
  fasta::Sequence s = fasta::make_sequence("WF---F", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_end(s, 2), 4);

  s = fasta::make_sequence("WF---", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_end(s, 2), 5);

  s = fasta::make_sequence("---WF", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_end(s, 0), 2);

}


BOOST_AUTO_TEST_CASE(test_find_gap_start)
{
  fasta::Sequence s = fasta::make_sequence("WF---WF", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_start(s, 4), 2);

  s = fasta::make_sequence("---WF", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_start(s, 2), -1);

  s = fasta::make_sequence("WF---", 1);
  BOOST_CHECK_EQUAL(optimizer::find_gap_start(s, 4), 2);
}

BOOST_AUTO_TEST_SUITE_END()
