#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "profile.h"
#include "types.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>


BOOST_AUTO_TEST_SUITE(test_profile)

BOOST_AUTO_TEST_CASE(test_create_profile)
{
  fasta::Sequence s1 = fasta::make_sequence(
      "d", "AzzzzzzMzzzzzzEzzzzzzLzzzzzz", 7);
  fasta::Sequence s2 = fasta::make_sequence(
      "d", "AzzzzzzEzzzzzzEzzzzzzKzzzzzz", 7);
  fasta::Sequence s3 = fasta::make_sequence(
      "d", "AzzzzzzMzzzzzzEzzzzzzLzzzzzz", 7);
  fasta::Sequence s4 = fasta::make_sequence(
      "d", "MzzzzzzMzzzzzzKzzzzzzLzzzzzz", 7);
  fasta::Sequence s5 = fasta::make_sequence(
      "d", "AzzzzzzMzzzzzzAzzzzzzLzzzzzz", 7);

  fasta::SequenceList sequences = {s1, s2, s3, s4, s5};

  ProfileMap p = create_profile(sequences);
  BOOST_CHECK_EQUAL(p.size(), 20);
  BOOST_CHECK_EQUAL(p['A'][0], 4);
  BOOST_CHECK_EQUAL(p['A'][1], 0);
  BOOST_CHECK_EQUAL(p['A'][2], 1);
  BOOST_CHECK_EQUAL(p['A'][3], 0);
  BOOST_CHECK_EQUAL(p['M'][0], 1);
  BOOST_CHECK_EQUAL(p['M'][1], 4);
  BOOST_CHECK_EQUAL(p['M'][2], 0);
  BOOST_CHECK_EQUAL(p['M'][3], 0);
  BOOST_CHECK_EQUAL(p['K'][0], 0);
  BOOST_CHECK_EQUAL(p['K'][1], 0);
  BOOST_CHECK_EQUAL(p['K'][2], 1);
  BOOST_CHECK_EQUAL(p['K'][3], 1);
  BOOST_CHECK_EQUAL(p['E'][0], 0);
  BOOST_CHECK_EQUAL(p['E'][1], 1);
  BOOST_CHECK_EQUAL(p['E'][2], 3);
  BOOST_CHECK_EQUAL(p['E'][3], 0);
  BOOST_CHECK_EQUAL(p['L'][0], 0);
  BOOST_CHECK_EQUAL(p['L'][1], 0);
  BOOST_CHECK_EQUAL(p['L'][2], 0);
  BOOST_CHECK_EQUAL(p['L'][3], 4);
}


BOOST_AUTO_TEST_CASE(test_create_score_profile)
{
  fasta::Sequence s1 = fasta::make_sequence(
      "d", "AzzzzzzMzzzzzzEzzzzzzLzzzzzz", 7);
  fasta::Sequence s2 = fasta::make_sequence(
      "d", "AzzzzzzEzzzzzzEzzzzzzKzzzzzz", 7);
  fasta::Sequence s3 = fasta::make_sequence(
      "d", "AzzzzzzMzzzzzzEzzzzzzLzzzzzz", 7);
  fasta::Sequence s4 = fasta::make_sequence(
      "d", "MzzzzzzMzzzzzzKzzzzzzLzzzzzz", 7);
  fasta::Sequence s5 = fasta::make_sequence(
      "d", "AzzzzzzMzzzzzzAzzzzzzLzzzzzz", 7);

  fasta::SequenceList sequences = {s1, s2, s3, s4, s5};

  ProfileMap p = create_score_profile(sequences);

  BOOST_CHECK_EQUAL(p.size(), 20);
  // TODO: Use BOOST_CHECK to check each value in the score profile to ensure
  //       it's correct.
  BOOST_CHECK_EQUAL(p['A'][0], 0);
  BOOST_CHECK_EQUAL(p['A'][1], 0);
  BOOST_CHECK_EQUAL(p['A'][2], 0);
  BOOST_CHECK_EQUAL(p['A'][3], 0);
  BOOST_CHECK_EQUAL(p['M'][0], 0);
  BOOST_CHECK_EQUAL(p['M'][1], 0);
  BOOST_CHECK_EQUAL(p['M'][2], 0);
  BOOST_CHECK_EQUAL(p['M'][3], 0);
  BOOST_CHECK_EQUAL(p['K'][0], 0);
  BOOST_CHECK_EQUAL(p['K'][1], 0);
  BOOST_CHECK_EQUAL(p['K'][2], 0);
  BOOST_CHECK_EQUAL(p['K'][3], 0);
  BOOST_CHECK_EQUAL(p['E'][0], 0);
  BOOST_CHECK_EQUAL(p['E'][1], 0);
  BOOST_CHECK_EQUAL(p['E'][2], 0);
  BOOST_CHECK_EQUAL(p['E'][3], 0);
  BOOST_CHECK_EQUAL(p['L'][0], 0);
  BOOST_CHECK_EQUAL(p['L'][1], 0);
  BOOST_CHECK_EQUAL(p['L'][2], 0);
  BOOST_CHECK_EQUAL(p['L'][3], 0);
}

BOOST_AUTO_TEST_SUITE_END()
