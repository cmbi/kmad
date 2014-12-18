#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "profile.h"
#include "types.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <map>
#include <vector>


BOOST_AUTO_TEST_SUITE(test_profile)

namespace {
  static const std::vector<char> ALPHABET = {'A','R','N','D','C','Q','E','G',
                                             'H','I','L','K','M','F','P','S',
                                             'T','W','Y','V'};
}

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

   
  ProfileMap expected_score_profile = {{'A', {3.0, -1.0, 0.0, -1.0}}, 
                                       {'R', {-1.0, -0.8, 0.2, -1.2}},
                                       {'N', {-2.0, -1.6, -0.4, -2.4}},
                                       {'D', {-2.2, -2.0, 0.6, -3.4}},
                                       {'C', {-0.2, -1.6, -3.0, -1.4}},
                                       {'Q', {-0.8, 0.4, 1.2, -1.4}},
                                       {'E', {-1.2, -0.6, 3.0, -2.2}},
                                       {'G', {-0.6, -2.8, -1.6, -3.6}},
                                       {'H', {-2.0, -1.6, -0.6, -2.6}},
                                       {'I', {-0.6, 0.2, -2.6, 1.0}},
                                       {'L', {-0.4, 1.0, -2.4, 2.8}},
                                       {'K', {-1.0, -0.6, 1.4, -0.6}},
                                       {'M', {0.2, 3.6, -1.6, 1.4}},
                                       {'F', {-1.6, -0.6, -2.8, -0.6}},
                                       {'P', {-1.2, -1.8, -1.0, -2.6}},
                                       {'S', {0.6, -0.8, 0.2, -1.6}},
                                       {'T', {-0.2, -1.0, -0.8, -1.0}},
                                       {'W', {-2.6, -1.4, -3.0, -2.2}},
                                       {'Y', {-1.8, -1.2, -2.0, -1.2}}, 
                                       {'V', {0.2, 0.4, -1.6, 0.4}}};

  BOOST_CHECK_EQUAL(p.size(), 20);
  // // TODO: Use BOOST_CHECK to check each value in the score profile to ensure
  // //       it's correct.
  for (auto& aa: ALPHABET){
    BOOST_CHECK_EQUAL_COLLECTIONS(p[aa].begin(), p[aa].end(), 
                                  expected_score_profile[aa].begin(),
                                  expected_score_profile[aa].end()); 
  }
}

BOOST_AUTO_TEST_SUITE_END()
