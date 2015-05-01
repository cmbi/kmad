#define BOOST_TEST_DYN_LINK

#include "src/fasta.h"
#include "src/profile.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE(test_profile)

namespace {
  static const std::vector<char> ALPHABET = {'A','R','N','D','C','Q','E','G',
                                             'H','I','L','K','M','F','P','S',
                                             'T','W','Y','V', 'B', 'Z', 'X'};
}

BOOST_AUTO_TEST_CASE(test_create_profile)
{
  fasta::Sequence s1 = fasta::make_sequence(
      "d", "AAAAAAAMAAAAAAEAAAAAALAAAAAABAAAAAAZAAAAAAXAAAAAA", 7);
  fasta::Sequence s2 = fasta::make_sequence(
      "d", "AAAAAAAEAAAAAAEAAAAAAKAAAAAABAAAAAAZAAAAAALAAAAAA", 7);
  fasta::Sequence s3 = fasta::make_sequence(
      "d", "AAAAAAAMAAAAAAEAAAAAALAAAAAABAAAAAAZAAAAAALAAAAAA", 7);
  fasta::Sequence s4 = fasta::make_sequence(
      "d", "MAAAAAAMAAAAAAKAAAAAALAAAAAANAAAAAAQAAAAAALAAAAAA", 7);
  fasta::Sequence s5 = fasta::make_sequence(
      "d", "AAAAAAAMAAAAAAAAAAAAALAAAAAALAAAAAALAAAAAALAAAAAA", 7);

  fasta::SequenceList sequences = {s1, s2, s3, s4, s5};

  profile::ProfileMap p = profile::create_profile(sequences);
  profile::ProfileMap expected_profile = {{'A', {4, 0, 1, 0,   0,   0, 0.05}}, 
                                          {'R', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'N', {0, 0, 0, 0, 2.5,   0, 0.05}},
                                          {'D', {0, 0, 0, 0, 1.5,   0, 0.05}},
                                          {'C', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'Q', {0, 0, 0, 0,   0, 2.5, 0.05}},
                                          {'E', {0, 1, 3, 0,   0, 1.5, 0.05}},
                                          {'G', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'H', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'I', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'L', {0, 0, 0, 4,   1,   1, 4.05}},
                                          {'K', {0, 0, 1, 1,   0,   0, 0.05}},
                                          {'M', {1, 4, 0, 0,   0,   0, 0.05}},
                                          {'F', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'P', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'S', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'T', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'W', {0, 0, 0, 0,   0,   0, 0.05}},
                                          {'Y', {0, 0, 0, 0,   0,   0, 0.05}}, 
                                          {'V', {0, 0, 0, 0,   0,   0, 0.05}}};
  BOOST_CHECK_EQUAL(p.size(), expected_profile.size());
  for (auto aa_it  = p.begin(); aa_it != p.end(); ++aa_it) {
    for (size_t i = 0; i < aa_it->second.size(); ++i) {
      BOOST_CHECK(std::abs(aa_it->second[i]
            - expected_profile[aa_it->first][i]) < 0.001);
    }
  }
}


BOOST_AUTO_TEST_CASE(test_create_score_profile)
{
  fasta::Sequence s1 = fasta::make_sequence(
      "d", "AAAAAAAMAAAAAAEAAAAAALAAAAAA", 7);
  fasta::Sequence s2 = fasta::make_sequence(
      "d", "AAAAAAAEAAAAAAEAAAAAAKAAAAAA", 7);
  fasta::Sequence s3 = fasta::make_sequence(
      "d", "AAAAAAAMAAAAAAEAAAAAALAAAAAA", 7);
  fasta::Sequence s4 = fasta::make_sequence(
      "d", "MAAAAAAMAAAAAAKAAAAAALAAAAAA", 7);
  fasta::Sequence s5 = fasta::make_sequence(
      "d", "AAAAAAAMAAAAAAAAAAAAALAAAAAA", 7);

  fasta::SequenceList sequences = {s1, s2, s3, s4, s5};
  std::string sbst_mat = "BLOSUM";
  profile::ProfileMap p = profile::create_score_profile(sequences, sbst_mat);

  profile::ProfileMap expected_score_profile; 
  expected_score_profile = {{'A', { 3.0, -1.0,  0.0, -1.0}}, 
                            {'R', {-1.0, -0.8,  0.2, -1.2}},
                            {'N', {-2.0, -1.6, -0.4, -2.4}},
                            {'D', {-2.2, -2.0,  0.6, -3.4}},
                            {'C', {-0.2, -1.6, -3.0, -1.4}},
                            {'Q', {-0.8,  0.4,  1.2, -1.4}},
                            {'E', {-1.2, -0.6,  3.0, -2.2}},
                            {'G', {-0.6, -2.8, -1.6, -3.6}},
                            {'H', {-2.0, -1.6, -0.6, -2.6}},
                            {'I', {-0.6,  0.2, -2.6,  1.0}},
                            {'L', {-0.4,  1.0, -2.4,  2.8}},
                            {'K', {-1.0, -0.6,  1.4, -0.6}},
                            {'M', { 0.2,  3.6, -1.6,  1.4}},
                            {'F', {-1.6, -0.6, -2.8, -0.6}},
                            {'P', {-1.2, -1.8, -1.0, -2.6}},
                            {'S', { 0.6, -0.8,  0.2, -1.6}},
                            {'T', {-0.2, -1.0, -0.8, -1.0}},
                            {'W', {-2.6, -1.4, -3.0, -2.2}},
                            {'Y', {-1.8, -1.2, -2.0, -1.2}}, 
                            {'V', { 0.2,  0.4, -1.6,  0.4}},
                            {'B', {-2.1, -1.8,  0.1, -2.9}},
                            {'Z', {-1.0, -0.1,  2.1, -1.8}},
                            {'X', {-0.77, -0.69, -0.84, -1.17}}};

  BOOST_CHECK_EQUAL(p.size(), expected_score_profile.size());
  for (auto& aa: ALPHABET){
    for (size_t i = 0; i < p[aa].size(); ++i) {
      BOOST_CHECK(std::abs(p[aa][i] - expected_score_profile[aa][i]) < 0.001);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
