#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "features_profile.h"
#include "types.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE(test_features_profile)

BOOST_AUTO_TEST_CASE(test_create_features_profile)
{
  // fasta::Sequence s1 = fasta::make_sequence(
  //     "d", "AAAAdaaMAAAAAAEAaaAAALAAAAAA", 7);
  // fasta::Sequence s2 = fasta::make_sequence(
  //     "d", "AAAAdaaEAAAZAAEAaaZAAKAAAZAA", 7);
  // fasta::Sequence s3 = fasta::make_sequence(
  //     "d", "AAAAZaaMAAAZAAEAacZAALAaaZaa", 7);
  // fasta::Sequence s4 = fasta::make_sequence(
  //     "d", "MAAAZabMAAAZAAKAaaZAALAAAAaa", 7);
  // fasta::Sequence s5 = fasta::make_sequence(
  //     "d", "AAAAAacMAAAZAAAAAAAAALAAAAAA", 7);
  fasta::Residue r1_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa"}));
  fasta::Residue r1_2("MAAAAAA", std::vector<std::string>({}));
  fasta::Residue r1_3("EAaadaa", std::vector<std::string>({"domain_aa"}));
  fasta::Residue r1_4("LAAAAAA", std::vector<std::string>({}));

  fasta::Residue r2_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa"}));
  fasta::Residue r2_2("EAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue r2_3("EAaaZAA", std::vector<std::string>({"ptm_phosph0",
                                                           "domain_aa"}));
  fasta::Residue r2_4("KAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue r3_1("AAAAZaa", std::vector<std::string>({"ptm_phosph0",
                                                           "motif_aa"}));
  fasta::Residue r3_2("MAAAZAA", std::vector<std::string>({"ptm_phosph0"}));


  fasta::SequenceList sequences = {s1, s2, s3, s4, s5};

  ProfileMap p = create_features_profile(sequences);
  // BOOST_CHECK_EQUAL(p.size(), 20);
  // BOOST_CHECK_EQUAL(p['A'][0], 4);
  // BOOST_CHECK_EQUAL(p['A'][1], 0);
  // BOOST_CHECK_EQUAL(p['A'][2], 1);
  // BOOST_CHECK_EQUAL(p['A'][3], 0);
  // BOOST_CHECK_EQUAL(p['M'][0], 1);
  // BOOST_CHECK_EQUAL(p['M'][1], 4);
  // BOOST_CHECK_EQUAL(p['M'][2], 0);
  // BOOST_CHECK_EQUAL(p['M'][3], 0);
  // BOOST_CHECK_EQUAL(p['K'][0], 0);
  // BOOST_CHECK_EQUAL(p['K'][1], 0);
  // BOOST_CHECK_EQUAL(p['K'][2], 1);
  // BOOST_CHECK_EQUAL(p['K'][3], 1);
  // BOOST_CHECK_EQUAL(p['E'][0], 0);
  // BOOST_CHECK_EQUAL(p['E'][1], 1);
  // BOOST_CHECK_EQUAL(p['E'][2], 3);
  // BOOST_CHECK_EQUAL(p['E'][3], 0);
  // BOOST_CHECK_EQUAL(p['L'][0], 0);
  // BOOST_CHECK_EQUAL(p['L'][1], 0);
  // BOOST_CHECK_EQUAL(p['L'][2], 0);
  // BOOST_CHECK_EQUAL(p['L'][3], 4);
}


// BOOST_AUTO_TEST_CASE(test_create_score_profile)
// {
//   fasta::Sequence s1 = fasta::make_sequence(
//       "d", "AzzzzzzMzzzzzzEzzzzzzLzzzzzz", 7);
//   fasta::Sequence s2 = fasta::make_sequence(
//       "d", "AzzzzzzEzzzzzzEzzzzzzKzzzzzz", 7);
//   fasta::Sequence s3 = fasta::make_sequence(
//       "d", "AzzzzzzMzzzzzzEzzzzzzLzzzzzz", 7);
//   fasta::Sequence s4 = fasta::make_sequence(
//       "d", "MzzzzzzMzzzzzzKzzzzzzLzzzzzz", 7);
//   fasta::Sequence s5 = fasta::make_sequence(
//       "d", "AzzzzzzMzzzzzzAzzzzzzLzzzzzz", 7);
// 
//   fasta::SequenceList sequences = {s1, s2, s3, s4, s5};
// 
//   ProfileMap p = create_score_profile(sequences);
// 
//    
//   ProfileMap expected_score_profile = {{'A', {3.0, -1.0, 0.0, -1.0}}, 
//                                        {'R', {-1.0, -0.8, 0.2, -1.2}},
//                                        {'N', {-2.0, -1.6, -0.4, -2.4}},
//                                        {'D', {-2.2, -2.0, 0.6, -3.4}},
//                                        {'C', {-0.2, -1.6, -3.0, -1.4}},
//                                        {'Q', {-0.8, 0.4, 1.2, -1.4}},
//                                        {'E', {-1.2, -0.6, 3.0, -2.2}},
//                                        {'G', {-0.6, -2.8, -1.6, -3.6}},
//                                        {'H', {-2.0, -1.6, -0.6, -2.6}},
//                                        {'I', {-0.6, 0.2, -2.6, 1.0}},
//                                        {'L', {-0.4, 1.0, -2.4, 2.8}},
//                                        {'K', {-1.0, -0.6, 1.4, -0.6}},
//                                        {'M', {0.2, 3.6, -1.6, 1.4}},
//                                        {'F', {-1.6, -0.6, -2.8, -0.6}},
//                                        {'P', {-1.2, -1.8, -1.0, -2.6}},
//                                        {'S', {0.6, -0.8, 0.2, -1.6}},
//                                        {'T', {-0.2, -1.0, -0.8, -1.0}},
//                                        {'W', {-2.6, -1.4, -3.0, -2.2}},
//                                        {'Y', {-1.8, -1.2, -2.0, -1.2}}, 
//                                        {'V', {0.2, 0.4, -1.6, 0.4}}};
// 
//   BOOST_CHECK_EQUAL(p.size(), expected_score_profile.size());
//   for (auto& aa: ALPHABET){
//     for (size_t i = 0; i < p['A'].size(); ++i) {
//       BOOST_CHECK(std::abs(p[aa][i] - expected_score_profile[aa][i]) < 0.1);
//     }
//   }
// }

BOOST_AUTO_TEST_SUITE_END()
