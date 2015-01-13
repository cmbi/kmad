#define BOOST_TEST_DYN_LINK

#include "fasta.h"
#include "f_config.h"
#include "seq_data.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <string>


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_config)
{

  f_config::FeatureSettingsMap test_map;
  f_config::FeatureSettings test_settings;
  f_config::FeaturePositions positions1;
  positions1.seq_no = 1;
  positions1.positions = {1, 2};
  f_config::FeaturePositions positions2;
  positions2.seq_no = 2;
  positions2.seq_no = {0};
  test_settings.positions = {positions1, positions2};


  // SEQUENCE S1
  fasta::Residue r1_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa"}));
  fasta::Residue r1_2("MAAAAAA", std::vector<std::string>({}));
  fasta::Residue r1_3("EAaadaa", std::vector<std::string>({"domain_aa"}));
  fasta::Residue r1_4("LAAAAAA", std::vector<std::string>({}));
  // SEQUENCE S2
  fasta::Residue r2_1("AAAAdaa", std::vector<std::string>({"ptm_phosphP",
                                                          "motif_aa"}));
  fasta::Residue r2_2("EAAAZAA", std::vector<std::string>({"ptm_phosph0"}));
  fasta::Residue r2_3("EAaaZAA", std::vector<std::string>({"ptm_phosph0",
                                                           "domain_aa"}));
  fasta::Residue r2_4("KAAAZAA", std::vector<std::string>({"ptm_phosph0"}));







  // std::string feat_name = "USR_feature1";
  // BOOST_CHECK(test_result.find(feat_name) != test_result.end());
  // BOOST_CHECK_EQUAL(test_result[feat_name].add_score, 5);
  // BOOST_CHECK_EQUAL(test_result[feat_name].subtract_score, 4);

  // std::vector<std::string> expected_list = {"USR_feature1", "USR_feature2",
  //                                           "USR_feature4"};
  // std::vector<std::string> result_list = test_result[feat_name].add_features;
  // BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
  //                               result_list.begin(), result_list.end());

  // expected_list = {"USR_feature3", "USR_feature6"};
  // result_list = test_result[feat_name].subtract_features;
  // BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
  //                               result_list.begin(), result_list.end());

  // std::vector<int> expected_int_list = {2, 3};
  // std::vector<int> result_int_list;
  // result_int_list = test_result[feat_name].positions[0].positions;
  // BOOST_CHECK_EQUAL_COLLECTIONS(expected_int_list.begin(),
  //                               expected_int_list.end(),
  //                               result_int_list.begin(),
  //                               result_int_list.end());
}


BOOST_AUTO_TEST_SUITE_END()
