#define BOOST_TEST_DYN_LINK

#include "f_config.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <map>
#include <vector>
#include <algorithm>
#include <string>


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_config)
{
  const char* filename = "tests/test_conf_file.cfg";

  f_config::FeatureSettingsMap test_result;
  test_result = f_config::ConfParser::parse_conf_file(filename);

  std::string feat_name = "feature1";
  BOOST_CHECK(test_result.find(feat_name) != test_result.end());
  BOOST_CHECK_EQUAL(test_result[feat_name].tag, "groupA");
  BOOST_CHECK_EQUAL(test_result[feat_name].add_score, 5);
  BOOST_CHECK_EQUAL(test_result[feat_name].subtract_score, 4);

  std::vector<std::string> expected_list = {"feature1", "feature2"};
  std::vector<std::string> result_list = test_result[feat_name].add_features;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());

  expected_list = {"groupA"};
  result_list = test_result[feat_name].add_tags;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());

  expected_list = {"feature3"};
  result_list = test_result[feat_name].add_exceptions;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());

  expected_list = {"feature3"};
  result_list = test_result[feat_name].subtract_features;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());


  expected_list = {"groupB"};
  result_list = test_result[feat_name].subtract_tags;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());

  expected_list = {"feature2"};
  result_list = test_result[feat_name].subtract_exceptions;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());

  std::vector<int> expected_int_list = {2, 3};
  std::vector<int> result_int_list;
  result_int_list = test_result[feat_name].positions[0].positions;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_int_list.begin(),
                                expected_int_list.end(),
                                result_int_list.begin(),
                                result_int_list.end());





/*
  BOOST_CHECK_EQUAL(sequences.size(), 2);
  BOOST_CHECK_EQUAL(single_codon.first, expected_codon.first);
  BOOST_CHECK_EQUAL_COLLECTIONS(single_codon.second.begin(),
                                single_codon.second.end(),
                                expected_codon.second.begin(),
                                expected_codon.second.end());

  BOOST_CHECK_EQUAL(prob_map.size(), 1);
  BOOST_CHECK_EQUAL(prob_map["ab"], 0.7);
  */

}


BOOST_AUTO_TEST_SUITE_END()
