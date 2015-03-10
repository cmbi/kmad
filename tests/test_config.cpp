#define BOOST_TEST_DYN_LINK

#include "f_config.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace lcg = libconfig;


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_config)
{
  const char* filename = "tests/test_conf_file.cfg";

  f_config::FeatureSettingsMap test_result;
  test_result = f_config::ConfParser::parse_conf_file(filename);

  std::string feat_name = "USR_feature1";
  BOOST_CHECK(test_result.find(feat_name) != test_result.end());
  BOOST_CHECK_EQUAL(test_result[feat_name].add_score, 5);
  BOOST_CHECK_EQUAL(test_result[feat_name].subtract_score, 4);

  std::vector<std::string> expected_list = {"USR_feature1", "USR_feature2",
                                            "USR_feature4"};
  std::vector<std::string> result_list = test_result[feat_name].add_features;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());

  expected_list = {"USR_feature3", "USR_feature6"};
  result_list = test_result[feat_name].subtract_features;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_list.begin(), expected_list.end(),
                                result_list.begin(), result_list.end());

  std::vector<int> expected_int_list = {1, 2};
  std::vector<int> result_int_list;
  result_int_list = test_result[feat_name].positions[0].positions;
  BOOST_CHECK_EQUAL_COLLECTIONS(expected_int_list.begin(),
                                expected_int_list.end(),
                                result_int_list.begin(),
                                result_int_list.end());
  filename = "tests/test_nonexistent.cfg";

  BOOST_CHECK_THROW(f_config::ConfParser::parse_conf_file(filename), 
                    lcg::FileIOException);
  filename = "tests/test_conf_file_wrongformat.cfg";
  BOOST_CHECK_THROW(f_config::ConfParser::parse_conf_file(filename), 
                    lcg::ParseException);
  filename = "tests/test_conf_file_settingnotfound.cfg";
  BOOST_CHECK_THROW(f_config::ConfParser::parse_conf_file(filename), 
                    lcg::SettingNotFoundException);

}


BOOST_AUTO_TEST_SUITE_END()
