#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Config

#include "f_config.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <map>


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_config)
{
  const char* filename = "tests/test_conf_file.cfg";
  
  f_config::UsrFeatureMap test_result; 
  test_result = f_config::ConfParser::parse_conf_file(filename);

  std::string expected_feat_name = "feature1";
  std::cout << test_result.size() << std::endl;
  BOOST_CHECK(test_result.find(expected_feat_name) != test_result.end());

  

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
