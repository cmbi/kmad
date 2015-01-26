#define BOOST_TEST_DYN_LINK

#include "feature_analysis.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>


BOOST_AUTO_TEST_SUITE(test_feature_analysis)


BOOST_AUTO_TEST_CASE(test_parse_mapfile) {
  BOOST_CHECK_THROW(feature_analysis::parse_mapfile(
        "tests/nonexistent.map"), std::invalid_argument);
  BOOST_CHECK_THROW(feature_analysis::parse_mapfile("tests/wrong_format.map"),
                    std::runtime_error);
  feature_analysis::CodesMap cm = feature_analysis::parse_mapfile(
      "tests/test.map");
  feature_analysis::CodesMap expected_cm = {{"motif_aa",
                                             {"SOME_MOTIF_NAME", "KR."}},
                                            {"domain_aa",
                                             {"SOME_DOMAIN_NAME"}},
                                            {"domain_ab",
                                             {"SOME_DOMAIN_NAME_2"}}};
  BOOST_CHECK_EQUAL(expected_cm.size(), cm.size());
  for (auto feat_it = expected_cm.begin(); feat_it != expected_cm.end();
       ++feat_it) {
    BOOST_CHECK(cm.find(feat_it->first) != cm.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(cm[feat_it->first].begin(),
                                  cm[feat_it->first].end(),
                                  feat_it->second.begin(),
                                  feat_it->second.end());
  }
}


BOOST_AUTO_TEST_SUITE_END()
