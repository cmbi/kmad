#define BOOST_TEST_DYN_LINK

#include "compare_doubles.h"

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(test_compare_doubles)

BOOST_AUTO_TEST_CASE(test_is_equal) {
  double one = 1;
  double two = 1.05;
  BOOST_CHECK(!compare_doubles::is_equal(one, two));
  one = two + 1e-08;
  BOOST_CHECK(compare_doubles::is_equal(one, two));
}

BOOST_AUTO_TEST_SUITE_END()
