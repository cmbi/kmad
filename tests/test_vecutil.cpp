#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main


#include "vecUtil.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <turtle/mock.hpp>

#include <sstream>

#include <vector>


BOOST_AUTO_TEST_SUITE(test_kman_suite)

BOOST_AUTO_TEST_CASE(test_sum)
{
  std::vector<double> vec = {2., 0., 2., 0.};
}


BOOST_AUTO_TEST_SUITE_END()
