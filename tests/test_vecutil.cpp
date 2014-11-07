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
  int vec_sum = vecUtil::sum(vec);

  BOOST_CHECK_EQUAL(vec_sum, 4);
}


MOCK_CLASS(mock_class){
  MOCK_METHOD(function_to_mock, 0, int(), meth1)
};


BOOST_AUTO_TEST_CASE(test_function_to_mock){
  mock_class c;
  MOCK_EXPECT(c.meth1).returns(4);
}

BOOST_AUTO_TEST_SUITE_END()
