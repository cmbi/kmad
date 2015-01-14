#include "vec_util.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>


ProfileMatrixColumn vec_util::ConvertIntVecToDoubleVec(SbstMatColumn& vec) {
  ProfileMatrixColumn result;
  for (auto &item : vec) {
    result.push_back(double(item));
  }
  return result;
}
