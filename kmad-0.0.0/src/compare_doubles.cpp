#include "compare_doubles.h"

#include <cmath>
#include <iostream>


bool compare_doubles::is_equal(double one, double two) {
 // allowed difference between the two doubles
 double epsilon = 1e-07;
 if (std::abs(one - two) < epsilon) {
   return true;
 }
 else {
   return false;
 }
}
