#ifndef VECUTIL_H
#define VECUTIL_H


#include "fasta.h"
#include "types.h"

#include <iostream>
#include <vector>
#include <string>


namespace vec_util{
  ProfileMatrixColumn ConvertIntVecToDoubleVec(SbstMatColumn& vec);
}

#endif /* VECUTIL_H */
