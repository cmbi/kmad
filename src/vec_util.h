#ifndef VECUTIL_H
#define VECUTIL_H


#include "fasta.h"
#include "types.h"

#include <iostream>
#include <vector>
#include <string>


namespace vec_util{
  ///
  /// transposes a vector
  ///
  void TransposeVec(ProfileMatrix& vec);
  ///
  /// divides vector by an integer
  ///
  void DivideVectorByAScalar(ProfileMatrixRow& vec, int scalar);
  ///
  /// divides vector by a double
  ///
  void DivideVectorByAScalar(ProfileMatrixRow& vec, double& scalar);
  ///
  /// mutliplies vector by a double
  ///
  void MultiplyVectorByAScalar(ProfileMatrixRow& vec, double& scalar);
  ///
  /// adds up scores in every column (result from each column is an element of
  /// the output vector)
  ///
  ProfileMatrixRow AddUp(Matrix2D& vec);
  ProfileMatrixColumn ConvertIntVecToDoubleVec(SbstMatColumn& vec);
  ///
  /// converts sequences of Residue object to string sequences
  ///
  std::vector<std::string> Flatten(const std::vector<fasta::Sequence>&);
  ///
  /// calcualtes an average for each column in the given vec
  ///
  ProfileMatrixRow Average(const SbstMatrixColumns& vec);
}

#endif /* VECUTIL_H */
