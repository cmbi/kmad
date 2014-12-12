#ifndef VECUTIL_H
#define VECUTIL_H

#include "types.h"

#include<iostream>
#include<vector>
#include<string>


class Residue;
namespace vec_util{
  ///
  /// returns index of the feature in the list of features
  ///
  int FindIndex(std::string& val, FeatureNamesList& vec);
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
  StringSequences Flatten(const SequenceList&);
  ///
  /// calcualtes an average for each column in the given vec
  ///
  ProfileMatrixRow Average(const SbstMatrixColumns& vec);
}

#endif /* VECUTIL_H */
