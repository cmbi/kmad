#ifndef VECUTIL_H
#define VECUTIL_H 

#include "types.h"
#include<iostream>
#include<vector>
#include<string>


class Residue;
namespace vecUtil{
  ///
  /// checks if vec contains string x
  ///
	bool contains(FeatureNamesList& vec, std::string& x);
  ///
  /// returns index of the feature in the list of features
  ///
	int findIndex(std::string& val, FeatureNamesList& vec);
  ///
  /// transposes a vector
  ///
	void transposeVec(ProfileMatrix& vec);
  ///
  /// divides vector by an integer
  ///
	void divideVectorByAScalar(ProfileMatrixRow& vec, int scalar);
  ///
  /// divides vector by a double
  ///
	void divideVectorByAScalar(ProfileMatrixRow& vec, double& scalar);
  ///
  /// mutliplies vector by a double
  ///
	void multiplyVectorByAScalar(ProfileMatrixRow& vec, double& scalar);
  ///
  /// adds up scores in every column (result from each column is an element of
  /// the output vector)
  ///
	ProfileMatrixRow addUp(Matrix2D & vec);
	ProfileMatrixColumn convertIntVectorToDoubleVector(SbstMatColumn&);
  ///
  /// converts sequences of Residue object to string sequences
  ///
	StringSequences flatten(const SequenceList&);
  ///
  /// adds an element to the front of the vector
  ///
	ResidueSequence push_front(ResidueSequence& seq, Residue newRes);
  ///
  /// calcualtes an average for each column in the given vec
  ///
	ProfileMatrixRow average(const SbstMatrixColumns& vec);
}

#endif /* VECUTIL_H */
