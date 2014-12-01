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
	void transposeVec(profile_matrix& vec);
  ///
  /// divides vector by an integer
  ///
	void divideVectorByAScalar(profileMatrixRow& vec, int scalar);
  ///
  /// divides vector by a double
  ///
	void divideVectorByAScalar(profileMatrixRow& vec, double& scalar);
  ///
  /// mutliplies vector by a double
  ///
	void multiplyVectorByAScalar(profileMatrixRow& vec, double& scalar);
  ///
  /// adds up scores in every column (result from each column is an element of
  /// the output vector)
  ///
	profileMatrixRow addUp(matrix2d & vec);
	profileMatrixColumn convertIntVectorToDoubleVector(sbstMatColumn&);
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
	profileMatrixRow average(const sbst_matrix_columns& vec);
}

#endif /* VECUTIL_H */
