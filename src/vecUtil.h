#ifndef VECUTIL_H
#define VECUTIL_H 

#include "types.h"
#include<iostream>
#include<vector>
#include<string>
class Residue;
namespace vecUtil{
	bool contains(featureNamesList& vec, std::string& x);
	int findIndex(std::string& val, featureNamesList& vec);
	void transposeVec(profile_matrix& vec);
	void divideVectorByAScalar(profileMatrixRow& vec, int scalar);
	void divideVectorByAScalar(profileMatrixRow& vec, double& scalar);
	void multiplyVectorByAScalar(profileMatrixRow& vec, double& scalar);
	profileMatrixRow addUp(matrix2d & vec);
	profileMatrixColumn convertIntVectorToDoubleVector(sbstMatColumn&);
	string_sequences flatten(const sequenceList&);
	sequence push_front(sequence&, Residue);
	profileMatrixRow average(const sbst_matrix_columns& vec);
}

#endif /* VECUTIL_H */
