#include "vec_util.h"
#include "residue.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


bool vecUtil::contains(FeatureNamesList& vec, std::string& x){
	if (std::find(vec.begin(),vec.end(),x) != vec.end()) return true;
	else return false;
}


void vecUtil::transposeVec(ProfileMatrix& vec){
  ProfileMatrix newVec;
  ProfileMatrixRow newRow;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		for (unsigned int j = 0; j < vec.size(); j++){
			newRow.push_back(vec[j][i]);	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
}


void vecUtil::divideVectorByAScalar(ProfileMatrixRow& vec, int scalar){
  ProfileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vecUtil::divideVectorByAScalar(ProfileMatrixRow& vec, double& scalar){
  ProfileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vecUtil::multiplyVectorByAScalar(ProfileMatrixRow& vec, double& scalar){
  ProfileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item*scalar);
	}
	vec = result;
}


ProfileMatrixRow vecUtil::addUp(Matrix2D& vec){
  ProfileMatrixRow newVec;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for (unsigned int j = 0; j < vec.size(); j++){
			sum += vec[j][i];
		}
		newVec.push_back(sum);
	}
	return newVec;
}


ProfileMatrixColumn vecUtil::convertIntVectorToDoubleVector(SbstMatColumn& vec){
  ProfileMatrixColumn result;
  for (auto &item: vec){
		result.push_back(double(item));
	}
	return result;
}


StringSequences vecUtil::flatten(const SequenceList& vec){
	StringSequences result;
  for (auto &row: vec){
		std::string newSeq = "";
    for(auto &item: row){
			newSeq += item.getCodon();
		}
		result.push_back(newSeq);
	}
	return result;
}


ResidueSequence vecUtil::push_front(ResidueSequence& seq, Residue newRes){
	reverse(seq.begin(),seq.end());
	seq.push_back(newRes);
	reverse(seq.begin(),seq.end());
	return seq;
}


ProfileMatrixRow vecUtil::average(const SbstMatrixColumns& vec){
  ProfileMatrixRow result;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for(unsigned int j = 0; j < vec.size();j++){
			sum += vec[j][i];
		}
		result.push_back(sum/vec.size());
	}
	return result;
}


int vecUtil::findIndex(std::string& val, FeatureNamesList& vec){
	int res = -1;
	for (unsigned int i = 0; i < vec.size(); i++){
		if (vec[i] == val){
			res = i;
			break;
		}
	}
	return res;
}
