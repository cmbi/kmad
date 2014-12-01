#include "vec_util.h"
#include "residue.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


bool vec_util::contains(FeatureNamesList& vec, std::string& x){
	if (std::find(vec.begin(),vec.end(),x) != vec.end()) return true;
	else return false;
}


void vec_util::transposeVec(ProfileMatrix& vec){
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


void vec_util::divideVectorByAScalar(ProfileMatrixRow& vec, int scalar){
  ProfileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vec_util::divideVectorByAScalar(ProfileMatrixRow& vec, double& scalar){
  ProfileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vec_util::multiplyVectorByAScalar(ProfileMatrixRow& vec, double& scalar){
  ProfileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item*scalar);
	}
	vec = result;
}


ProfileMatrixRow vec_util::addUp(Matrix2D& vec){
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


ProfileMatrixColumn vec_util::convertIntVectorToDoubleVector(SbstMatColumn& vec){
  ProfileMatrixColumn result;
  for (auto &item: vec){
		result.push_back(double(item));
	}
	return result;
}


StringSequences vec_util::flatten(const SequenceList& vec){
	StringSequences result;
  for (auto &row: vec){
		std::string newSeq = "";
    for(auto &item: row){
			newSeq += item.get_codon();
		}
		result.push_back(newSeq);
	}
	return result;
}


ResidueSequence vec_util::push_front(ResidueSequence& seq, Residue newRes){
	reverse(seq.begin(),seq.end());
	seq.push_back(newRes);
	reverse(seq.begin(),seq.end());
	return seq;
}


ProfileMatrixRow vec_util::average(const SbstMatrixColumns& vec){
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


int vec_util::findIndex(std::string& val, FeatureNamesList& vec){
	int res = -1;
	for (unsigned int i = 0; i < vec.size(); i++){
		if (vec[i] == val){
			res = i;
			break;
		}
	}
	return res;
}
