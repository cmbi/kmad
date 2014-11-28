#include "vecUtil.h"
#include "Residue.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


//checks if the vector of strings vec contains the string x
bool vecUtil::contains(featureNamesList& vec, std::string& x){
	if (std::find(vec.begin(),vec.end(),x) != vec.end()) return true;
	else return false;
}


void vecUtil::transposeVec(profile_matrix& vec){
  profile_matrix newVec;
  profileMatrixRow newRow;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		for (unsigned int j = 0; j < vec.size(); j++){
			newRow.push_back(vec[j][i]);	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
}


void vecUtil::divideVectorByAScalar(profileMatrixRow& vec, int scalar){
  profileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vecUtil::divideVectorByAScalar(profileMatrixRow& vec, double& scalar){
  profileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vecUtil::multiplyVectorByAScalar(profileMatrixRow& vec, double& scalar){
  profileMatrixRow result;
  for (auto &item: vec){
		result.push_back(item*scalar);
	}
	vec = result;
}


//function addUp - takes 2D matrix, adds up elements from each column, returns a 1D vector
profileMatrixRow vecUtil::addUp(matrix2d& vec){
  profileMatrixRow newVec;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for (unsigned int j = 0; j < vec.size(); j++){
			sum += vec[j][i];
		}
		newVec.push_back(sum);
	}
	return newVec;
}


//convertIntVectorToDoubleVector
profileMatrixColumn vecUtil::convertIntVectorToDoubleVector(sbstMatColumn& vec){
  profileMatrixColumn result;
  for (auto &item: vec){
		result.push_back(double(item));
	}
	return result;
}

//flatten a vector of vectors of residues to a vector of strings
string_sequences vecUtil::flatten(const sequenceList& vec){
	string_sequences result;
  for (auto &row: vec){
		std::string newSeq = "";
    for(auto &item: row){
			newSeq += item.getCodon();
		}
		result.push_back(newSeq);
	}
	return result;
}

//add an element to the beginning of the vector
sequence vecUtil::push_front(sequence& seq, Residue newElement){
	reverse(seq.begin(),seq.end());
	seq.push_back(newElement);
	reverse(seq.begin(),seq.end());
	return seq;
}


//calculate average for every column in vector 
profileMatrixRow vecUtil::average(const sbst_matrix_columns& vec){
  profileMatrixRow result;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for(unsigned int j = 0; j < vec.size();j++){
			sum += vec[j][i];
		}
		result.push_back(sum/vec.size());
	}
	return result;
}


//returns index of the first occurence of val in vec
int vecUtil::findIndex(std::string& val, featureNamesList& vec){
	int res = -1;
	for (unsigned int i = 0; i < vec.size(); i++){
		if (vec[i] == val){
			res = i;
			break;
		}
	}
	return res;
}
