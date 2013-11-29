#include "vecUtil.h"
#include <iostream>
#include <vector>
#include <string>
//function transposeVec
void vecUtil::transposeVec(std::vector< std::vector<int> >& vec){
	std::vector< std::vector<int> > newVec;
	std::vector<int> newRow;
	for (int i = 0; i < vec.at(0).size(); i++){
		for (int j = 0; j < vec.size(); j++){
			newRow.push_back(vec.at(j).at(i));	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
	//return newVec;
}
//function transposeVec  - transposes vector< vector<double> > and returns vector< vector<double> >
void vecUtil::transposeVec(std::vector< std::vector<double> >& vec){
	std::vector< std::vector<double> > newVec;
	std::vector<double> newRow;
	for (int i = 0; i < vec.at(0).size(); i++){
		for (int j = 0; j < vec.size(); j++){
			newRow.push_back(vec.at(j).at(i));	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
}
void vecUtil::divideVectorByAScalar(std::vector<double>& vec, int scalar){
	std::vector<double> result;
	for (int i = 0; i < vec.size(); i++){
		result.push_back(vec[i]/scalar);	
	}
	vec = result;
}
void vecUtil::multiplyVectorByAScalar(std::vector<double>& vec, double scalar){
	std::vector<double> result;
	for (int i = 0; i < vec.size(); i++){
		result.push_back(vec[i]*scalar);
	}
	vec = result;
}
//function addUp - takes 2D matrix, adds up elements from each column, returns a 1D vector
std::vector<double> vecUtil::addUp(std::vector< std::vector<double> > vec){
	std::vector<double> newVec;
	for (int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for (int j = 0; j < vec.size(); j++){
			sum += vec.at(j).at(i);
		}
		newVec.push_back(sum);
	}
	return newVec;
}
//convertIntVectorToDoubleVector
std::vector<double> vecUtil::convertIntVectorToDoubleVector(std::vector<int> vec){
	std::vector<double> result;
	for (int i = 0; i < vec.size();i++){
		result.push_back(double(vec.at(i)));
	}
	return result;
}
//function printDoubleVector
void vecUtil::printDoubleVector(const std::vector<double>& vec){
	for (int i = 0; i < vec.size(); i++){
		std::cout << vec[i];
	}
	std::cout << "\n";
}
std::vector<std::string> vecUtil::flatten(const std::vector<std::vector<std::string> >& vec){
	std::vector<std::string> result;
	for (int i = 0; i < vec.size();i++){
		std::string newSeq = "";
		for(int j = 0; j < vec[i].size(); j++){
			char newChar = vec[i][j][0];
		//	newSeq.append(&newChar);
		//	newSeq.append(&vec[i][j][0]);
			newSeq+=newChar;
			std::cout << vec[i][j][0] << "\n";
		}
		result.push_back(newSeq);
		std::cout << newSeq << "\n";
	}
	return result;
}
std::vector<std::string> vecUtil::push_front(std::vector<std::string> vec, std::string newElement){
	reverse(vec.begin(),vec.end());
	vec.push_back(newElement);
	reverse(vec.begin(),vec.end());
	return vec;
}
