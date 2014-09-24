#include "vecUtil.h"
#include "Residue.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
//checks if the vector of strings vec contains the string x
bool vecUtil::contains(std::vector<std::string>& vec, std::string x){
	if (std::find(vec.begin(),vec.end(),x) != vec.end()) return true;
	else return false;
}
void vecUtil::transposeVec(std::vector< std::vector<int> >& vec){
	std::vector< std::vector<int> > newVec;
	std::vector<int> newRow;
	for (int i = 0; i < vec[0].size(); i++){
		for (int j = 0; j < vec.size(); j++){
			newRow.push_back(vec[j][i]);	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
}
void vecUtil::transposeVec(std::vector< std::vector<double> >& vec){
	std::vector< std::vector<double> > newVec;
	std::vector<double> newRow;
	for (int i = 0; i < vec[0].size(); i++){
		for (int j = 0; j < vec.size(); j++){
			newRow.push_back(vec[j][i]);	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
}
void vecUtil::divideVectorByAScalar(std::vector<double>& vec, double scalar){
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
			sum += vec[j][i];
		}
		newVec.push_back(sum);
	}
	return newVec;
}
//convertIntVectorToDoubleVector
std::vector<double> vecUtil::convertIntVectorToDoubleVector(std::vector<int> vec){
	std::vector<double> result;
	for (int i = 0; i < vec.size();i++){
		result.push_back(double(vec[i]));
	}
	return result;
}
//function printDoubleVector
void vecUtil::printVector(const std::vector<int>& vec){
	for (int i = 0; i < vec.size(); i++){
		std::cout << vec[i] << " ";
	}
	std::cout << "\n";
}
void vecUtil::printVector(const std::vector<double>& vec){
	for (int i = 0; i < vec.size(); i++){
		std::cout << vec[i] << " ";
	}
	std::cout << "\n";
}
void vecUtil::printVector(const std::vector<std::string>& vec){
	for (int i = 0; i < vec.size(); i++){
		std::cout << vec[i] << " ";
	}
	std::cout << "\n";
}
//takes a set of encoded sequences (=vector of vectors of strings) and returns a vector of nonencoded sequences
std::vector<std::string> vecUtil::flattenWithoutFeatures(const std::vector<std::vector<std::string> >& vec){
	std::vector<std::string> result;
	for (int i = 0; i < vec.size();i++){
		std::string newSeq = "";
		for(int j = 0; j < vec[i].size(); j++){
			newSeq+=vec[i][j][0];
		}
		result.push_back(newSeq);
	}
	return result;
}
//flatten a vector of residue vectors to a vector of strings
std::vector<std::string> vecUtil::flatten(const std::vector<std::vector<Residue> > & vec){
	std::vector<std::string> result;
	for (int i = 0; i < vec.size();i++){
		std::string newSeq = "";
		for(int j = 0; j < vec[i].size(); j++){
			newSeq+=vec[i][j].getCodon();
		}
		result.push_back(newSeq);
	}
	return result;
}
std::vector<std::string> vecUtil::flatten(const std::vector<std::vector<std::string> >& vec){
	std::vector<std::string> result;
	for (int i = 0; i < vec.size();i++){
		std::string newSeq = "";
		for(int j = 0; j < vec[i].size(); j++){
			newSeq+=vec[i][j];
		}
		result.push_back(newSeq);
	}
	return result;
}
//add an element to the beginning of the vector
std::vector<Residue> vecUtil::push_front(std::vector<Residue> vec, Residue newElement){
	reverse(vec.begin(),vec.end());
	vec.push_back(newElement);
	reverse(vec.begin(),vec.end());
	return vec;
}
//calc the sum of elements in the vector
double vecUtil::sum(const std::vector<double>& vec){
	double sum = 0;
	for (int i = 0; i < vec.size(); i++){
		sum += vec[i];
	}
	return sum;
}
//calc the average of the elements in the vector
std::vector<double> vecUtil::average(std::vector< std::vector<double> > vec){
	std::vector<double> result;
	for (int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for(int j = 0; j < vec.size();j++){
			sum += vec[j][i];
		}
		result.push_back(sum/vec.size());
	}
	return result;
}
//calculate average for every column in vector 
std::vector<double> vecUtil::average(std::vector< std::vector<int> > vec){
	std::vector<double> result;
	for (int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for(int j = 0; j < vec.size();j++){
			sum += vec[j][i];
		}
		result.push_back(sum/vec.size());
	}
	return result;
}
int vecUtil::countTrueValuesInVector(const std::vector<bool>& vec){
	int result = 0;
	for (int i = 0; i < vec.size(); i++){
		if (vec[i]){
			result++;	
		}
	}
	return result;
}
void vecUtil::printSequence(std::vector<Residue> vec){
	for (int i = 0; i < vec.size(); i++){
		std::cout << vec[i].getAA();
	}
	std::cout << std::endl;
}
//returns index of the first occurence of val in vec
int vecUtil::findIndex(std::string val, std::vector<std::string> vec){
	int res = -1;
	for (int i = 0; i < vec.size(); i++){
		if (vec[i] == val){
			res = i;
			break;
		}
	}
	return res;
}
