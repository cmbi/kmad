#include "vecUtil.h"
#include "Residue.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>


//checks if the vector of strings vec contains the string x
bool vecUtil::contains(std::vector<std::string>& vec, std::string& x){
	if (std::find(vec.begin(),vec.end(),x) != vec.end()) return true;
	else return false;
}


void vecUtil::transposeVec(std::vector< std::vector<int> >& vec){
	std::vector< std::vector<int> > newVec;
	std::vector<int> newRow;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		for (unsigned int j = 0; j < vec.size(); j++){
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
	for (unsigned int i = 0; i < vec[0].size(); i++){
		for (unsigned int j = 0; j < vec.size(); j++){
			newRow.push_back(vec[j][i]);	
		}
		newVec.push_back(newRow);
		newRow.clear();
	}
	vec = newVec;
}


void vecUtil::divideVectorByAScalar(std::vector<double>& vec, int scalar){
	std::vector<double> result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vecUtil::divideVectorByAScalar(std::vector<double>& vec, double& scalar){
	std::vector<double> result;
  for (auto &item: vec){
		result.push_back(item/scalar);	
	}
	vec = result;
}


void vecUtil::multiplyVectorByAScalar(std::vector<double>& vec, int scalar){
	std::vector<double> result;
  for (auto &item: vec){
		result.push_back(item*scalar);
	}
	vec = result;
}


void vecUtil::multiplyVectorByAScalar(std::vector<double>& vec, double& scalar){
	std::vector<double> result;
  for (auto &item: vec){
		result.push_back(item*scalar);
	}
	vec = result;
}


//function addUp - takes 2D matrix, adds up elements from each column, returns a 1D vector
std::vector<double> vecUtil::addUp(std::vector< std::vector<double> >& vec){
	std::vector<double> newVec;
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
std::vector<double> vecUtil::convertIntVectorToDoubleVector(std::vector<int>& vec){
	std::vector<double> result;
  for (auto &item: vec){
		result.push_back(double(item));
	}
	return result;
}


//function printDoubleVector
void vecUtil::printVector(const std::vector<int>& vec){
  for (auto &item: vec){
		std::cout << item << " ";
	}
	std::cout << "\n";
}


void vecUtil::printVector(const std::vector<double>& vec){
  for (auto &item: vec){
		std::cout << item << " ";
	}
	std::cout << "\n";
}


void vecUtil::printVector(const std::vector<std::string>& vec){
  for (auto &item: vec){
		std::cout << item << " ";
	}
	std::cout << "\n";
}


//takes a set of encoded sequences (=vector of vectors of strings) and returns a vector of nonencoded sequences
std::vector<std::string> vecUtil::flattenWithoutFeatures(const std::vector<std::vector<std::string> >& vec){
	std::vector<std::string> result;
  for (auto &row: vec){
		std::string newSeq = "";
    for (auto &item: row){
			newSeq += item[0];
		}
		result.push_back(newSeq);
	}
	return result;
}


//flatten a vector of residue vectors to a vector of strings
std::vector<std::string> vecUtil::flatten(const sequenceList& vec){
	std::vector<std::string> result;
  for (auto &row: vec){
		std::string newSeq = "";
    for(auto &item: row){
			newSeq += item.getCodon();
		}
		result.push_back(newSeq);
	}
	return result;
}


std::vector<std::string> vecUtil::flatten(const std::vector<std::vector<std::string> >& vec){
	std::vector<std::string> result;
  for (auto &row: vec){
		std::string newSeq = "";
    for (auto &item: row){
			newSeq += item;
		}
		result.push_back(newSeq);
	}
	return result;
}


//add an element to the beginning of the vector
std::vector<Residue> vecUtil::push_front(sequence& seq, Residue newElement){
	reverse(seq.begin(),seq.end());
	seq.push_back(newElement);
	reverse(seq.begin(),seq.end());
	return seq;
}


//calc the sum of elements in the vector
double vecUtil::sum(const std::vector<double>& vec){
	double sum = 0;
  for (auto &item: vec){
		sum += item;
	}
	return sum;
}


//calc the average of the elements in the vector
std::vector<double> vecUtil::average(std::vector< std::vector<double> >& vec){
	std::vector<double> result;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for(unsigned int j = 0; j < vec.size();j++){
			sum += vec[j][i];
		}
		result.push_back(sum/vec.size());
	}
	return result;
}


//calculate average for every column in vector 
std::vector<double> vecUtil::average(const std::vector< std::vector<int> >& vec){
	std::vector<double> result;
	for (unsigned int i = 0; i < vec[0].size(); i++){
		double sum = 0;
		for(unsigned int j = 0; j < vec.size();j++){
			sum += vec[j][i];
		}
		result.push_back(sum/vec.size());
	}
	return result;
}


int vecUtil::countTrueValuesInVector(const std::vector<bool>& vec){
	int result = 0;
  for (const auto &item: vec){
		if (item){
			result++;	
		}
	}
	return result;
}


void vecUtil::printSequence(sequence& seq){
  for (auto &res: seq){
		std::cout << res.getAA();
	}
	std::cout << std::endl;
}


//returns index of the first occurence of val in vec
int vecUtil::findIndex(std::string& val, std::vector<std::string>& vec){
	int res = -1;
	for (unsigned int i = 0; i < vec.size(); i++){
		if (vec[i] == val){
			res = i;
			break;
		}
	}
	return res;
}
