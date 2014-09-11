#include "vecUtil.h"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
//function contains
bool vecUtil::contains(std::vector<std::string>& vec, std::string x){
	if (std::find(vec.begin(),vec.end(),x) != vec.end()) return true;
	else return false;
}
//function transposeVec
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
	//return newVec;
}
//function transposeVec  - transposes vector< vector<double> > and returns vector< vector<double> >
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
void vecUtil::printDoubleVector(const std::vector<double>& vec){
	for (int i = 0; i < vec.size(); i++){
		std::cout << vec[i];
	}
	std::cout << "\n";
}
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
std::vector<std::string> vecUtil::push_front(std::vector<std::string> vec, std::string newElement){
	reverse(vec.begin(),vec.end());
	vec.push_back(newElement);
	reverse(vec.begin(),vec.end());
	return vec;
}
double vecUtil::sum(const std::vector<double>& vec){
	double sum = 0;
	for (int i = 0; i < vec.size(); i++){
		sum += vec[i];
	}
	return sum;
}
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
//claculate average for vector of vectors of doubles
double vecUtil::singleAverage(std::vector<std::vector<double> > vec){
	double sum = 0;
	for (int i = 0; i < vec.size();i++){
		for (int j = 0; j < vec[0].size();j++){
			sum += vec[i][j];
		}
	}
	sum = sum / (vec.size()*vec[0].size());
	return sum;
}
double vecUtil::max(std::vector<std::vector<double> > vec){
	double max = -1000;
	for (int i = 0; i < vec.size();i++){
		for (int j = 0; j < vec[0].size();j++){
			if (max < vec[i][j]){
				max = vec[i][j];
			}
		}
	}
	return max;
}
double vecUtil::min(std::vector<std::vector<double> > vec){
	double min = 1000000;
	for (int i = 0; i < vec.size();i++){
		for (int j = 0; j < vec[0].size();j++){
			if (min > vec[i][j]){
				min = vec[i][j];
			}
		}
	}
	return min;
}
double vecUtil::median(std::vector<std::vector<double> > vec){
	std::vector<double> sorted;
	for (int i = 0; i < vec.size();i++){
		for(int j = 0; j < vec[0].size();j++){
			sorted.push_back(vec[i][j]);
		}
	}	
	double tmp;
	for (int i = 0; i < sorted.size();i++){
		for(int j = i; j < sorted.size();j++){
			if (sorted[j]> sorted[j+1]){
				tmp = sorted[j];
				sorted[j] = sorted[j+1];
				sorted[j+1] = tmp;
			}	
		}
	}
	return sorted[int(sorted.size()/2)];
}
