#ifndef VECUTIL_H
#define VECUTIL_H 

#include<iostream>
#include<vector>
#include<string>
class Residue;
namespace vecUtil{
	bool contains(std::vector<std::string>&, std::string);
	int findIndex(std::string, std::vector<std::string>&);
	void transposeVec(std::vector< std::vector<int> >&);
	void transposeVec(std::vector< std::vector<double> >&);
	void divideVectorByAScalar(std::vector<double>&, double);
	void multiplyVectorByAScalar(std::vector<double>&, double);
	std::vector<double> addUp(std::vector< std::vector<double> >&);
	double sum(const std::vector<double>&);
	std::vector<double> convertIntVectorToDoubleVector(std::vector<int>&);
	void printVector(const std::vector<double>&);
	void printVector(const std::vector<int>&);
	void printVector(const std::vector<std::string>&);
	std::vector<std::string> flattenWithoutFeatures(const std::vector<std::vector<std::string> >&);
	std::vector<std::string> flatten(const std::vector<std::vector<std::string> >&);
	std::vector<std::string> flatten(const std::vector<std::vector<Residue> >&);
	std::vector<Residue> push_front(std::vector<Residue>&, Residue);
	std::vector<double> average(std::vector<std::vector<double> >&);
	//std::vector<double> average(std::vector<std::vector<int> >&);
	std::vector<double> average(const std::vector<std::vector<int> >&);
	int countTrueValuesInVector(const std::vector<bool>&);
	void printSequence(std::vector<Residue>&);
}

#endif /* VECUTIL_H */
