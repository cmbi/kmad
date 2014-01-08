#include<iostream>
#include<vector>
#include<string>

namespace vecUtil{
	bool contains(std::vector<std::string>&, std::string);
	void transposeVec(std::vector< std::vector<int> >&);
	void transposeVec(std::vector< std::vector<double> >&);
	void divideVectorByAScalar(std::vector<double>&, int);
	void divideVectorByAScalar(std::vector<double>&, double);
	void multiplyVectorByAScalar(std::vector<double>&, double);
	std::vector<double> addUp(std::vector< std::vector<double> >);
	double sum(const std::vector<double>&);
	std::vector<double> convertIntVectorToDoubleVector(std::vector<int>);
	void printDoubleVector(const std::vector<double>&);
	std::vector<std::string> flattenWithoutFeatures(const std::vector<std::vector<std::string> >&);
	std::vector<std::string> flatten(const std::vector<std::vector<std::string> >&);
	std::vector<std::string> push_front(std::vector<std::string>,std::string);
}
