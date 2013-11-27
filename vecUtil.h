#include<iostream>
#include<vector>
#include<string>

namespace vecUtil{
	void transposeVec(std::vector< std::vector<int> >&);
	void transposeVec(std::vector< std::vector<double> >&);
	void divideVectorByAScalar(std::vector<double>&, int);
	void multiplyVectorByAScalar(std::vector<double>&, double);
	std::vector<double> addUp(std::vector< std::vector<double> >);
	std::vector<double> convertIntVectorToDoubleVector(std::vector<int>);
	void printDoubleVector(const std::vector<double>&);
}
