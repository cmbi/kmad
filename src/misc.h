#include <iostream>
#include <vector>
class Residue;
namespace misc{
//	std::vector<std::string> pseudoSequence(int,int);
	Residue gapRes(int);
	std::vector<Residue> pseudoResidueSequence(int,int);
	void printEncodedSeq(const std::vector<std::string>&);
	bool checkParameters(int, int, int, int, double, double,bool,double);
};
