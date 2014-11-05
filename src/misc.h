#ifndef MISC_H
#define MISC_H

#include <iostream>
#include <vector>
class Residue;
namespace misc{
	Residue gapRes(int);
	std::vector<Residue> pseudoResidueSequence(int, int);
	void printEncodedSeq(const std::vector<std::string>&);
	bool checkParameters(int, int, int, int, double, double, bool, double);
}

#endif /* MISC_H */
