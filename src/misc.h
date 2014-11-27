#ifndef MISC_H
#define MISC_H

#include "types.h"
#include <iostream>
#include <vector>
#include <fstream>

class Residue;
namespace misc{
	Residue gapRes(int);
  std::vector<Residue> pseudoResidueSequence(int, int);
	bool checkParameters(int, int, int, int, double, double, bool, double);
  bool file_exists(const std::string*);
}

#endif /* MISC_H */
