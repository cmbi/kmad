#ifndef TXTPROC_H
#define TXTPROC_H

#include "Sequences.h"

#include <iostream>
#include <string>
#include <vector>
class FeaturesProfile;
class Sequences;
namespace txtProc{
	int convertStringToInt(std::string&);
	double convertStringToDouble(std::string&);
	std::string charToString(char);
	std::string charToString(char,char);
	std::istream& safeGetline(std::istream&, std::string&);
	std::vector<std::string> split(const std::string&, char);
	bool acceptedChar(char);
	Sequences read_fasta(std::string, int, std::vector<std::string>*,
                       std::vector<double>*);
  //std::vector<std::vector<std::vector<std::string> > > read_fasta(std::string, int, std::vector<std::string>*,
  //                     std::vector<double>*);
	void writeAlignmentToFile(std::vector<std::string>&, 
                            std::vector<std::string>&,
                            std::string);
	void writeAlignmentWithoutCodeToFile(std::vector<std::string>&, 
                                       std::vector< std::vector< std::vector<std::string> > >&, 
                                       std::string);
	void writeVector(std::vector<std::vector<double> >&, std::string);
	void process_conf_file(std::string, FeaturesProfile&, Sequences&);
	std::vector<int> unfold(std::string, std::vector<std::string>&);
}

#endif /* TXTPROC_H */
