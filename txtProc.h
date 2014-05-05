#include <iostream>
#include <string>
#include <vector>
namespace txtProc{
	int convertStringToInt(std::string);
	double convertStringToDouble(std::string);
	std::string charToString(char);
	std::string gapCode(int);
	std::istream& safeGetline(std::istream&, std::string&);
	std::vector<std::string> &split(const std::string &, char, std::vector<std::string> &);
	std::vector<std::string> split(const std::string &, char);
	std::vector< std::vector<std::string> > processFASTA(std::string);
	bool acceptableChar(char);
	std::vector< std::vector< std::vector<std::string> > >processFASTA(std::string,int, std::vector<std::string>*, std::vector<double>*);
	void writeAlignmentToFile(std::vector<std::string>,std::vector< std::vector<std::string> >, std::string);
	void writeAlignmentToFile(std::vector<std::string>,std::vector< std::vector< std::vector<std::string> > >, std::string);
	void writeAlignmentWithoutCodeToFile(std::vector<std::string>,std::vector< std::vector< std::vector<std::string> > >, std::string);
	void writeVector(std::vector<std::vector<double> >,std::string);
}
