#include <iostream>
#include <string>
#include <vector>
namespace txtProc{
	int convertStringToInt(std::string);
	std::vector< std::vector<std::string> > processFASTA(std::string);
	std::vector< std::vector< std::vector<std::string> > >processFASTA(std::string,int);
	void writeAlignmentToFile(std::vector<std::string>,std::vector< std::vector<std::string> >, std::string);
	void writeAlignmentToFile(std::vector<std::string>,std::vector< std::vector< std::vector<std::string> > >, std::string);
}
