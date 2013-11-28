#include <iostream>
#include <vector>
namespace substitutionMatrix{
	void printSbstMatrix();
	std::vector< std::vector<double> > convertToProfileFormat(std::string);
	std::vector< std::vector<double> > convertToProfileFormat(std::vector<std::string>);
	//getters
	char getLetter(int);
	int getElement(char,char);
	int getElement(int,int);
	std::vector<int> getColumn(int);
	int findAminoAcidsNo(char);
};
