#ifndef SBSTMATRIX_H
#define SBSTMATRIX_H

#include <iostream>
#include <vector>
class Residue;
namespace substitutionMatrix{
	void printSbstMatrix();
	std::vector< std::vector<double> > convertToProfileFormat(std::vector<Residue>);
	//getters
	char getLetter(int);
	int getElement(char,char);
	int getElement(int,int);
	std::vector<int> getColumn(int);
	int findAminoAcidsNo(char);
}

#endif /* SBSTMARTIX_H */
