#ifndef SBSTMATRIX_H
#define SBSTMATRIX_H

#include <iostream>
#include <vector>
class Residue;
namespace substitutionMatrix{
	std::vector< std::vector<double> > convertToProfileFormat(std::vector<Residue>&);
	//getters
	int getElement(char, char);
	int getElement(int, int);
	void getColumn(unsigned int&, std::vector<int>&);
	int findAminoAcidsNo(char);
}

#endif /* SBSTMARTIX_H */
