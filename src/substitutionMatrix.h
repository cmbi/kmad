#ifndef SBSTMATRIX_H
#define SBSTMATRIX_H

#include "types.h"
#include <iostream>
#include <vector>
class Residue;
namespace substitutionMatrix{
	profile_matrix convertToProfileFormat(sequence& seq);
	//getters
	int getElement(char char1, char char2);
	int getElement(int i, int j);
	void getColumn(unsigned int& columnNo, sbstMatColumn& column_int);
	int findAminoAcidsNo(char aa);
}

#endif /* SBSTMARTIX_H */
