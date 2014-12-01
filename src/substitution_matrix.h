#ifndef SBSTMATRIX_H
#define SBSTMATRIX_H

#include "types.h"
#include <iostream>
#include <vector>
class Residue;
namespace substitution_matrix{
  ///
  /// Creates a profile matrix (20 x sequence length)
  /// based on a single sequence
  ///
	ProfileMatrix convertToProfileFormat(ResidueSequence& seq);
  ///
  /// gets score from the sbst matrix for the given amino acids
  /// @param char1 amino acid code
  /// @param char2 amino acid code
  ///
	int getElement(char char1, char char2);
  ///
  /// gets score from the sbst matrix for the given amino acids
  /// @param i amino acid index
  /// @param j amino acid index
  ///
	int getElement(int i, int j);
  /// 
  /// assigns a column from the susbtitution matrix to the 'column_int'
  ///
	void getColumn(unsigned int& columnNo, SbstMatColumn& column_int);
  ///
  /// returns index of the given amino acid
  ///
	int findAminoAcidsIndex(char aa);
}

#endif /* SBSTMARTIX_H */
