#ifndef SBSTMATRIX_H
#define SBSTMATRIX_H


#include "fasta.h"
#include "types.h"

#include <iostream>
#include <vector>


class Residue;


namespace substitution_matrix
{
  ///
  /// Creates a profile matrix (20 x sequence length)
  /// based on a single sequence
  ///
  ProfileMatrix ConvertToProfileFormat(std::vector<fasta::Residue>& seq);
  ///
  /// gets score from the sbst matrix for the given amino acids
  /// @param char1 amino acid code
  /// @param char2 amino acid code
  ///
  int get_element(char char1, char char2);
  ///
  /// gets score from the sbst matrix for the given amino acids
  /// @param i amino acid index
  /// @param j amino acid index
  ///
  int get_element(int i, int j);
  ///
  /// assigns a column from the susbtitution matrix to the 'column_int'
  ///
  void get_column(unsigned int& columnNo, SbstMatColumn& column_int);
  ///
  /// returns a column from the sbst matrix for the aa amino acid
  ///
  std::vector<double> get_column(const char&);
  ///
  /// returns index of the given amino acid
  ///
  int FindAminoAcidsIndex(char aa);
}

#endif /* SBSTMARTIX_H */
