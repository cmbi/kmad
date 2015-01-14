#ifndef SBSTMATRIX_H
#define SBSTMATRIX_H


#include "fasta.h"

#include <iostream>
#include <vector>


namespace substitution_matrix
{
  ///
  /// returns a column from the sbst matrix for the aa amino acid
  ///
  std::vector<double> get_column(const char&);
}

#endif /* SBSTMARTIX_H */
