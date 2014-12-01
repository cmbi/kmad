#ifndef MISC_H
#define MISC_H

#include "types.h"
#include <iostream>
#include <vector>
#include <fstream>

class Residue;
namespace misc{
  ///
  /// creates a gap Residue object
  ///
	Residue gapRes(int codon_length);
  ///
  /// creates a polyA sequence of length seqLength, where each resdue is coded 
  /// by a codon of length codon_length
  ///
  sequence pseudoResidueSequence(int seqLength, int codon_length);
  ///
  /// looks for mistakes in the given command line arguments
  ///
	bool checkParameters(int codonLength, int phosph, int domain, int motif, 
                       double gep, double gop, double endPenalty);
  ///
  /// returns True if the file 'name' exists
  ///
  bool file_exists(const std::string* name);
}

#endif /* MISC_H */
