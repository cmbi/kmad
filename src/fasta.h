#ifndef FASTA_H
#define FASTA_H

#include "sequences.h"
#include "types.h"

#include <iostream>
#include <fstream>
#include <vector>

class Sequences;

namespace fasta{
  ///
  /// parsess a 'fasta' file, returns a Sequences object, 
  /// writes motif ids and probabilities to ids and probs
  ///
  Sequences parse_fasta(std::string filename, int codonLength, 
                        IDsList* ids, ProbsList* probs);
  
}


#endif /* FASTA_H */
