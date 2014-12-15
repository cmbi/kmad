#ifndef TXTPROC_H
#define TXTPROC_H

#include "sequences.h"

#include <iostream>
#include <string>
#include <vector>


class FeaturesProfile;
class Sequences;
namespace txtproc{
  ///
  /// writes alignment to file (encoded)
  ///
  void WriteAlignmentToFile(StringSequences& sequences,
                            SeqNames& sequence_names,
                            std::string filename);
  ///
  /// writes alignment to a file in the regular fasta format (decoded)
  ///
  void WriteAlignmentWithoutCodeToFile(StringSequences& sequences,
                                       SeqNames& sequence_names,
                                       std::string filename, int codon_length);
}

#endif /* TXTPROC_H */
