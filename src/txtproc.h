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
  /// writes a line of text from the istream to string t
  /// (removes newline characters)
  ///
  std::istream& SafeGetline(std::istream& is, std::string& t);
  ///
  /// checks if character from the fasta file is on the list of accepetd
  /// characters
  ///
  bool AcceptedChar(char my_char);
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
