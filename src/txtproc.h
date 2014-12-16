#ifndef TXTPROC_H
#define TXTPROC_H


#include <iostream>
#include <string>
#include <vector>


class FeaturesProfile;


namespace txtproc{
  ///
  /// writes alignment to file (encoded)
  ///
  void WriteAlignmentToFile(std::vector<std::string>& sequences,
                            std::vector<std::string>& sequence_names,
                            std::string filename);
  ///
  /// writes alignment to a file in the regular fasta format (decoded)
  ///
  void WriteAlignmentWithoutCodeToFile(std::vector<std::string>& sequences,
                                       std::vector<std::string>& sequence_names,
                                       std::string filename, int codon_length);
}

#endif /* TXTPROC_H */
